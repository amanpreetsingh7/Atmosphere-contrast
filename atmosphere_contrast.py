"""Calculate thermal planet/star contrast from atmosphere model spectra.

The public workflow is load_spectrum -> Response/load_nprime ->
contrast_from_spectra. Use contrast_from_models for direct file
mappings,
or planet_contrast when wavelength/flux arrays are already available.

All internal spectra are emergent surface F_lambda in W m^-2 um^-1.
Photon contrast is (Rp/Rs)^2 times the ratio of integrals of
F_lambda * filter * sky * wavelength / (h*c). The integral coordinate is
micrometres; only the photon-energy factor uses wavelength in metres.

This module includes no app, plots, bundled predictions, catalog search,
model interpolation, stellar heating, PSF generation, or source
injection.
Native model files and transmission curves are external inputs. Missing
coverage, incompatible units, and unsupported formats raise exceptions.

NumPy is required. Astropy and netCDF4 are optional imports at module
scope, required only by their respective native readers. See README.md
for supported formats, physical assumptions, and reproducible examples.
"""

__version__ = "0.1.0"
import gzip
import hashlib
import io
import json
import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

try:
    from astropy import units as u
    from astropy.io import fits
    from astropy.table import Table
except ImportError:
    u = None
    fits = None
    Table = None
try:
    from netCDF4 import Dataset
except ImportError:
    Dataset = None
H = 6.62607015e-34
PC_M = 3.085677581491367e16
_NUMBER = "([+-]?\\d+(?:\\.\\d+)?)"
MAX_TEMPERATURE_OFFSET_K = 1.0
MAX_GRAVITY_LABEL_OFFSET_M_S2 = 1.0
C = 299792458.0
FLUX_UNIT = "W m^-2 um^-1"
WAVE_FACTORS = {
    "um": 1.0,
    "nm": 0.001,
    "angstrom": 0.0001,
    "m": 1000000.0,
    "cm": 10000.0,
}
PYNOMIC_COMMIT = "e6c63c740a6fe1ca74e57228fafa3d7da8b6e67a"


@dataclass
class Spectrum:
    """Store an emergent surface spectrum and source metadata.

    Parameters
    ----------
    wavelength_um : array-like
        Positive wavelengths in micrometres; sorted on construction.
    surface_flux : array-like
        Matching nonnegative F_lambda in W m^-2 um^-1. This is not
        an
        observed flux at Earth and not a specific intensity.
    metadata : dict, optional
        Source files/hashes, parameters, and air/vacuum convention.

    Raises
    ------
    ValueError
        For mismatched arrays, nonfinite values, or duplicate
        wavelengths.
    """

    wavelength_um: np.ndarray
    surface_flux: np.ndarray
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        """Validate and sort spectral arrays; reject duplicate wavelengths."""
        x, y = (
            np.asarray(self.wavelength_um, float),
            np.asarray(self.surface_flux, float),
        )
        if x.ndim != 1 or y.shape != x.shape or len(x) < 3:
            raise ValueError(
                (
                    "A spectrum needs matching "
                    "one-dimensional arrays with "
                    "at least three samples."
                )
            )
        if (
            not np.all(np.isfinite(x))
            or np.any(x <= 0)
            or (not np.all(np.isfinite(y)))
            or np.any(y < 0)
        ):
            raise ValueError(
                (
                    "Wavelengths must be positive "
                    "and finite; flux must be "
                    "finite and nonnegative."
                )
            )
        # Sort coordinates and flux together, then reject
        # ambiguous duplicates.
        order = np.argsort(x)
        x, y = (x[order], y[order])
        if np.any(np.diff(x) <= 0):
            raise ValueError(
                (
                    "Duplicate wavelengths are "
                    "ambiguous. Resolve them in "
                    "the source file."
                )
            )
        self.wavelength_um, self.surface_flux = (x, y)


@dataclass
class Response:
    """Represent separate filter and sky transmission grids.

    Parameters
    ----------
    filter_wave_um, sky_wave_um : array-like
        Strictly increasing wavelengths in micrometres.
    filter_transmission, sky_transmission : array-like
        Matching dimensionless transmissions in [0,1]. Neither is an
        already photon-weighted response or a sky-emission spectrum.
    limits_um : tuple of float
        Integration lower/upper wavelengths. Both grids must cover
        them.
    metadata : dict
        Response provenance and wavelength_frame (vacuum or air).
        The
        frame is an adopted convention when source files do not
        state it.

    Raises
    ------
    ValueError
        If components are invalid or do not cover the requested
        interval.

    Notes
    -----
    The response is the product of independently interpolated filter and
    sky curves. Missing telescope optics/QE and observation-specific sky
    conditions are not inferred. A shared response scale cancels in
    contrast.
    """

    filter_wave_um: np.ndarray
    filter_transmission: np.ndarray
    sky_wave_um: np.ndarray
    sky_transmission: np.ndarray
    limits_um: tuple
    metadata: dict

    def __post_init__(self):
        """Validate both transmission grids and require full band coverage."""
        self.filter_wave_um, self.filter_transmission = _validate_transmission(
            self.filter_wave_um, self.filter_transmission, "Filter"
        )
        self.sky_wave_um, self.sky_transmission = _validate_transmission(
            self.sky_wave_um, self.sky_transmission, "Sky"
        )
        lo, hi = map(float, self.limits_um)
        if not np.isfinite([lo, hi]).all() or lo >= hi:
            raise ValueError("Response limits must be finite and increasing.")
        # Both response components must independently cover
        # the entire band.
        for x in (self.filter_wave_um, self.sky_wave_um):
            if lo < x[0] or hi > x[-1]:
                raise ValueError(
                    (
                        "The filter and sky must both "
                        "cover the entire chosen band; "
                        "extrapolation is disabled."
                    )
                )
        self.limits_um = (lo, hi)

    def at(self, wave):
        """Return filter times sky transmission at query wavelengths in um.

        wave may be an array. Return dimensionless values interpolated
        linearly on each separate original grid. Raise ValueError
        outside
        the declared band; no extrapolation or partial-band clipping.
        """
        q = np.asarray(wave, float)
        lo, hi = self.limits_um
        if not np.isfinite(q).all() or np.any(q < lo) or np.any(q > hi):
            raise ValueError("Response query outside the declared band.")
        return np.interp(
            q, self.filter_wave_um, self.filter_transmission
        ) * np.interp(q, self.sky_wave_um, self.sky_transmission)


def to_surface_flambda(
    wavelength, flux, wavelength_unit, flux_unit, convention="surface_flux"
):
    """Convert declared spectral units to um and W m^-2 um^-1.

    Parameters
    ----------
    wavelength, flux : array-like
        Input spectral samples. Spectrum validates their final
        shapes.
    wavelength_unit : {"um", "nm", "angstrom", "m", "cm"}
        Wavelength unit token; no air/vacuum conversion is
        performed.
    flux_unit : str
        One of W_m2_um, W_m2_m, erg_cm2_s_um, erg_cm2_s_A,
        erg_cm2_s_cm, erg_cm2_s_Hz, or W_m2_Hz.
    convention : {"surface_flux", "eddington_flux"}, optional
        Eddington H_lambda is converted to emergent flux by 4*pi.

    Returns
    -------
    tuple of numpy.ndarray
        Wavelength in um and emergent surface F_lambda in W m^-2
        um^-1.

    Raises
    ------
    ValueError
        If unit/convention is unsupported, or F_nu wavelengths
        nonpositive.

    Notes
    -----
    F_lambda(per um) = F_nu * c/lambda_m^2 * 1e-6. A cgs F_nu
    first needs a 1e-3 conversion to W m^-2 Hz^-1. Observer fluxes and
    specific intensities require different treatment and are not
    accepted
    as conventions. Surface flux needs no additional pi multiplier.
    """
    if wavelength_unit not in WAVE_FACTORS:
        raise ValueError(
            "Supported wavelength units: " + ", ".join(WAVE_FACTORS)
        )
    # Convert the coordinate first; flux-density conversion uses
    # its SI value.
    wave = np.asarray(wavelength, float) * WAVE_FACTORS[wavelength_unit]
    flux = np.asarray(flux, float)
    # These factors include energy/area units and the per-
    # wavelength Jacobian.
    factors = {
        "W_m2_um": 1.0,
        "W_m2_m": 1e-06,
        "erg_cm2_s_um": 0.001,
        "erg_cm2_s_A": 10.0,
        "erg_cm2_s_cm": 1e-07,
    }
    if flux_unit in factors:
        converted = flux * factors[flux_unit]
    elif flux_unit in ("erg_cm2_s_Hz", "W_m2_Hz"):
        if np.any(wave <= 0):
            raise ValueError(
                ("Positive wavelengths are " "required for F_nu conversion.")
            )
        converted = (
            flux
            * (0.001 if flux_unit == "erg_cm2_s_Hz" else 1.0)
            * C
            / (wave * 1e-06) ** 2
            * 1e-06
        )
    else:
        raise ValueError("Unknown flux unit: " + flux_unit)
    # Only Eddington H_lambda needs 4*pi. Emergent surface flux
    # needs none.
    if convention == "eddington_flux":
        converted = converted * (4 * math.pi)
    elif convention != "surface_flux":
        raise ValueError(
            (
                "Use surface_flux or "
                "eddington_flux. Observed flux "
                "and intensity require "
                "different treatment."
            )
        )
    return (wave, converted)


def load_spectrum(spec, base):
    """Load a native or explicitly mapped ASCII surface spectrum.

    Parameters
    ----------
    spec : dict
        Required: file and format. Optional native settings include
        parameters, label, source_url, wavelength_frame,
        gravity_column,
        and wavelength_file (required for PHOENIX high resolution).
        Format ascii additionally requires wavelength_unit,
        flux_unit,
        and flux_convention; skiprows, delimiter, wavelength_column,
        and
        flux_column may select an exported table explicitly.
    base : str or pathlib.Path
        Parent directory used to resolve relative file paths.

    Returns
    -------
    Spectrum
        Surface F_lambda in W m^-2 um^-1 and a source metadata
        dictionary.
        Header-verified parameters are kept distinct from
        declarations.

    Raises
    ------
    FileNotFoundError
        If a selected model or separate wavelength file is missing.
    ImportError
        If the selected native reader needs astropy or netCDF4.
    ValueError
        If native layout, units, coordinates, or conventions
        disagree.

    Notes
    -----
    No downloads, catalog lookup, interpolation, blackbody substitution,
    or sigma*T^4 renormalization occur. Supported format tokens are
    listed
    in README.md. Each adapter applies the original release's
    conventions;
    converted data require an explicit ASCII mapping.
    """
    # Resolve only the selected file. Never replace a missing
    # atmosphere.
    path = (Path(base) / spec["file"]).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            (
                "Model file missing: %s. "
                "Install or explicitly select "
                "the required model; no "
                "substitute is inserted."
            )
            % path
        )
    # Dispatch by an explicit release token, not by guessing from
    # the file suffix.
    fmt = spec["format"]
    if fmt in {
        "ck04_stsci_fits",
        "phoenix_hires",
        "atmo2020_stsci_fits",
        "bt_settl_svo",
        "sonora_elf_owl_v2",
        "sonora_flame_skimmer_v1_r3k",
    }:
        _require_optional("astropy", u)
    # Keep declarations separate from values actually verified in
    # a header.
    metadata = {
        "label": spec.get("label", path.name),
        "format": fmt,
        "parameters_declared": spec.get("parameters", {}),
        ("parameters_verified_from_heade" "r"): {},
        "files": [_file_record(path)],
        "wavelength_frame": spec.get("wavelength_frame", "vacuum"),
        "source_url": spec.get("source_url"),
        "flux_output_unit": FLUX_UNIT,
    }
    if metadata["wavelength_frame"] not in ("air", "vacuum"):
        raise ValueError("Declare wavelength_frame as air or vacuum.")
    # ASCII exports need declared units; native Bobcat units come
    # from its release.
    if (
        fmt
        in (
            "phoenix_hires",
            "sonora_bobcat",
            "sonora_cholla",
            "sonora_diamondback",
            "atmo2020_ascii",
            "atmo2020_stsci_fits",
            "sonora_elf_owl_v2",
            "sonora_flame_skimmer_v1_r3k",
        )
        and metadata["wavelength_frame"] != "vacuum"
    ):
        raise ValueError(
            (
                "This native emission release "
                "uses vacuum wavelengths. A "
                "converted export needs an "
                "explicit ASCII mapping, not a "
                "changed native-format label."
            )
        )
    # ASCII exports need declared units; native Bobcat units come
    # from its release.
    if fmt in ("ascii", "sonora_bobcat"):
        text = (
            gzip.open(path, "rt").read()
            if path.suffix.lower() == ".gz"
            else path.read_text()
        )
        text = text.replace("D+", "E+").replace("D-", "E-")
        if fmt == "sonora_bobcat":
            lines = text.splitlines()
            skip = None
            for i, line in enumerate(lines[:30]):
                tokens = line.split()
                if len(tokens) == 2:
                    try:
                        a, b = map(float, tokens)
                    except ValueError:
                        continue
                    if a > 0 and b >= 0:
                        skip = i
                        break
            if skip is None:
                raise ValueError(
                    (
                        "No two-column Bobcat spectrum "
                        "found after its header."
                    )
                )
        else:
            skip = spec.get("skiprows", 0)
        data = np.loadtxt(
            io.StringIO(text),
            skiprows=skip,
            delimiter=spec.get("delimiter"),
            usecols=(
                spec.get("wavelength_column", 0),
                spec.get("flux_column", 1),
            ),
            ndmin=2,
        )
        wave, flux = (data[:, 0], data[:, 1])
        if fmt == "sonora_bobcat":
            wave_unit, flux_unit, convention = (
                "um",
                "erg_cm2_s_Hz",
                "surface_flux",
            )
            metadata["original_header"] = text.splitlines()[:skip]
            metadata["header_lines_skipped"] = skip
            unit_header = "".join(
                " ".join(metadata["original_header"][1:]).lower().split()
            ).replace("^", "")
            if unit_header and (
                not ("micron" in unit_header and "erg/cm2/s/hz" in unit_header)
            ):
                raise ValueError(
                    (
                        "Bobcat native unit header "
                        "must specify microns and "
                        "erg/cm2/s/Hz. Use an explicit "
                        "ASCII mapping for a converted "
                        "file."
                    )
                )
            metadata["native_units_verified_from_header"] = bool(unit_header)
            if not unit_header:
                metadata["unit_note"] = (
                    "No unit line in this legacy "
                    "export; units follow the "
                    "explicitly selected Bobcat "
                    "release convention."
                )
            if skip >= 1:
                try:
                    vals = [float(v) for v in text.splitlines()[0].split()[:8]]
                except ValueError:
                    vals = []
                if len(vals) >= 8 and vals[0] > 0 and (vals[1] > 0):
                    measured = {
                        "teff_K": vals[0],
                        "logg_cgs": math.log10(vals[1] * 100),
                        "metallicity_dex": vals[5],
                        "co_relative_solar": vals[6],
                    }
                    _check_parameters(spec.get("parameters", {}), measured)
                    metadata["parameters_verified_from_header"] = measured
                    coordinates = _bobcat_coordinates(path.name, measured)
                    if coordinates is not None:
                        metadata["grid_coordinates"] = coordinates
        else:
            wave_unit, flux_unit = (spec["wavelength_unit"], spec["flux_unit"])
            convention = spec["flux_convention"]
    elif fmt == "ck04_stsci_fits":
        if metadata["wavelength_frame"] != "vacuum":
            raise ValueError(
                (
                    "This integration adopts "
                    "vacuum CK04 coordinates; "
                    "converted files need an "
                    "explicit adapter."
                )
            )
        wave, flux, extra = _read_ck04(path, spec)
        _check_parameters(
            spec.get("parameters", {}),
            extra[("parameters_verified_from_heade" "r")],
        )
        metadata.update(extra)
        wave_unit, flux_unit, convention = (
            "angstrom",
            "erg_cm2_s_A",
            "surface_flux",
        )
    elif fmt == "phoenix_hires":
        wave_path = (
            (Path(base) / spec["wavelength_file"]).expanduser().resolve()
        )
        with fits.open(path, mode="readonly", memmap=False) as hdus:
            flux = np.asarray(hdus[0].data, float)
            header = hdus[0].header
            measured = {
                k: float(header[h])
                for k, h in (
                    ("teff_K", "PHXTEFF"),
                    ("logg_cgs", "PHXLOGG"),
                    ("metallicity_dex", "PHXM_H"),
                    ("alpha_dex", "PHXALPHA"),
                )
                if h in header
            }
        with fits.open(wave_path, mode="readonly", memmap=False) as hdus:
            wave = np.asarray(hdus[0].data, float)
            wave_header = hdus[0].header
        flux_unit_header, wave_unit_header = (
            header.get("BUNIT"),
            wave_header.get("UNIT"),
        )
        if not flux_unit_header or not wave_unit_header:
            raise ValueError(
                (
                    "Native PHOENIX high-resolution"
                    " files need flux BUNIT and "
                    "wavelength UNIT headers. Use "
                    "an explicit ASCII mapping for "
                    "converted data."
                )
            )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*multiple slashes.*",
                category=u.UnitsWarning,
            )
            flux_factor = u.Unit(flux_unit_header).to(
                u.erg / u.s / u.cm**2 / u.cm
            )
            wave_factor = u.Unit(wave_unit_header).to(u.AA)
        if not (
            math.isclose(flux_factor, 1.0, rel_tol=1e-12)
            and math.isclose(wave_factor, 1.0, rel_tol=1e-12)
        ):
            raise ValueError(
                (
                    "PHOENIX high-resolution units "
                    "must be erg/s/cm2/cm and "
                    "Angstrom. This file has a "
                    "different convention."
                )
            )
        metadata["unit_headers_checked"] = {
            "BUNIT": flux_unit_header,
            "wavelength_UNIT": wave_unit_header,
        }
        metadata["files"].append(_file_record(wave_path))
        _check_parameters(spec.get("parameters", {}), measured)
        metadata["parameters_verified_from_header"] = measured
        wave_unit, flux_unit, convention = (
            "angstrom",
            "erg_cm2_s_cm",
            "surface_flux",
        )
    elif fmt in (
        "sonora_cholla",
        "sonora_diamondback",
        "atmo2020_ascii",
        "atmo2020_stsci_fits",
        "bt_settl_svo",
        "sonora_elf_owl_v2",
        "sonora_flame_skimmer_v1_r3k",
    ):
        wave, flux, wave_unit, flux_unit, extra = _load_native_format(
            path, fmt, spec
        )
        _check_parameters(
            extra[("parameters_verified_from_filen" "ame")],
            extra[("parameters_verified_from_heade" "r")],
        )
        _check_parameters(
            spec.get("parameters", {}),
            extra[("parameters_verified_from_filen" "ame")],
        )
        _check_parameters(
            spec.get("parameters", {}),
            extra[("parameters_verified_from_heade" "r")],
        )
        metadata.update(extra)
        convention = "surface_flux"
    else:
        raise ValueError(
            (
                "Unknown format %s. Use an "
                "explicit ASCII adapter for "
                "other exported grids."
            )
            % fmt
        )
    metadata.update(
        input_wavelength_unit=wave_unit,
        input_flux_unit=flux_unit,
        input_flux_convention=convention,
    )
    wave, flux = to_surface_flambda(
        wave, flux, wave_unit, flux_unit, convention
    )
    # All reader branches converge to the same surface-flux
    # representation.
    result = Spectrum(wave, flux, metadata)
    result.metadata["wavelength_range_um"] = [
        float(result.wavelength_um[0]),
        float(result.wavelength_um[-1]),
    ]
    result.metadata["native_samples"] = len(result.wavelength_um)
    return result


def load_nprime(directory, cutoff=None):
    """Load the reviewed pyNOMIC Nprime filter and reference Mt Graham sky.

    Parameters
    ----------
    directory : str or pathlib.Path
        Contains Nprime_filter_curve.txt and
        MtGrahamTransmittance.npz
        from PYNOMIC_COMMIT. Files are external, not downloaded
        here.
    cutoff : float or None, optional
        Explicit positive filter threshold below 1. None retains the
        full published filter interval (9-14 um in the reviewed
        table).

    Returns
    -------
    Response
        Separately tabulated filter/sky curves, limits, and SHA
        provenance.

    Raises
    ------
    OSError
        If either required external response file cannot be read.
    ValueError
        If files differ from reviewed hashes, normalization is
        unexpected,
        or cutoff leaves insufficient support.

    Notes
    -----
    The reviewed tiny negative filter artifacts are floored to zero;
    values below -0.001 are rejected and positive values not
    renormalized.
    Filter normalization is relative, not absolute throughput. The
    adopted
    vacuum frame and reference sky applicability are recorded, not
    measured.
    """
    directory = Path(directory)
    # Pin the reviewed response bytes so a changed table cannot go
    # unnoticed.
    manifest = {
        "repository": ("https://github.com/vivekv7323/" "pyNOMIC"),
        "commit": ("e6c63c740a6fe1ca74e57228fafa3d" "7da8b6e67a"),
        "sha256": {
            "Nprime_filter_curve.txt": (
                "8039385921be84c56bcbaf82c4f717"
                "4d3cd0dd0f323dfa94eb187ab8e7b6"
                "c9d6"
            ),
            "MtGrahamTransmittance.npz": (
                "36c83823d5abe2267e0b07b3c7b433"
                "6a4405c6acc4ffeb9b05d391240a1a"
                "eda2"
            ),
        },
        ("reference_function_source_sha2" "56"): (
            "cfa31331ac8454cf2239f605a19046"
            "ac6e41560e7221dc078c257df6b1a0"
            "197b"
        ),
        "license": (
            "MIT; Copyright (c) 2026 Vivek " "Vijayakumar; LICENSE included."
        ),
        "source_prefix": (
            "https://raw.githubusercontent."
            "com/vivekv7323/pyNOMIC/e6c63c7"
            "40a6fe1ca74e57228fafa3d7da8b6e"
            "67a/filterdata/"
        ),
    }
    if manifest["commit"] != PYNOMIC_COMMIT:
        raise ValueError("Unexpected pyNOMIC response release.")
    # Verify both source hashes before accepting the external
    # response curves.
    for name, expected in manifest["sha256"].items():
        if _sha256_file(directory / name) != expected:
            raise ValueError("Pinned response file has changed: " + name)
    data = np.loadtxt(directory / "Nprime_filter_curve.txt", skiprows=1)
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError("Nprime must have two native columns.")
    wave, raw = (data[:, 0], data[:, 1])
    if (
        not np.isfinite(raw).all()
        or np.min(raw) < -0.001
        or np.max(raw) > 1
        or (not math.isclose(float(np.max(raw)), 1.0, abs_tol=1e-10))
    ):
        raise ValueError(
            (
                "Nprime normalization or "
                "negative values differ from "
                "the reviewed table."
            )
        )
    # Clip only the reviewed negative artifacts; keep positive
    # values intact.
    transmission = np.maximum(raw, 0.0)
    with np.load(
        directory / "MtGrahamTransmittance.npz", allow_pickle=False
    ) as sky:
        sx, st = (np.array(sky["arr_0"], float), np.array(sky["arr_1"], float))
    lo, hi = (float(wave[0]), float(wave[-1]))
    if cutoff is not None:
        cutoff = float(cutoff)
        if not 0 < cutoff < 1:
            raise ValueError(
                ("An optional filter cutoff " "must be between zero and one.")
            )
        support = wave[transmission > cutoff]
        if len(support) < 2:
            raise ValueError("Cutoff leaves insufficient response support.")
        lo, hi = (float(support[0]), float(support[-1]))
    # Record adopted conventions and reference-sky limitations
    # with the response.
    metadata = {
        "name": ("pyNOMIC Nprime + supplied Mt " "Graham sky"),
        "repository": ("https://github.com/vivekv7323/" "pyNOMIC"),
        "commit": PYNOMIC_COMMIT,
        "files_sha256": manifest["sha256"],
        "filter_normalization": (
            "Published peak=1 shape; " "absolute throughput not " "supplied."
        ),
        ("negative_filter_samples_floore" "d"): int(np.count_nonzero(raw < 0)),
        "original_minimum_filter_value": float(np.min(raw)),
        "cleanup": (
            "Floor negative transmission "
            "samples to zero; reject "
            "values below -0.001; no "
            "positive-value renormalization"
            "."
        ),
        "cutoff": cutoff,
        "limits_um": [lo, hi],
        "default_detector": "photon",
        "wavelength_frame": "vacuum",
        "wavelength_frame_note": (
            "Shared vacuum convention "
            "adopted; source response "
            "files do not independently "
            "state air/vacuum."
        ),
        "sky_conditions": (
            "Provided reference table; "
            "PWV/airmass and observation-sp"
            "ecific applicability "
            "unverified."
        ),
        "calibration_scope": (
            "Filter shape times reference "
            "atmospheric transmission; no "
            "telescope area, optics, QE or "
            "background emission model."
        ),
    }
    return Response(wave, transmission, sx, st, (lo, hi), metadata)


def integrate_band(spectrum, response, detector="photon", rule="gauss3"):
    """Integrate one surface spectrum through filter times sky
    transmission.

    Parameters
    ----------
    spectrum : Spectrum
        Surface F_lambda in W m^-2 um^-1 on a complete wavelength
        grid.
    response : Response
        Separate transmission grids and explicit integration limits.
    detector : {"photon", "energy"}, optional
        Include wavelength/(h*c) only for photon detection.
    rule : {"gauss3", "boole"}, optional
        Three-point Gauss integration, or independent five-point
        Boole
        integration on every interval of the combined knot sequence.

    Returns
    -------
    float
        Response-weighted surface signal: photons m^-2 s^-1 for
        photon
        detection or W m^-2 for energy detection. It is not a
        calibrated
        telescope count rate or an observer flux without radius
        scaling.

    Raises
    ------
    ValueError
        For unsupported detector/rule, incomplete coverage, or frame
        mismatch.

    Notes
    -----
    Piecewise-linear F_lambda*filter*sky is cubic on each merged
    interval;
    multiplying by wavelength for photons makes it quartic. Both rules
    integrate this representation exactly up to floating-point rounding.
    That numerical property does not validate an atmosphere's physics.
    """
    if detector not in ("energy", "photon"):
        raise ValueError("Choose photon or energy detection.")
    # Preserve every spectrum/filter/sky breakpoint before integrating.
    nodes = _integration_nodes(spectrum, response)
    # No slope changes remain within any of these wavelength intervals.
    left, right = (nodes[:-1], nodes[1:])
    # Three-point Gauss is exact for the piecewise quartic photon
    # integrand.
    if rule == "gauss3":
        mid, half = ((left + right) / 2, (right - left) / 2)
        a = math.sqrt(3 / 5)
        return float(
            np.sum(
                half
                * (
                    5
                    / 9
                    * _band_integrand(
                        spectrum, response, mid - a * half, detector
                    )
                    + 8
                    / 9
                    * _band_integrand(spectrum, response, mid, detector)
                    + 5
                    / 9
                    * _band_integrand(
                        spectrum, response, mid + a * half, detector
                    )
                )
            )
        )
    if rule == "boole":
        step = (right - left) / 4
        return float(
            np.sum(
                (right - left)
                / 90
                * (
                    7 * _band_integrand(spectrum, response, left, detector)
                    + 32
                    * _band_integrand(
                        spectrum, response, left + step, detector
                    )
                    + 12
                    * _band_integrand(
                        spectrum, response, left + 2 * step, detector
                    )
                    + 32
                    * _band_integrand(
                        spectrum, response, left + 3 * step, detector
                    )
                    + 7 * _band_integrand(spectrum, response, right, detector)
                )
            )
        )
    raise ValueError("Unknown quadrature rule.")


def contrast_from_spectra(
    star, planet, star_radius_m, planet_radius_m, response, detector="photon"
):
    """Return the dimensionless thermal companion/star band contrast.

    Parameters
    ----------
    star, planet : Spectrum
        Emergent surface spectra in W m^-2 um^-1, on wavelength
        grids in um.
        Each must cover the complete response and use its wavelength
        frame.
    star_radius_m, planet_radius_m : float
        Positive emitting radii in metres. They are imposed
        independently
        of the atmosphere's temperature, gravity, mass, or
        evolutionary age.
    response : Response
        Filter and sky transmission, their grids, limits, and
        provenance.
    detector : {"photon", "energy"}, optional
        Photon weighting multiplies by wavelength/(h*c). Energy
        weighting
        integrates transmitted F_lambda without that factor.

    Returns
    -------
    float
        Companion/star band flux ratio, not percent or delta
        magnitude.
        A common observer distance cancels. This is not a detection
        limit.

    Raises
    ------
    ValueError
        If radii/signals are nonpositive or nonfinite, wavelength
        coverage
        is incomplete, the response is invalid, or frames disagree.

    Notes
    -----
    No extra pi factor applies to emergent surface flux. Telescope area
    and any common response normalization cancel in the ratio. Separate
    sky/filter grids are retained during integration; no model is fitted
    or renormalized. See integrate_band for the quadrature convention.
    """
    # Integrate each object independently through exactly the same
    # response.
    star_signal = integrate_band(star, response, detector)
    planet_signal = integrate_band(planet, response, detector)
    if not all(
        (
            math.isfinite(value) and value > 0
            for value in (
                star_radius_m,
                planet_radius_m,
                star_signal,
                planet_signal,
            )
        )
    ):
        raise ValueError(
            ("Radii and both band signals " "must be positive and finite.")
        )
    # Apply emitting-area scaling. The common observer distance cancels.
    contrast = (
        (planet_radius_m / star_radius_m) ** 2 * planet_signal / star_signal
    )
    if not math.isfinite(contrast) or contrast <= 0:
        raise ValueError(
            ("Contrast is outside the " "representable positive range.")
        )
    return float(contrast)


def contrast_from_models(
    star_model,
    planet_model,
    star_radius_m,
    planet_radius_m,
    response,
    base=".",
    detector="photon",
):
    """Load two explicitly selected models and return their band contrast.

    Parameters
    ----------
    star_model, planet_model : dict
        Spectrum mappings accepted by load_spectrum. File paths are
        resolved against base; this function performs no catalog
        lookup,
        nearest-neighbor selection, or interpolation.
    star_radius_m, planet_radius_m : float
        Positive emitting radii in metres.
    response : Response
        Validated filter/sky curves and integration interval.
    base : str or pathlib.Path, optional
        Parent directory for relative native model paths; default is
        ".".
    detector : {"photon", "energy"}, optional
        Detector weighting passed to contrast_from_spectra.

    Returns
    -------
    float
        Dimensionless companion/star thermal flux ratio.

    Raises
    ------
    FileNotFoundError
        If a selected model or wavelength file is not installed.
    ImportError
        If a selected native reader needs an unavailable dependency.
    ValueError
        If units, coordinates, coverage, or radii fail validation.

    Notes
    -----
    Call load_spectrum separately to inspect the returned metadata,
    including hashes and header-verified versus declared parameters.
    No network access or model downloads occur here.
    """
    # Native readers validate source conventions before any ratio
    # is formed.
    star = load_spectrum(star_model, base)
    planet = load_spectrum(planet_model, base)
    return contrast_from_spectra(
        star, planet, star_radius_m, planet_radius_m, response, detector
    )


def planet_contrast(
    star_wave,
    star_flux,
    planet_wave,
    planet_flux,
    *,
    star_radius_m,
    planet_radius_m,
    response,
    wavelength_unit,
    flux_unit,
    wavelength_frame="vacuum",
    flux_convention="surface_flux",
    detector="photon",
    planet_wavelength_unit=None,
    planet_flux_unit=None,
    planet_flux_convention=None,
):
    """Return a scalar contrast from already available surface-flux arrays.

    Parameters
    ----------
    star_wave, planet_wave : array-like
        Positive wavelengths in the explicitly declared input units.
    star_flux, planet_flux : array-like
        Nonnegative surface flux densities matching the wavelength
        arrays.
    star_radius_m, planet_radius_m : float
        Positive emitting radii in metres.
    response : Response
        Validated filter and sky transmission curves.
    wavelength_unit : str
        Stellar wavelength unit: um, nm, angstrom, m, or cm.
    flux_unit : str
        Stellar flux unit token accepted by to_surface_flambda.
    wavelength_frame : {"vacuum", "air"}, optional
        Common coordinate convention of both spectra. This function
        does
        not convert frames; it must agree with the response
        metadata.
    flux_convention : {"surface_flux", "eddington_flux"}, optional
        Stellar emergent surface flux or Eddington H_lambda
        convention.
    detector : {"photon", "energy"}, optional
        Detector weighting used by the integral.
    planet_wavelength_unit, planet_flux_unit : str or None, optional
        Planet unit overrides; None uses the corresponding stellar
        unit.
    planet_flux_convention : str or None, optional
        Planet convention override; None uses flux_convention.

    Returns
    -------
    float
        Dimensionless companion/star contrast for a downstream
        injector.

    Raises
    ------
    ValueError
        If inputs are invalid, units unsupported, signals
        nonpositive,
        wavelength frames inconsistent, or band coverage incomplete.

    Notes
    -----
    This preserves the original array interface. No PSF normalization,
    image injection, or observational recovery is performed here.
    """
    # Convert both bodies explicitly; the planet may use different
    # input units.
    star_wave_um, star_flambda = to_surface_flambda(
        star_wave, star_flux, wavelength_unit, flux_unit, flux_convention
    )
    planet_wave_um, planet_flambda = to_surface_flambda(
        planet_wave,
        planet_flux,
        planet_wavelength_unit or wavelength_unit,
        planet_flux_unit or flux_unit,
        planet_flux_convention or flux_convention,
    )
    # Both input arrays must actually use the declared common
    # wavelength frame.
    metadata = {"wavelength_frame": wavelength_frame}
    star = Spectrum(star_wave_um, star_flambda, dict(metadata))
    planet = Spectrum(planet_wave_um, planet_flambda, dict(metadata))
    return contrast_from_spectra(
        star, planet, star_radius_m, planet_radius_m, response, detector
    )


def _require_optional(package, imported_object):
    """Raise ImportError when a selected reader's dependency is missing."""
    if imported_object is None:
        raise ImportError(
            f"The selected native format requires {package}. "
            "Install it in the interpreter running this module."
        )


def _header_quantity(params, key, target_unit):
    """Convert one native quantity record to an explicit target unit.

    params is a header dictionary of {value, unit} records. Return float
    in target_unit; malformed records raise ValueError or unit errors.
    """
    # Require an explicit quantity record rather than inferring a unit.
    value = params[key]
    if not isinstance(value, dict) or not {"value", "unit"}.issubset(value):
        raise ValueError("Missing native quantity/unit for " + key)
    return float(
        (float(value["value"]) * u.Unit(value["unit"])).to_value(target_unit)
    )


def _band_integrand(spectrum, response, q, detector):
    """Evaluate transmitted surface F_lambda with explicit detector
    weighting.

    spectrum is a Spectrum, response a Response, and q wavelengths in
    um.
    detector is photon or energy (validated by integrate_band). Return
    photons m^-2 s^-1 um^-1 or W m^-2 um^-1; only the photon-energy
    factor
    converts q to metres. Both transmission components stay separate.
    """
    # Integrate surface F_lambda times two separately interpolated
    # transmissions.
    v = np.interp(
        q, spectrum.wavelength_um, spectrum.surface_flux
    ) * response.at(q)
    return v * q * 1e-06 / (H * C) if detector == "photon" else v


def _filename_parameters(name, fmt):
    """Decode native Sonora temperature, gravity, and composition labels.

    Parameters
    ----------
    name : str
        Published filename, not a renamed export.
    fmt : str
        Explicit Sonora Elf Owl or Flame Skimmer native format
        token.

    Returns
    -------
    dict
        Temperature in K, log10 gravity in cm/s^2, metallicity in
        dex,
        relative-to-solar C/O, and chemistry/mixing labels where
        present.

    Raises
    ------
    ValueError
        If the filename is inconsistent with the selected release.
    """
    if fmt == "sonora_flame_skimmer_v1_r3k":
        match = re.fullmatch(
            "spectra_(?:logkzz_"
            + _NUMBER
            + "_)?teff_"
            + _NUMBER
            + "_grav_"
            + _NUMBER
            + "_mh_"
            + _NUMBER
            + "_co_"
            + _NUMBER
            + "_R3000\\.nc",
            name,
        )
    else:
        match = re.fullmatch(
            "spectra_logzz_"
            + _NUMBER
            + "_teff_"
            + _NUMBER
            + "_grav_"
            + _NUMBER
            + "_mh_"
            + _NUMBER
            + "_co_"
            + _NUMBER
            + "\\.nc",
            name,
        )
    if not match:
        raise ValueError(
            (
                "Retain the published filename "
                "for the selected Sonora "
                "release and resolution."
            )
        )
    k, t, g, m, co = match.groups()
    t, g, m, co = map(float, (t, g, m, co))
    if t <= 0 or g <= 0 or co <= 0:
        raise ValueError("Invalid native atmosphere coordinates.")
    # Native gravity is m/s^2; multiply by 100 before taking
    # log10(g_cgs).
    params = {
        "teff_K": t,
        "logg_cgs": math.log10(g * 100),
        "metallicity_dex": m,
        "co_relative_solar": co,
        "chemistry": "equilibrium" if k is None else "disequilibrium",
    }
    if k is not None:
        params["logKzz_cgs"] = float(k)
    return params


def _load_sonora_netcdf(path, fmt):
    """Read native Sonora NetCDF surface emission and its provenance.

    path is a Path; fmt selects Elf Owl v2 or Flame Skimmer v1 R3000.
    Require matching 1D wavelength/flux variables and documented native
    units. Return wavelength, flux, unit tokens, and header/filename
    metadata. Missing dependencies raise ImportError; invalid native
    layouts raise ValueError. C/O is relative to solar, not absolute.
    """
    _require_optional("netCDF4", Dataset)
    _require_optional("astropy", u)
    # Cross-check filename coordinates against the native
    # parameter metadata.
    filename = _filename_parameters(path.name, fmt)
    flame = fmt == "sonora_flame_skimmer_v1_r3k"
    flux_name = "flux_emission" if flame else "flux"
    with Dataset(path, "r") as ds:
        if not {"wavelength", flux_name}.issubset(ds.variables):
            raise ValueError(
                ("Missing native wavelength or " "emission flux variable.")
            )
        wvar, fvar = (ds.variables["wavelength"], ds.variables[flux_name])
        if wvar.dimensions != fvar.dimensions or len(wvar.dimensions) != 1:
            raise ValueError(
                (
                    "Native spectrum must use "
                    "matching 1D wavelength "
                    "coordinates."
                )
            )
        wu = getattr(wvar, "units", "")
        fu = getattr(fvar, "units", "")
        try:
            wfactor = u.Unit(wu).to(u.um)
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=".*multiple slashes.*",
                    category=u.UnitsWarning,
                )
                ffactor = u.Unit(fu).to(u.erg / u.s / u.cm**2 / u.cm)
        except (ValueError, TypeError) as exc:
            raise ValueError(
                (
                    "NetCDF wavelength and flux "
                    "units are absent or "
                    "incompatible with this "
                    "release."
                )
            ) from exc
        if not math.isclose(wfactor, 1, rel_tol=1e-12) or not math.isclose(
            ffactor, 1, rel_tol=1e-12
        ):
            raise ValueError(
                (
                    "Native NetCDF adapter "
                    "requires microns and "
                    "erg/s/cm²/cm. Use an explicit "
                    "ASCII adapter for converted "
                    "units."
                )
            )
        w, f = (wvar[:], fvar[:])
        if np.ma.is_masked(w) or np.ma.is_masked(f):
            raise ValueError(
                (
                    "Masked native spectral "
                    "samples cannot be silently "
                    "discarded."
                )
            )
        raw = getattr(ds, "planet_params", None)
        if not isinstance(raw, str):
            raise ValueError("Native planet_params JSON is required.")
        params = json.loads(raw)
        if not isinstance(params, dict):
            raise ValueError("planet_params must be a JSON object.")
        headers = {
            "wavelength_units": wu,
            "flux_units": fu,
            "planet_params": params,
        }
        if flame:
            resolution = getattr(ds, "rebinned_resolution", "")
            if resolution != "R=3000":
                raise ValueError(
                    (
                        "Select the R=3000 quickstart "
                        "release; do not mix native "
                        "resolutions in one family."
                    )
                )
            headers["rebinned_resolution"] = resolution
        # Retain the release-specific meaning of each native
        # header quantity.
        measured = {}
        if flame:
            measured.update(
                teff_K=_header_quantity(params, "effective_temp", u.K),
                logg_cgs=math.log10(
                    _header_quantity(params, "gravity", u.cm / u.s**2)
                ),
                metallicity_dex=float(params["logmh"]),
            )
            co = params.get("cto", {})
            if (
                not isinstance(co, dict)
                or "relative" not in str(co.get("unit", "")).lower()
                or "0.458" not in str(co["unit"])
            ):
                raise ValueError(
                    (
                        "Flame Skimmer quickstart C/O "
                        "must be explicitly relative "
                        "to solar (=0.458)."
                    )
                )
            measured["co_relative_solar"] = float(co["value"])
        else:
            measured.update(
                teff_K=_header_quantity(params, "teff", u.K),
                logg_cgs=math.log10(
                    _header_quantity(params, "logg", u.cm / u.s**2)
                ),
                metallicity_dex=float(params["mh"]),
                co_relative_solar=float(params["cto"]),
                logKzz_cgs=float(params["logkzz"]),
            )
        if "logKzz_cgs" in filename and flame:
            for key in ("logkzz", "logKzz"):
                if key in params and isinstance(params[key], (float, int)):
                    measured["logKzz_cgs"] = float(params[key])
    # Preserve C/O normalization and source identity for
    # downstream inspection.
    extra = {
        ("parameters_verified_from_filen" "ame"): filename,
        ("parameters_verified_from_heade" "r"): measured,
        "native_headers": headers,
        ("native_units_verified_from_hea" "der"): True,
        "co_convention": (
            "C/O divided by solar C/O "
            "(0.458); not the absolute "
            "atomic ratio."
        ),
        "source_release": (
            ("Flame Skimmer v1 R3000 " "quickstart")
            if flame
            else (
                "Elf Owl v2 (user-selected "
                "release; the file itself does "
                "not independently establish "
                "version)"
            )
        ),
        "source_record_url": (
            ("https://zenodo.org/records/200" "30439")
            if flame
            else ("https://zenodo.org/records/151" "50874")
        ),
        "native_resolving_power": 3000 if flame else None,
    }
    return (
        np.asarray(w, float),
        np.asarray(f, float),
        "um",
        "erg_cm2_s_cm",
        extra,
    )


def _read_two_column_text(path):
    """Read two native numeric columns below their original text header.

    path is a Path, optionally gzip compressed. Return a 2D numeric
    array
    and preceding header lines. Fortran D exponents are accepted. Raise
    ValueError if no valid two-column start exists in the first 40
    lines.
    """
    # Gzip and Fortran exponent handling must not alter numerical
    # values.
    text = (
        gzip.open(path, "rt").read()
        if path.suffix == ".gz"
        else path.read_text()
    )
    text = text.replace("D+", "E+").replace("D-", "E-")
    lines = text.splitlines()
    for i, line in enumerate(lines[:40]):
        tokens = line.split()
        if len(tokens) == 2:
            try:
                a, b = map(float, tokens)
            except ValueError:
                continue
            if a > 0 and b >= 0:
                return (
                    np.loadtxt(io.StringIO(text), skiprows=i, ndmin=2),
                    lines[:i],
                )
    raise ValueError("No two-column spectrum found below the native header.")


def _load_native_format(path, fmt, spec):
    """Dispatch additional native text, FITS, VOTable, or NetCDF formats.

    path is a Path; fmt and spec explicitly select a release/convention.
    Return wavelength, flux, wavelength-unit token, flux-unit token, and
    provenance. Do not convert a different release by changing its
    label.
    Invalid units/layouts raise ValueError; unavailable readers raise
    ImportError through load_spectrum's dependency checks.
    """
    filename = {}
    header_parameters = {}
    extra = {}
    if fmt in ("sonora_elf_owl_v2", "sonora_flame_skimmer_v1_r3k"):
        return _load_sonora_netcdf(path, fmt)
    if fmt in ("sonora_cholla", "sonora_diamondback", "atmo2020_ascii"):
        data, header = _read_two_column_text(path)
        wave, flux = (data[:, 0], data[:, 1])
        # Inspect the native header before assigning a
        # release-specific flux unit.
        compact = "".join(" ".join(header).lower().split()).replace("^", "")
        if fmt in ("sonora_cholla", "sonora_diamondback"):
            unit_text = re.sub(
                "\\s*/\\s*", "/", " ".join(header).lower().replace("^", "")
            )
            flux_units = re.findall(
                ("(?<![a-z0-9])w/(?:m2|m²)/([a-z" "0-9µμ/]+)"), unit_text
            )
            if (
                not ("micron" in compact or "wavelength(um)" in compact)
                or not flux_units
                or any((unit != "m" for unit in flux_units))
            ):
                raise ValueError(
                    (
                        "Expected native micron and "
                        "W/m2/m emission header. "
                        "Supply an explicit ASCII "
                        "adapter for a different "
                        "release."
                    )
                )
            # Return explicit unit tokens to the shared converter below.
            wave_unit, flux_unit = ("um", "W_m2_m")
            name = (
                path.name.removesuffix(".gz")
                if hasattr(str, "removesuffix")
                else path.name[:-3] if path.name.endswith(".gz") else path.name
            )
            if fmt == "sonora_cholla":
                match = re.fullmatch(
                    ("(\\d+)K_(\\d+(?:\\.\\d+)?)g_logkzz" "(\\d+)\\.spec"),
                    name,
                )
                if not match:
                    raise ValueError(
                        (
                            "Unrecognized Cholla filename; "
                            "retain its published model "
                            "filename."
                        )
                    )
                t, g, k = map(float, match.groups())
                filename = {
                    "teff_K": t,
                    "logg_cgs": math.log10(g * 100),
                    "logKzz_cgs": k,
                    "metallicity_dex": 0.0,
                }
            else:
                match = re.fullmatch(
                    (
                        "t(\\d+)g(\\d+(?:\\.\\d+)?)(f\\d+(?:"
                        "\\.\\d+)?|nc)_m([+-]?\\d+(?:\\.\\d+"
                        ")?)_co(\\d+(?:\\.\\d+)?)(?:\\.[^.]"
                        "+)?"
                    ),
                    name,
                )
                if not match:
                    raise ValueError(
                        (
                            "Unrecognized Diamondback "
                            "filename; retain the "
                            "published t/g/fsed/metal/C-O "
                            "filename."
                        )
                    )
                t, g, cloud, m, co = match.groups()
                filename = {
                    "teff_K": float(t),
                    "logg_cgs": math.log10(float(g) * 100),
                    "fsed": (
                        "cloud_free" if cloud == "nc" else float(cloud[1:])
                    ),
                    "metallicity_dex": float(m),
                    "co_relative_solar": float(co),
                }
                extra.update(
                    co_convention=(
                        "C/O relative to solar (solar " "C/O=0.458)."
                    ),
                    source_release="Sonora Diamondback",
                    source_record_url=(
                        "https://zenodo.org/records/127" "35103"
                    ),
                )
        else:
            wave_unit, flux_unit = ("um", "W_m2_um")
            if header and (
                not (
                    ("micron" in compact or "wavelength[um]" in compact)
                    and ("w/m2/micron" in compact or "w/m2/um" in compact)
                )
            ):
                raise ValueError(
                    (
                        "Original ATMO emission header "
                        "must specify micron and "
                        "W/m2/micron. Use explicit "
                        "ASCII units for a converted "
                        "release."
                    )
                )
            match = re.fullmatch(
                (
                    "spec_t(\\d+)_lg([\\d.]+)_(ceq|ne"
                    "q_weak|neq_strong)\\.txt"
                ),
                path.name,
                re.IGNORECASE,
            )
            if match:
                filename = {
                    "teff_K": float(match[1]),
                    "logg_cgs": float(match[2]),
                    "chemistry": match[3].lower(),
                }
            extra["native_units_verified_from_header"] = bool(header)
            extra["unit_authority"] = (
                "Original ATMO 2020 emission "
                "layout documented by STScI; "
                "user must select this "
                "original release explicitly."
            )
        extra["original_header"] = header
    # STScI FITS units are read from the table, including the FLAM
    # convention.
    elif fmt == "atmo2020_stsci_fits":
        with fits.open(path, memmap=False) as hdus:
            tables = [
                h
                for h in hdus
                if isinstance(h, (fits.BinTableHDU, fits.TableHDU))
                and h.data is not None
            ]
            if len(tables) != 1:
                raise ValueError("Expected one ATMO STScI spectral table.")
            table = tables[0]
            names = {n.lower(): n for n in table.columns.names}
            if not {"wavelength", "flux"}.issubset(names):
                raise ValueError("ATMO STScI table needs WAVELENGTH and FLUX.")
            wname, fname = (names["wavelength"], names["flux"])
            wu = table.columns[wname].unit
            fu = table.columns[fname].unit
            if not wu or not fu:
                raise ValueError("FITS spectral column units are required.")
            wunit = u.AA if wu.strip().upper() == "ANGSTROMS" else u.Unit(wu)
            funit = (
                u.erg / u.s / u.cm**2 / u.AA
                if fu.strip().upper() == "FLAM"
                else u.Unit(fu)
            )
            wave = (np.asarray(table.data[wname], float) * wunit).to_value(
                u.um
            )
            flux = (np.asarray(table.data[fname], float) * funit).to_value(
                u.W / u.m**2 / u.um
            )
            extra["fits_column_units"] = {"wavelength": wu, "flux": fu}
            for key, candidates in {
                "teff_K": ("TEFF", "PHXTEFF"),
                "logg_cgs": ("LOGG", "PHXLOGG"),
            }.items():
                for h in hdus:
                    for candidate in candidates:
                        if candidate in h.header:
                            header_parameters[key] = float(h.header[candidate])
        wave_unit, flux_unit = ("um", "W_m2_um")
        match = re.fullmatch(
            ("spec_t(\\d+)_lg([\\d.]+)_(ceq|ne" "q_weak|neq_strong)\\.fits"),
            path.name,
        )
        if match:
            filename = {
                "teff_K": float(match[1]),
                "logg_cgs": float(match[2]),
                "metallicity_dex": 0.0,
                "chemistry": match[3],
            }
    elif fmt == "bt_settl_svo":
        if spec.get("wavelength_frame") != "air":
            raise ValueError(
                (
                    "SVO BT-Settl CIFIST uses air "
                    "wavelengths. Choose air "
                    "explicitly; this adapter does "
                    "not convert frames."
                )
            )
        table = Table.read(path, format="votable")
        names = {n.lower(): n for n in table.colnames}
        if not {"wavelength", "flux"}.issubset(names):
            raise ValueError(
                (
                    "Expected SVO emission VOTable "
                    "Wavelength and Flux columns."
                )
            )
        wcol, fcol = (table[names["wavelength"]], table[names["flux"]])
        if wcol.unit is None or fcol.unit is None:
            raise ValueError(
                ("SVO wavelength and " "surface-flux units are " "required.")
            )
        wave = u.Quantity(wcol).to_value(u.um)
        flux = u.Quantity(fcol).to_value(u.W / u.m**2 / u.um)
        wave_unit, flux_unit = ("um", "W_m2_um")
        extra["unit_authority"] = (
            "SVO VOTable column units; "
            "BT-Settl CIFIST air "
            "convention. Atmospheric "
            "parameters are user "
            "declarations unless "
            "separately verified."
        )
    else:
        raise ValueError("Unsupported additional adapter.")
    extra.update(
        parameters_verified_from_filename=filename,
        parameters_verified_from_header=header_parameters,
    )
    return (wave, flux, wave_unit, flux_unit, extra)


def _bobcat_coordinates(filename, reported):
    """Verify Bobcat filename adjacency labels against reported
    coordinates.

    filename is the published basename; reported contains header values.
    Return connectivity metadata, or None for an unrecognized filename.
    Permit only the original 1 K and 1 m/s^2 label offsets. Gravity and
    temperature remain at reported physical values, without snapping.
    Raise ValueError for conflicting coordinates/composition.
    """
    match = re.fullmatch(
        (
            "sp_t(\\d+)g(\\d+(?:\\.\\d+)?)nc_m("
            "[+-]?\\d+(?:\\.\\d+)?)(?:_co(\\d+("
            "?:\\.\\d+)?))?(?:\\.gz)?"
        ),
        filename,
    )
    if not match:
        return None
    nominal_temperature, filename_gravity, metallicity = map(
        float, match.groups()[:3]
    )
    relative_co = float(match[4]) if match[4] is not None else 1.0
    actual_temperature = float(reported["teff_K"])
    # Compare in m/s^2 before deciding whether label rounding is
    # acceptable.
    actual_gravity = 10.0 ** float(reported["logg_cgs"]) / 100.0
    if (
        abs(actual_temperature - nominal_temperature)
        > MAX_TEMPERATURE_OFFSET_K + 1e-09
    ):
        raise ValueError(
            "Bobcat "
            + filename
            + (
                ": filename temperature %g K, "
                "reported header %g K; allowed "
                "identity offset is at most %g "
                "K."
            )
            % (
                nominal_temperature,
                actual_temperature,
                MAX_TEMPERATURE_OFFSET_K,
            )
        )
    if (
        abs(actual_gravity - filename_gravity)
        > MAX_GRAVITY_LABEL_OFFSET_M_S2 + 1e-09
    ):
        raise ValueError(
            "Bobcat "
            + filename
            + (
                ": filename gravity %g m/s2, "
                "reported header %g m/s2; the "
                "label/header difference "
                "exceeds %g m/s2."
            )
            % (filename_gravity, actual_gravity, MAX_GRAVITY_LABEL_OFFSET_M_S2)
        )
    for key, value in (
        ("metallicity_dex", metallicity),
        ("co_relative_solar", relative_co),
    ):
        if not math.isclose(
            float(reported[key]), value, rel_tol=1e-08, abs_tol=1e-08
        ):
            raise ValueError(
                "Bobcat "
                + filename
                + ": filename "
                + key
                + ("=%g disagrees with reported " "header %g.")
                % (value, reported[key])
            )
    return {
        "nominal_teff_K": nominal_temperature,
        "filename_gravity_m_s2": filename_gravity,
        "reported_teff_K": actual_temperature,
        "reported_gravity_m_s2": actual_gravity,
        ("temperature_header_minus_filen" "ame_K"): actual_temperature
        - nominal_temperature,
        "meaning": (
            "Filename temperature is a "
            "connectivity index only. "
            "Interpolation and displayed "
            "physical parameters retain "
            "reported header temperature "
            "and gravity."
        ),
    }


def _read_ck04(path, spec):
    """Read one exact native STScI CK04 FITS gravity column.

    path is a Path; spec may select gravity_column (default g45).
    Return Angstrom wavelengths, FLAM emergent surface flux, and native
    parameter/unit provenance. Require compatible filename and
    TEFF/LOG_Z
    headers. Reject missing or all-zero gravity columns with ValueError.
    """
    match = re.fullmatch("ck([mp])(\\d{2})_(\\d+)\\.fits", path.name)
    if not match:
        raise ValueError(
            ("Keep the published CK04 " "cksmh_temperature.fits " "filename.")
        )
    sign, metal, temperature = match.groups()
    named = {
        "teff_K": float(temperature),
        "metallicity_dex": (1 if sign == "p" else -1) * float(metal) / 10,
    }
    # CK04 stores exact gravity choices as separate spectral columns.
    column = spec.get("gravity_column", "g45")
    if not re.fullmatch("g\\d{2}", column):
        raise ValueError(
            ("Choose an exact CK04 gravity " "column, for example g45.")
        )
    with fits.open(path, memmap=False) as hdus:
        h = hdus[0].header
        if "TEFF" not in h or "LOG_Z" not in h:
            raise ValueError(
                ("CK04 native TEFF and LOG_Z " "headers are required.")
            )
        # The chosen gravity column defines log g;
        # temperature/metallicity use headers.
        measured = {
            "teff_K": float(h["TEFF"]),
            "metallicity_dex": float(h["LOG_Z"]),
            "logg_cgs": float(column[1:]) / 10,
        }
        for key, value in named.items():
            if not math.isclose(
                value, measured[key], rel_tol=1e-10, abs_tol=1e-10
            ):
                raise ValueError(
                    ("CK04 filename/header " "parameter mismatch for ") + key
                )
        tables = [x for x in hdus if isinstance(x, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected one CK04 native binary table.")
        table = tables[0]
        names = {n.lower(): n for n in table.columns.names}
        if "wavelength" not in names or column not in names:
            raise ValueError(
                (
                    "Missing wavelength or "
                    "requested native gravity "
                    "column."
                )
            )
        wname, fname = (names["wavelength"], names[column])
        wu, fu = (table.columns[wname].unit, table.columns[fname].unit)
        if not wu or not fu:
            raise ValueError("CK04 wavelength and flux units are required.")
        wunit = u.AA if wu.upper() in ("ANGSTROM", "ANGSTROMS") else u.Unit(wu)
        funit = (
            u.erg / u.s / u.cm**2 / u.AA
            if fu.upper() == "FLAM"
            else u.Unit(fu)
        )
        if not math.isclose(wunit.to(u.AA), 1.0, rel_tol=1e-12):
            raise ValueError("Native CK04 wavelengths must be Angstroms.")
        if not math.isclose(
            funit.to(u.erg / u.s / u.cm**2 / u.AA), 1.0, rel_tol=1e-12
        ):
            raise ValueError("Native CK04 flux must be FLAM surface flux.")
        # These native arrays remain Angstrom/FLAM until the
        # common conversion step.
        wave = np.asarray(table.data[wname], float)
        flux = np.asarray(table.data[fname], float)
        if not np.any(flux > 0):
            raise ValueError(
                ("This CK04 gravity column is " "unavailable (all zero).")
            )
    return (
        wave,
        flux,
        {
            ("parameters_verified_from_filen" "ame"): named,
            ("parameters_verified_from_heade" "r"): measured,
            ("native_units_verified_from_hea" "der"): True,
            "selected_gravity_column": column,
            "source_release": "STScI CK04 ATLAS9",
            "fits_column_units": {"wavelength": wu, "flux": fu},
            "surface_flux_authority": (
                "STScI native HISTORY and CK04 "
                "atlas documentation; scale "
                "observed flux by (R/D)^2."
            ),
            "wavelength_frame_authority": (
                "Vacuum convention adopted for "
                "this integration; the FITS "
                "table does not independently "
                "specify air/vacuum."
            ),
            "source_record_url": (
                "https://www.stsci.edu/hst/inst"
                "rumentation/reference-data-for"
                "-calibration-and-tools/astrono"
                "mical-catalogs/castelli-and-ku"
                "rucz-atlas"
            ),
        },
    )


def _file_record(path):
    """Return absolute path, byte size, and SHA-256 for a native file.

    path must be a Path to an installed file. File-system failures
    propagate. Hashes record source identity, not physical accuracy.
    """
    # Hash in chunks so provenance does not require a second full-
    # memory copy.
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(part)
    return {
        "path": str(path.resolve()),
        "sha256": digest.hexdigest(),
        "bytes": path.stat().st_size,
    }


def _check_parameters(declared, measured):
    """Reject conflicts between overlapping declared and measured keys.

    Dictionaries may contain numbers or categorical strings. Compare
    numbers at the original relative/absolute tolerances and strings
    exactly. Unmeasured declarations do not become independently
    verified.
    """
    for key in declared.keys() & measured.keys():
        equal = (
            math.isclose(
                float(declared[key]),
                float(measured[key]),
                rel_tol=1e-06,
                abs_tol=0.002,
            )
            if isinstance(measured[key], (int, float))
            else declared[key] == measured[key]
        )
        if not equal:
            raise ValueError(
                ("Declared %s=%s disagrees with " "model header value %s.")
                % (key, declared[key], measured[key])
            )


def _sha256_file(path):
    """Return the streamed SHA-256 hex digest of an installed file.

    path is a str or Path; missing/unreadable files raise OSError.
    """
    out = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            out.update(chunk)
    return out.hexdigest()


def _validate_transmission(wave, transmission, label):
    """Validate a dimensionless transmission component on a wavelength
    grid.

    wave is in um; transmission must be finite in [0,1], with nonzero
    support. Return float arrays. label identifies the component in a
    ValueError for invalid shapes, values, or wavelength ordering.
    """
    x, t = (np.asarray(wave, float), np.asarray(transmission, float))
    if x.ndim != 1 or t.shape != x.shape or len(x) < 3:
        raise ValueError(
            label
            + (": expected matching 1D arrays " "with at least three samples.")
        )
    if (
        not np.isfinite(x).all()
        or not np.isfinite(t).all()
        or np.any(x <= 0)
        or np.any(np.diff(x) <= 0)
    ):
        raise ValueError(
            label
            + (": wavelengths must be finite, " "positive and increasing.")
        )
    if np.any(t < 0) or np.any(t > 1) or (not np.any(t > 0)):
        raise ValueError(
            label
            + (": transmission must be in " "[0,1] with positive support.")
        )
    return (x, t)


def _integration_nodes(spectrum, response):
    """Return the union of spectrum/filter/sky knots and band edges in um.

    spectrum and response must share a wavelength frame and cover the
    entire declared interval. Raise ValueError otherwise. Preserving all
    knots prevents a coarse plotting grid from changing the band
    integral.
    """
    lo, hi = response.limits_um
    x = spectrum.wavelength_um
    if lo < x[0] or hi > x[-1]:
        raise ValueError(
            (
                "Spectrum must cover %.6g–%.6g "
                "um; this file covers "
                "%.6g–%.6g um. No "
                "extrapolation."
            )
            % (lo, hi, x[0], x[-1])
        )
    if spectrum.metadata.get(
        "wavelength_frame", "vacuum"
    ) != response.metadata.get("wavelength_frame", "vacuum"):
        raise ValueError(
            ("Spectrum and response " "wavelength conventions differ.")
        )
    # Merge all native knots instead of resampling onto an
    # arbitrary coarse grid.
    arrays = (x, response.filter_wave_um, response.sky_wave_um)
    return np.unique(
        np.concatenate(([lo, hi], *[a[(a > lo) & (a < hi)] for a in arrays]))
    )
