# Atmosphere contrast: API and native-model guide

Run commands in this guide from the repository root. This guide preserves
the original review module's interface and scientific validation scope.
For the repository overview and additional file examples, see
[README.md](../README.md).

This extraction supplies the **thermal companion/star brightness ratio**
needed by a downstream injector. It follows the request to separate the
brightness calculation from Atmosphere Lab and make that code readable.

The source is the supplied `0.8.5-single-file.1` standalone export. This is
an independent review module, version `0.1.0`, not an update to the installed
app and not a claim that the app is version 0.1.0. It retains the selected
export's native unit/header rules and numerical band integration.

## Files

- `atmosphere_contrast.py`: the importable calculation and native readers.
- `example_usage.py`: self-contained analytic example, no models required.
- `test_contrast.py`: manufactured-data tests, including native layouts.
- `VALIDATION.json`: recorded preparation checks and their limited scope.
- `requirements.txt`: NumPy for the array interface and ASCII/Bobcat readers.
- `requirements_native.txt`: add Astropy and netCDF4 for all native readers.
- `requirements_dev.txt`: add the formatting/style tools used in review.

There is no Streamlit dependency. All imports are at module scope. Classes
and functions are at the top level, except ordinary class methods. There are
no builder namespaces or nested functions. Each class, function, and method
has a docstring. Block comments explain important numerical/physical steps.

The extraction intentionally excludes plotting, bundled old results, catalog
selection, interpolation, heating, PSFs, sensitivity import, and image
injection/recovery. It predicts brightness; it does not produce a calibrated
instrument sensitivity or a detection probability.

## Install and run

Use Python 3.9 or newer in a clean environment:

```sh
python -m pip install -r requirements.txt
python example_usage.py
python -m unittest -v
```

The example prints `0.05`, the analytically expected ratio for deliberately
constant spectra. It is not a calculated atmosphere-grid prediction.

For FITS, VOTable, and NetCDF spectra, install the optional dependencies:

```sh
python -m pip install -r requirements_native.txt
python -m unittest -v
```

The native fixtures cover the readers whose optional dependencies are
installed. The preparation run installed both optional packages and exercised
all 11 format tokens listed below. Syntax was checked against Python 3.9;
the preparation interpreter was Python 3.12. No Python 3.9 runtime execution
or real native-grid recomputation is claimed by these preparation tests.

## The scalar interface for an injector

```python
from atmosphere_contrast import planet_contrast

contrast = planet_contrast(
    star_wave,
    star_surface_flux,
    planet_wave,
    planet_surface_flux,
    star_radius_m=star_radius_m,
    planet_radius_m=planet_radius_m,
    response=response,
    wavelength_unit="um",
    flux_unit="W_m2_um",
    wavelength_frame="vacuum",
)
```

The result is one positive, dimensionless **planet/star flux ratio**, not a
percentage or magnitude. Supply it to the downstream injector according to
that injector's existing brightness/PSF normalization convention. This module
does not call pyNOMIC or modify its source.

The planet may use different units through `planet_wavelength_unit`,
`planet_flux_unit`, and `planet_flux_convention`. Both spectra must use the
declared common air/vacuum wavelength frame; frame conversion is not automatic.

## Actual atmosphere files and the supplied Nprime response

Model files and transmission curves remain external. This module does not
download them, and the ZIP does not contain model grids or response tables.
The paths below are placeholders to replace with existing local files:

```python
from atmosphere_contrast import contrast_from_models, load_nprime

response = load_nprime("/path/to/pinned/filterdata")

star_model = {
    "file": "/path/to/ckp00_5000.fits",
    "format": "ck04_stsci_fits",
    "gravity_column": "g45",
    "wavelength_frame": "vacuum",
    "parameters": {
        "teff_K": 5000.0,
        "logg_cgs": 4.5,
        "metallicity_dex": 0.0,
    },
}
planet_model = {
    "file": "/path/to/sp_t1200g100nc_m0.0.gz",
    "format": "sonora_bobcat",
    "wavelength_frame": "vacuum",
}

contrast = contrast_from_models(
    star_model,
    planet_model,
    star_radius_m=6.957e8,
    planet_radius_m=7.1492e7,
    response=response,
)
print(contrast)
```

For Bobcat, inspect the native header before adding optional `parameters`:
the reported physical temperature/gravity can differ slightly from the
filename labels. The loader retains the original offset rules and checks
overlapping declarations against available header measurements. It does not
silently choose a nearest atmosphere or snap a physical coordinate.

`load_nprime` requires the reviewed response bytes at pyNOMIC commit
`e6c63c740a6fe1ca74e57228fafa3d7da8b6e67a`:

- `Nprime_filter_curve.txt`
- `MtGrahamTransmittance.npz`

The file hashes are checked. Tiny reviewed negative filter artifacts are
floored to zero and recorded. The full default filter interval is 9-14 um.
This is a relative filter shape and reference atmospheric transmission;
observation-specific PWV/airmass, telescope optics, QE, and sky emission are
not inferred. If the upstream curves change, review the new data rather than
bypassing their hashes.

Source: <https://github.com/vivekv7323/pyNOMIC/tree/main/filterdata>.

For a different response, construct a `Response` from explicit filter and
sky arrays instead. Wavelengths are in um, transmissions are fractions in
[0, 1], and both grids must span `limits_um`:

```python
from atmosphere_contrast import Response

response = Response(
    filter_wave_um=filter_wave_um,
    filter_transmission=filter_transmission,
    sky_wave_um=sky_wave_um,
    sky_transmission=sky_transmission,
    limits_um=(lower_um, upper_um),
    metadata={"wavelength_frame": "vacuum"},
)
```

Use `load_spectrum(mapping, base)` separately to inspect original header
records, declared versus verified coordinates, and source hashes. Pass the
resulting `Spectrum` objects to `contrast_from_spectra(...)` if already loaded.

## Supported native mappings

| `format` | Input convention / optional dependencies |
| --- | --- |
| `ascii` | Explicit units, convention, and optional table columns |
| `sonora_bobcat` | Native um and cgs F_nu; text/gzip |
| `ck04_stsci_fits` | Native Angstrom/FLAM and exact gravity column; Astropy |
| `phoenix_hires` | Flux FITS plus separate wavelength FITS; Astropy |
| `sonora_cholla` | Native text um and W/m2/m |
| `sonora_diamondback` | Native text um and W/m2/m |
| `atmo2020_ascii` | Original text um and W/m2/um |
| `atmo2020_stsci_fits` | STScI spectral-table units; Astropy |
| `bt_settl_svo` | SVO VOTable, explicitly air wavelengths; Astropy |
| `sonora_elf_owl_v2` | Native NetCDF layout; Astropy + netCDF4 |
| `sonora_flame_skimmer_v1_r3k` | Native R=3000 NetCDF; Astropy + netCDF4 |

These are specific release layouts, not a promise that every similarly named
download is compatible. Source unit checks and adopted wavelength-frame
limitations from the supplied export are retained. In particular, CK04 and
the pinned response adopt vacuum coordinates without independently encoding
that convention in the original tables. The newer app's BT-Settl STScI FITS
adapter is not part of this older export extraction; do not relabel that file
as `bt_settl_svo`.

For an explicit converted ASCII file, declare the convention:

```python
mapping = {
    "file": "converted_surface_spectrum.txt",
    "format": "ascii",
    "wavelength_unit": "um",
    "flux_unit": "W_m2_um",
    "flux_convention": "surface_flux",
    "wavelength_frame": "vacuum",
    "skiprows": 0,
}
```

## Physics and numerical method

Internally, every spectrum is surface F_lambda in W m^-2 um^-1. For a photon
response the surface band signal is

```text
S = integral F_lambda(lambda) * filter(lambda) * sky(lambda)
             * (lambda_metres / (h*c)) d(lambda_um)
contrast = (R_planet / R_star)^2 * S_planet / S_star
```

For energy weighting, omit lambda/(h*c). Emergent surface flux needs no
additional pi. Explicit Eddington H_lambda is multiplied by 4*pi during
conversion. Observed flux and specific intensity are not surface-flux input
conventions. A common distance and shared multiplicative response scale
cancel in the ratio.

The integrator merges all spectrum/filter/sky knots and the band edges.
Each factor is interpolated linearly on its own source grid. Their product
is cubic per interval; photon weighting makes it quartic. Three-point Gauss
and independent Boole quadrature integrate this tabulated representation to
floating-point accuracy. That is a numerical property, not a certificate of
physical accuracy across atmosphere grids.

No extrapolation or partial-band truncation is implicit. Radii are imposed
independently of atmospheric gravity, mass, and age. Reflected light,
irradiation/thermal lag, evolution, and detection probability are outside
this extraction's scope.

## Formatting and review

Use Python 3.12 for the pinned development tools. Python 3.9 remains the
minimum target for running the calculation module and examples.

```sh
python -m pip install -r requirements_dev.txt
python -m black --check --line-length 79 atmosphere_contrast.py test_contrast.py example_usage.py
python -m pycodestyle --max-line-length=79 atmosphere_contrast.py test_contrast.py example_usage.py
```

Conventions follow [PEP 8](https://peps.python.org/pep-0008/) and
[PEP 257](https://peps.python.org/pep-0257/). The extraction keeps all imports
at module scope, including guarded optional imports, and provides NumPy-style
parameter/return documentation for its public API.

No GitHub publication or pyNOMIC merge has been performed. Before merging,
use the destination repository's packaging, attribution, license, and test
conventions. Existing model/transmission data remain subject to their own
source documentation and attribution requirements.
