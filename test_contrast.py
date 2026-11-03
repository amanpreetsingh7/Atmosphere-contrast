"""Check analytic brightness identities and manufactured native-file readers.

Run ``python -m unittest -v`` from this directory. Native-reader tests use
temporary synthetic files, never downloaded atmosphere grids. Optional
native layouts are included when their dependencies are installed.
"""

import gzip
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import atmosphere_contrast as ac


def fixture_models(folder):
    """Write synthetic spectra in the supported native release layouts.

    Return explicit spectrum mappings. File headers/units test adapters;
    these data are manufactured constants, not scientific model spectra.
    """
    folder = Path(folder)
    wave = np.array([8.0, 9.0, 10.5, 12.0, 14.0, 15.0])
    fixtures = {}
    ascii_path = folder / "spectrum.txt"
    np.savetxt(ascii_path, np.column_stack([wave, np.ones(len(wave))]))
    fixtures["ascii"] = {
        "file": ascii_path.name,
        "format": "ascii",
        "wavelength_unit": "um",
        "flux_unit": "W_m2_um",
        "flux_convention": "surface_flux",
    }
    bobcat = folder / "sp_t1200g100nc_m0.0.gz"
    with gzip.open(bobcat, "wt") as stream:
        stream.write("1200 100 .28 0 1e5 0 1 1\n")
        stream.write("microns Flux (erg/cm^2/s/Hz)\n")
        np.savetxt(stream, np.column_stack([wave, np.ones(len(wave)) * 1e-12]))
    fixtures["sonora_bobcat"] = {
        "file": bobcat.name,
        "format": "sonora_bobcat",
    }
    for fmt, name, flux, header in (
        (
            "sonora_cholla",
            "1200K_100g_logkzz4.spec",
            1e6,
            "microns Flux (W/m2/m)\n",
        ),
        (
            "sonora_diamondback",
            "t1200g100f2_m0.0_co1.0.spec",
            1e6,
            "microns Flux (W/m2/m)\n",
        ),
        (
            "atmo2020_ascii",
            "spec_t1200_lg4.0_ceq.txt",
            1.0,
            "microns Flux (W/m2/micron)\n",
        ),
    ):
        path = folder / name
        with path.open("w") as stream:
            stream.write(header)
            np.savetxt(stream, np.column_stack([wave, wave * 0 + flux]))
        fixtures[fmt] = {"file": name, "format": fmt}

    # FITS/VOTable fixtures preserve native unit labels and gravity columns.
    if ac.fits is not None:
        fits = ac.fits
        columns = [
            fits.Column(
                name="WAVELENGTH",
                unit="Angstrom",
                format="D",
                array=wave * 1e4,
            ),
            fits.Column(
                name="g45",
                unit="FLAM",
                format="D",
                array=np.ones(len(wave)) * 0.1,
            ),
        ]
        primary = fits.PrimaryHDU()
        primary.header["TEFF"] = 5000.0
        primary.header["LOG_Z"] = 0.0
        path = folder / "ckp00_5000.fits"
        fits.HDUList(
            [primary, fits.BinTableHDU.from_columns(columns)]
        ).writeto(path)
        fixtures["ck04_stsci_fits"] = {
            "file": path.name,
            "format": "ck04_stsci_fits",
            "gravity_column": "g45",
        }
        path = folder / "phoenix.fits"
        primary = fits.PrimaryHDU(np.ones(len(wave)) * 1e7)
        primary.header["BUNIT"] = "erg/s/cm2/cm"
        primary.header["PHXTEFF"] = 5000.0
        primary.header["PHXLOGG"] = 4.5
        primary.header["PHXM_H"] = 0.0
        primary.writeto(path)
        wavelength_path = folder / "wavelength.fits"
        primary = fits.PrimaryHDU(wave * 1e4)
        primary.header["UNIT"] = "Angstrom"
        primary.writeto(wavelength_path)
        fixtures["phoenix_hires"] = {
            "file": path.name,
            "format": "phoenix_hires",
            "wavelength_file": wavelength_path.name,
        }
        path = folder / "spec_t1200_lg4.0_ceq.fits"
        table = fits.BinTableHDU.from_columns(
            [
                columns[0],
                fits.Column(
                    name="FLUX",
                    unit="FLAM",
                    format="D",
                    array=np.ones(len(wave)) * 0.1,
                ),
            ]
        )
        fits.HDUList([fits.PrimaryHDU(), table]).writeto(path)
        fixtures["atmo2020_stsci_fits"] = {
            "file": path.name,
            "format": "atmo2020_stsci_fits",
        }
        path = folder / "bt_settl.xml"
        table = ac.Table()
        table["Wavelength"] = wave * ac.u.um
        table["Flux"] = np.ones(len(wave)) * ac.u.W / ac.u.m**2 / ac.u.um
        table.write(path, format="votable")
        fixtures["bt_settl_svo"] = {
            "file": path.name,
            "format": "bt_settl_svo",
            "wavelength_frame": "air",
        }

    if ac.Dataset is not None and ac.u is not None:
        for flame in (False, True):
            fmt = (
                "sonora_flame_skimmer_v1_r3k" if flame else "sonora_elf_owl_v2"
            )
            name = (
                "spectra_teff_1200_grav_100_mh_0.0_co_1.0_R3000.nc"
                if flame
                else "spectra_logzz_4_teff_1200_grav_100_mh_0.0_co_1.0.nc"
            )
            params = (
                {
                    "effective_temp": {"value": 1200, "unit": "K"},
                    "gravity": {"value": 1e4, "unit": "cm/s2"},
                    "logmh": 0.0,
                    "cto": {"value": 1.0, "unit": "relative solar 0.458"},
                }
                if flame
                else {
                    "teff": {"value": 1200, "unit": "K"},
                    "logg": {"value": 1e4, "unit": "cm/s2"},
                    "mh": 0.0,
                    "cto": 1.0,
                    "logkzz": 4.0,
                }
            )
            with ac.Dataset(folder / name, "w") as dataset:
                dataset.createDimension("wavelength", len(wave))
                wavelength = dataset.createVariable(
                    "wavelength", "f8", ("wavelength",)
                )
                flux = dataset.createVariable(
                    "flux_emission" if flame else "flux", "f8", ("wavelength",)
                )
                wavelength.units = "um"
                flux.units = "erg / (s cm2 cm)"
                wavelength[:] = wave
                flux[:] = 1e7
                dataset.planet_params = json.dumps(params)
                if flame:
                    dataset.rebinned_resolution = "R=3000"
            fixtures[fmt] = {"file": name, "format": fmt}
    return fixtures


class ContrastTests(unittest.TestCase):
    """Verify independent physical identities and invalid-input behavior."""

    def setUp(self):
        """Prepare manufactured constant spectra and a unit response."""
        self.wave = np.array([8.0, 10.0, 12.0, 15.0])
        self.star = ac.Spectrum(self.wave, np.ones(4) * 10)
        self.planet = ac.Spectrum(self.wave, np.ones(4) * 2)
        self.response = ac.Response(
            np.array([9.0, 11.0, 14.0]),
            np.ones(3),
            np.array([9.0, 11.0, 14.0]),
            np.ones(3),
            (9.0, 14.0),
            {"wavelength_frame": "vacuum"},
        )

    def test_analytic_contrast_and_both_radii(self):
        """Check emitting-area scaling independently of quadrature details."""
        for star_radius, planet_radius in ((2, 1), (2, 2), (4, 1)):
            actual = ac.contrast_from_spectra(
                self.star,
                self.planet,
                star_radius,
                planet_radius,
                self.response,
            )
            self.assertAlmostEqual(
                actual, 0.2 * (planet_radius / star_radius) ** 2
            )

    def test_photon_integral(self):
        """Integrate constant F_lambda through a unit photon response."""
        expected = 10 * (14**2 - 9**2) / 2 * 1e-6 / (ac.H * ac.C)
        for rule in ("gauss3", "boole"):
            np.testing.assert_allclose(
                ac.integrate_band(self.star, self.response, rule=rule),
                expected,
                rtol=1e-13,
            )

    def test_energy_integral(self):
        """Energy detection omits the photon-energy factor."""
        self.assertAlmostEqual(
            ac.integrate_band(self.star, self.response, detector="energy"), 50
        )

    def test_polynomial_with_separate_filter_and_sky(self):
        """Check quartic integration against an analytic antiderivative."""
        wave = np.array([9.0, 11.0, 14.0])
        spectrum = ac.Spectrum(wave, 2 + 0.1 * wave)
        response = ac.Response(
            wave,
            0.1 + 0.02 * wave,
            np.array([9.0, 10.0, 12.0, 14.0]),
            0.9 - 0.01 * np.array([9.0, 10.0, 12.0, 14.0]),
            (9.0, 14.0),
            {},
        )
        polynomial = np.polynomial.Polynomial([2, 0.1])
        polynomial *= np.polynomial.Polynomial([0.1, 0.02])
        polynomial *= np.polynomial.Polynomial([0.9, -0.01])
        polynomial *= np.polynomial.Polynomial([0, 1e-6 / (ac.H * ac.C)])
        antiderivative = polynomial.integ()
        expected = antiderivative(14) - antiderivative(9)
        for rule in ("gauss3", "boole"):
            np.testing.assert_allclose(
                ac.integrate_band(spectrum, response, rule=rule),
                expected,
                rtol=1e-13,
            )

    def test_unit_conversions(self):
        """Check F_nu, FLAM, and Eddington conversions against definitions."""
        wave = np.array([10.0, 11.0, 12.0])
        _, converted = ac.to_surface_flambda(
            wave, np.ones(3), "um", "erg_cm2_s_Hz"
        )
        expected = 1e-3 * ac.C / (wave * 1e-6) ** 2 * 1e-6
        np.testing.assert_allclose(converted, expected, rtol=1e-14)
        _, converted = ac.to_surface_flambda(
            wave * 1e4, np.ones(3), "angstrom", "erg_cm2_s_A"
        )
        np.testing.assert_array_equal(converted, np.ones(3) * 10)
        _, converted = ac.to_surface_flambda(
            wave, np.ones(3), "um", "W_m2_um", "eddington_flux"
        )
        np.testing.assert_allclose(converted, np.ones(3) * 4 * math.pi)

    def test_mixed_array_units(self):
        """Check equivalent companion spectra expressed in different units."""
        actual = ac.planet_contrast(
            self.wave,
            self.star.surface_flux,
            self.wave * 1e4,
            self.planet.surface_flux / 10,
            star_radius_m=2,
            planet_radius_m=1,
            response=self.response,
            wavelength_unit="um",
            flux_unit="W_m2_um",
            planet_wavelength_unit="angstrom",
            planet_flux_unit="erg_cm2_s_A",
        )
        self.assertAlmostEqual(actual, 0.05)

    def test_invalid_inputs(self):
        """Reject unsafe ratios, unsupported weighting, and frame conflicts."""
        for value in (0, -1, np.inf, np.nan):
            with self.subTest(radius=value), self.assertRaises(ValueError):
                ac.contrast_from_spectra(
                    self.star, self.planet, 2, value, self.response
                )
        for kwargs in ({"detector": "unknown"}, {"rule": "unknown"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                ac.integrate_band(self.star, self.response, **kwargs)
        star = ac.Spectrum(self.wave, np.ones(4), {"wavelength_frame": "air"})
        with self.assertRaises(ValueError):
            ac.integrate_band(star, self.response)
        with self.assertRaises(ValueError):
            ac.integrate_band(
                ac.Spectrum([10, 11, 12], [1, 1, 1]), self.response
            )
        with self.assertRaises(ValueError):
            ac.Spectrum([9, 10, 10], [1, 1, 1])
        with self.assertRaises(ValueError):
            ac.contrast_from_spectra(
                self.star,
                ac.Spectrum(self.wave, np.zeros(4)),
                2,
                1,
                self.response,
            )

    def test_common_response_scaling(self):
        """A common throughput scale cancels in the contrast ratio."""
        scaled = ac.Response(
            self.response.filter_wave_um,
            np.ones(3) * 0.2,
            self.response.sky_wave_um,
            np.ones(3),
            (9.0, 14.0),
            {},
        )
        self.assertAlmostEqual(
            ac.contrast_from_spectra(self.star, self.planet, 2, 1, scaled),
            0.05,
        )

    def test_native_readers_and_direct_file_api(self):
        """Read manufactured native layouts and verify a file-based ratio."""
        with tempfile.TemporaryDirectory() as temporary:
            fixtures = fixture_models(temporary)
            for fmt, mapping in fixtures.items():
                with self.subTest(format=fmt):
                    spectrum = ac.load_spectrum(mapping, temporary)
                    self.assertGreater(len(spectrum.wavelength_um), 2)
                    self.assertEqual(spectrum.metadata["format"], fmt)
                    self.assertEqual(
                        len(spectrum.metadata["files"][0]["sha256"]), 64
                    )
                    if fmt not in ("sonora_bobcat",):
                        np.testing.assert_allclose(spectrum.surface_flux, 1.0)
            result = ac.contrast_from_models(
                fixtures["ascii"],
                fixtures["ascii"],
                2,
                1,
                self.response,
                base=temporary,
            )
            self.assertAlmostEqual(result, 0.25)

    def test_missing_model_and_changed_response(self):
        """Missing inputs raise errors instead of inserting substitutes."""
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(FileNotFoundError):
                ac.load_spectrum(
                    {"file": "missing", "format": "ascii"}, temporary
                )
            Path(temporary, "Nprime_filter_curve.txt").write_text("changed")
            with self.assertRaises(ValueError):
                ac.load_nprime(temporary)


if __name__ == "__main__":
    unittest.main()
