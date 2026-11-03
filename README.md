# Atmosphere Contrast

Calculate a **thermal planet-to-star flux ratio** from stellar and planetary
atmosphere spectra, emitting radii, and filter/sky transmission curves.
The output is a single dimensionless contrast for a downstream source
injector that uses the same planet/star brightness convention.

Developed by **Amanpreet Singh**. Initial review module: **0.1.0**.

This repository separates the brightness calculation from Atmosphere Lab.
It contains the documented calculation, native readers, examples, and tests.
It has no Streamlit dependency. The full app, plotting, stellar heating,
orbit dynamics, model interpolation, evolution, PSFs, sensitivity estimates,
and injection/recovery are outside this extraction's scope.

## Quick start: no model downloads

Use Python 3.9 or newer. Open a terminal in this repository folder:

```sh
python -m pip install -r requirements.txt
python example_usage.py
python -m unittest discover -v
```

The example prints a contrast of **0.05**. It uses deliberately constant
surface spectra and a unit-transmission test band, not real atmospheres.
Tests also use manufactured inputs. Optional native-layout tests are
included only when their required packages are installed.

To run the file-based interface with bundled manufactured inputs:

```sh
python example_from_files.py \
  --star examples/manufactured/star.json \
  --planet examples/manufactured/planet.json \
  --response-config examples/manufactured/response.json \
  --star-radius-m 2 \
  --planet-radius-m 1
```

This also returns **0.05**. The arbitrary radii are a mathematical check,
not an astronomical system. Add `--output results/contrast.json` to save
a run; use `--help` for arguments.

## Use in another program

```python
from atmosphere_contrast import planet_contrast

contrast = planet_contrast(
    star_wavelength_um,
    star_surface_flux,
    planet_wavelength_um,
    planet_surface_flux,
    star_radius_m=star_radius_m,
    planet_radius_m=planet_radius_m,
    response=response,
    wavelength_unit="um",
    flux_unit="W_m2_um",
    wavelength_frame="vacuum",
)
```

The result is planet flux divided by star flux, **not a percent, magnitude,
or detection limit**. Inputs must be emergent surface fluxes, not observed
fluxes at Earth or specific intensities. `Response` stores separate filter
and sky curves. Use `contrast_from_models` for native file mappings or
`contrast_from_spectra` for loaded spectra. See the
[API and native-format guide](docs/API_AND_MODELS.md).

## Actual atmosphere files

```sh
python -m pip install -r requirements_native.txt
```

Model files and response tables remain external; nothing is downloaded
implicitly. The [native-model walkthrough](docs/NATIVE_EXAMPLE.md) explains
the CK04/Bobcat configurations in `examples/native/` and the pinned Nprime
response. These template configurations need existing matching files.
For other bands, supply explicit `Response` arrays.

## Physics and conventions

For photon weighting:

$$
C = \left(\frac{R_p}{R_\star}\right)^2
\frac{\int F_{\lambda,p}(\lambda)T_f(\lambda)T_s(\lambda)
\lambda/(hc)\,d\lambda}
{\int F_{\lambda,\star}(\lambda)T_f(\lambda)T_s(\lambda)
\lambda/(hc)\,d\lambda}.
$$

Here `Rp` and `Rstar` are emitting radii, `F_lambda` is emergent surface
spectral flux density, `Tf` is filter transmission, and `Ts` is sky
transmission. `h` and `c` are Planck's constant and the speed of light.
Energy weighting omits the photon factor `lambda/(h*c)`.

Internal wavelengths are in micrometres; surface F_lambda is in
W m^-2 um^-1. Only the photon-energy factor uses wavelength in metres.
Emergent surface flux needs no extra pi. Explicit Eddington H is converted
with 4*pi. Common observer distance and common response scale cancel.

The integrator retains independent spectrum/filter/sky knots and interpolates
each input linearly. Three-point Gaussian quadrature integrates that
tabulated representation; independent Boole quadrature is available.
There is no extrapolation or silent partial-band integration.

### Interpretation limits

- Radii are imposed independently; mass/age/gravity consistency is not inferred.
- Models are selected explicitly; no parameter-grid interpolation is performed.
- This extraction does not calculate reflected light, heating, or thermal lag.
- Wavelength frames must agree; air/vacuum conversion is not automatic.
- Reference sky transmission does not establish individual observing conditions
  or include sky emission. Full telescope optics and detector QE are not inferred.
- Brightness is not an instrument sensitivity or a detection probability.

## Repository contents

| File | Purpose |
| --- | --- |
| `atmosphere_contrast.py` | Contrast calculation and native readers |
| `example_usage.py` | Analytic array example |
| `example_from_files.py` | Explicit model/response JSON command-line interface |
| `examples/manufactured/` | Synthetic files; expected contrast 0.05 |
| `examples/native/` | CK04/Bobcat local-file templates |
| `test_contrast.py` | Original analytic and reader checks |
| `test_file_example.py` | File workflow and input-error checks |
| `docs/API_AND_MODELS.md` | Full interface, units, formats, and provenance |
| `docs/NATIVE_EXAMPLE.md` | Native-model setup and command-line example |
| `VALIDATION.json` | Original module preparation record |
| `REPOSITORY_VALIDATION.json` | Checks for this repository package |
| `.github/workflows/tests.yml` | Configured automated tests and style checks |
| `CONTRIBUTING.md` | Development and review instructions |
| `ACKNOWLEDGMENTS.md` | Guidance and external-source attribution |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | MIT license for code/documentation |
| `UPLOAD_TO_GITHUB.md` | Website and Desktop upload instructions |

## Validation status

The core module is byte-for-byte identical to the prepared review ZIP,
extracted from the supplied `0.8.5-single-file.1` standalone export. Original
checks included 128 manufactured scalar regression cases with identical
results, analytic/unit tests, and 11 manufactured native layouts with
optional dependencies installed. See `VALIDATION.json` for the exact scope.

Repository preparation adds a file example, tests, documentation, and GitHub
settings. Local checks use Python 3.12 and check Python 3.9 syntax. The
workflow is configured to test Python 3.9 and 3.12; GitHub execution is only
established once it runs there. No real atmosphere archive was recomputed
during this preparation. Numerical tests do not establish physical accuracy
across full grids. See `REPOSITORY_VALIDATION.json` for package checks.

For formatting, use Python 3.12 and the pinned tools:

```sh
python -m pip install -r requirements_dev.txt
python -m black --check --line-length 79 .
python -m pycodestyle --max-line-length=79 .
```

## License, citation, and integration

Code and documentation use the [MIT license](LICENSE). External atmosphere
spectra, transmission tables, and dependencies retain their own terms.
Cite the model grids and source data used for scientific results as well
as this software. See [acknowledgments](ACKNOWLEDGMENTS.md) and `CITATION.cff`.

This repository does not modify or merge into pyNOMIC. Its scalar output is
intended for a separately reviewed integration. The installed Atmosphere Lab
app and its downloaded models remain separate.

