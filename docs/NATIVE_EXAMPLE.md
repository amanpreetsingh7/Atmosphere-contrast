# Run with native atmosphere spectra

This uses explicit CK04 stellar and Sonora Bobcat planetary files plus the
reviewed Nprime/sky response. Native model files are external. This recipe
is not a claim that actual atmospheres were recomputed during repository
preparation; run it with your existing matching model files.

## 1. Install readers

From the repository folder, in your chosen Python environment:

```sh
python -m pip install -r requirements_native.txt
python -m unittest discover -v
```

## 2. Configure existing spectra

Edit `examples/native/star_ck04.json` to point to the CK04 STScI FITS file
for a 5000 K solar-metallicity star and its `g45` gravity column. The
template filename is `ckp00_5000.fits`. The declarations, gravity column,
and actual file/header must match; changing a label does not change a spectrum.

Edit `examples/native/planet_bobcat.json` to point to an original Bobcat
file such as `sp_t1200g100nc_m0.0.gz`. The loader converts its native cgs
per-frequency flux to surface F_lambda. Inspect the header rather than
assuming rounded filename coordinates are exact physical parameters.

Template paths expect files under the repository's ignored `models/`
folder. Absolute paths to existing external storage also work. Each
relative path is resolved from the JSON file's parent directory, not the
launch directory. No files need to be moved out of the installed app.

Other grids require their exact documented release layout; see
[API_AND_MODELS.md](API_AND_MODELS.md). The newer app's BT-Settl STScI
FITS adapter is not part of this extraction and cannot be relabeled as SVO.

## 3. Supply the reviewed response

Pinned source directory:
<https://github.com/vivekv7323/pyNOMIC/tree/e6c63c740a6fe1ca74e57228fafa3d7da8b6e67a/filterdata>

Place these two original files together in any local directory:

- `Nprime_filter_curve.txt`
- `MtGrahamTransmittance.npz`

They can remain in existing app storage or the ignored `filterdata/`
directory. The loader checks these SHA-256 values:

```text
Nprime_filter_curve.txt
8039385921be84c56bcbaf82c4f7174d3cd0dd0f323dfa94eb187ab8e7b6c9d6

MtGrahamTransmittance.npz
36c83823d5abe2267e0b07b3c7b4336a4405c6acc4ffeb9b05d391240a1aeda2
```

If a file differs, obtain the reviewed version or review the new data.
Do not disable the check just to obtain an answer. A GitHub HTML page is
not a transmission table. No observation dataset or injection test is needed.

## 4. Calculate

With actual model paths configured and response files in `filterdata/`:

```sh
python example_from_files.py \
  --star examples/native/star_ck04.json \
  --planet examples/native/planet_bobcat.json \
  --nprime-dir filterdata \
  --star-radius-m 695700000 \
  --planet-radius-m 71492000 \
  --output results/native_contrast.json
```

The radii are one nominal solar radius and one nominal Jupiter equatorial
radius. They are explicit emitting-area assumptions, not radii derived from
atmospheric gravity or evolutionary mass/age models.

The output contains the dimensionless `planet_star_contrast`, detector
weighting, radii, band limits, module version, and configuration/response
paths. Keep inputs and inspect `Spectrum.metadata` separately for original
spectral hashes and verified model parameters when reporting scientific work.

## 5. Use the result

Call `contrast_from_models`, `contrast_from_spectra`, or `planet_contrast`
from another program. Supply its scalar according to the downstream
injector's brightness convention. This module does not fit a PSF, patch
pyNOMIC, or run an injection.

Both spectra must cover the whole band and match the response wavelength
frame. Incompatible units, missing coverage, and mismatched frames raise
errors rather than selecting substitute models.

For observational validation, choose compatible parameters, record source
grids and conditions, and compare with independently calculated fluxes or
matched-band measurements. Manufactured tests alone do not establish the
model physics or validate the selected atmosphere for a particular target.

