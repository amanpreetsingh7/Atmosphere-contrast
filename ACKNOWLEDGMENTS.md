# Acknowledgments and source data

Developed by Amanpreet Singh from the contrast calculation in Atmosphere Lab.
The scientific goals were discussed with Kevin Wagner. Vivek Vijayakumar
provided guidance on the contrast interface, filter/sky transmission inputs,
and preparation of the code for review and possible use with pyNOMIC.
These acknowledgments do not imply an upstream merge or endorsement.

## Transmission curves

The optional Nprime loader expects external files from the reviewed pyNOMIC
commit `e6c63c740a6fe1ca74e57228fafa3d7da8b6e67a`:

- `Nprime_filter_curve.txt`
- `MtGrahamTransmittance.npz`

Source: <https://github.com/vivekv7323/pyNOMIC/tree/e6c63c740a6fe1ca74e57228fafa3d7da8b6e67a/filterdata>

The curves are not redistributed here. The loader checks their hashes to
identify the exact reviewed inputs. Reference sky transmission is not a
measurement of conditions at an individual observation.

## Atmosphere grids

Native spectra are supplied separately by the user. Cite the specific model
family, release, and data source used for each scientific result. Reader
support is not a license to redistribute source archives. Record the native
headers and adopted wavelength/flux conventions when reporting results.

## Software and license

The Python implementation uses NumPy, with optional Astropy and netCDF4
readers. Cite these packages as appropriate to the work performed. See their
projects for their own licenses and citation guidance.

The repository's MIT license applies to its code and documentation. It does
not relicense external atmosphere spectra, transmission curves, or dependencies.

