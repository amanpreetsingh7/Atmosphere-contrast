# Contributing

Open an issue describing the scientific or software change before a large
refactor. For a small fix, a pull request with the reason and relevant checks
is sufficient. For pyNOMIC integration, coordinate the destination interface
and attribution with its maintainers; this repository does not modify it.

## Local checks

Run the core tests and examples with Python 3.9 or newer:

```sh
python -m pip install -r requirements_native.txt
python -m unittest discover -v
python example_usage.py
```

Use Python 3.12 for the pinned development formatter:

```sh
python -m pip install -r requirements_dev.txt
python -m black --check --line-length 79 .
python -m pycodestyle --max-line-length=79 .
```

Keep imports at module scope, document classes/functions and physical units,
and explain non-obvious numerical steps in comments. Avoid nested functions
unless a clear reason justifies one. Keep examples explicit about whether
inputs are manufactured or native scientific spectra.

For numerical changes, test a physical or mathematical identity, not just a
copy of the implementation. Check independent quadrature, known radius
scaling, unit conversions, missing coverage, and wavelength conventions.
Do not weaken a source-data hash check without reviewing the new input.

Do not commit downloaded model grids, machine-specific configurations,
generated results, credentials, virtual environments, or instrument datasets.
For new native readers, document the exact release layout, unit authority,
and parameter provenance. Manufactured layout fixtures are necessary but
do not establish physical accuracy across the full model grid.

