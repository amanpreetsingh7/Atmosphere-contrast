# Changelog

## 0.1.0 - initial review module

- Extract the thermal contrast calculation and native readers from the supplied
  `0.8.5-single-file.1` standalone export of Atmosphere Lab.
- Expose array, loaded-spectrum, and model-file interfaces returning a scalar
  dimensionless planet/star flux ratio.
- Keep explicit surface-flux units, wavelength frames, full-band coverage,
  response provenance, and emitting-radius scaling.
- Move imports and functions to module scope and add docstrings and comments.
- Include analytic/manufactured tests and the original validation record.
- Add repository documentation, file-based examples, an MIT license,
  contribution instructions, and a GitHub Actions workflow.

The calculation module is byte-for-byte identical to the previously prepared
review ZIP. Repository preparation does not update the installed app or add
heating, orbit dynamics, interpolation, plotting, or injection/recovery.

