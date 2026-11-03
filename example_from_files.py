"""Calculate one contrast from explicit model and response configuration files.

Run with ``--help`` for arguments. Relative model filenames are resolved
from the directory containing their JSON configuration. No models or
transmission tables are downloaded or silently substituted.
"""

import argparse
import json
from pathlib import Path

from atmosphere_contrast import (
    Response,
    __version__,
    contrast_from_spectra,
    load_nprime,
    load_spectrum,
)


def read_mapping(filename):
    """Read a JSON object and return it with its absolute configuration path.

    Parameters
    ----------
    filename : str or pathlib.Path
        UTF-8 JSON configuration file.

    Returns
    -------
    tuple of dict and pathlib.Path
        The mapping and its resolved path. Model paths remain relative
        to the returned path's parent directory.

    Raises
    ------
    OSError, ValueError
        If the file cannot be read or does not contain a JSON object.
    """
    path = Path(filename).expanduser().resolve()
    with path.open(encoding="utf-8") as stream:
        mapping = json.load(stream)
    if not isinstance(mapping, dict):
        raise ValueError(f"Expected a JSON object in {path.name}.")
    return mapping, path


def calculate(args):
    """Load declared inputs and return the scalar contrast and input settings.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed paths, emitting radii in metres, and detector weighting.

    Returns
    -------
    dict
        A dimensionless planet/star contrast and reproducibility settings.
        Source paths describe this run; no source spectra are copied.

    Notes
    -----
    This wrapper makes no catalog selections, heating adjustments, or
    evolutionary inferences. The core validates units and band coverage.
    """
    star_mapping, star_path = read_mapping(args.star)
    planet_mapping, planet_path = read_mapping(args.planet)

    # Resolve each model relative to its own configuration, so stellar
    # and planetary files can live in different storage locations.
    star = load_spectrum(star_mapping, star_path.parent)
    planet = load_spectrum(planet_mapping, planet_path.parent)

    # Nprime uses the reviewed external tables. A custom response instead
    # declares its own wavelength samples and transmissions explicitly.
    if args.nprime_dir is not None:
        response = load_nprime(Path(args.nprime_dir).expanduser())
        response_source = str(Path(args.nprime_dir).expanduser().resolve())
    else:
        response_mapping, response_path = read_mapping(args.response_config)
        response = Response(**response_mapping)
        response_source = str(response_path)

    # Radii scale the emitting areas; the output is a ratio, not percent.
    value = contrast_from_spectra(
        star,
        planet,
        args.star_radius_m,
        args.planet_radius_m,
        response,
        detector=args.detector,
    )
    return {
        "module_version": __version__,
        "planet_star_contrast": value,
        "contrast_unit": "dimensionless flux ratio",
        "detector_weighting": args.detector,
        "star_radius_m": args.star_radius_m,
        "planet_radius_m": args.planet_radius_m,
        "star_configuration": str(star_path),
        "planet_configuration": str(planet_path),
        "response_source": response_source,
        "band_limits_um": [float(v) for v in response.limits_um],
    }


def build_parser():
    """Construct the command-line parser with explicit required inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", required=True, help="Stellar model JSON.")
    parser.add_argument("--planet", required=True, help="Planet model JSON.")
    response_group = parser.add_mutually_exclusive_group(required=True)
    response_group.add_argument(
        "--nprime-dir", help="Directory with the pinned Nprime/sky tables."
    )
    response_group.add_argument(
        "--response-config", help="JSON containing explicit Response arrays."
    )
    parser.add_argument(
        "--star-radius-m", type=float, required=True, help="Star radius [m]."
    )
    parser.add_argument(
        "--planet-radius-m",
        type=float,
        required=True,
        help="Planet radius [m].",
    )
    parser.add_argument(
        "--detector", choices=("photon", "energy"), default="photon"
    )
    parser.add_argument(
        "--output", help="Optional JSON output path; printed to stdout too."
    )
    return parser


def main(argv=None):
    """Parse arguments, calculate a contrast, and print or save the report.

    Parameters
    ----------
    argv : list of str or None, optional
        Arguments for programmatic use; None uses the command line.

    Returns
    -------
    int
        Zero on success. Invalid inputs exit through argparse with code 2.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = calculate(args)
        text = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output is not None:
            destination = Path(args.output).expanduser()
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(text, encoding="utf-8")
    except (OSError, ValueError, TypeError, ImportError) as error:
        parser.error(str(error))
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
