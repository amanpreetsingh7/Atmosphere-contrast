"""Run a self-contained manufactured example of the contrast interface.

Run ``python example_usage.py``. No native model files are downloaded or
used; the constant spectra make the expected answer exactly calculable.
See README.md for the actual atmosphere-file workflow.
"""

import numpy as np

from atmosphere_contrast import Response, planet_contrast


def main():
    """Calculate the analytic example and print its dimensionless contrast."""
    wavelengths_um = np.array([8.0, 10.0, 12.0, 15.0])
    band_wavelengths_um = np.array([9.0, 11.0, 14.0])

    # This is a unit-transmission test response, not the real Nprime band.
    response = Response(
        filter_wave_um=band_wavelengths_um,
        filter_transmission=np.ones(3),
        sky_wave_um=band_wavelengths_um,
        sky_transmission=np.ones(3),
        limits_um=(9.0, 14.0),
        metadata={"wavelength_frame": "vacuum"},
    )

    # Constant planet/star surface flux ratio is 2/10. With Rp/Rs=1/2,
    # the band contrast must be (1/2)^2 * (2/10) = 0.05.
    contrast = planet_contrast(
        wavelengths_um,
        np.full(4, 10.0),
        wavelengths_um,
        np.full(4, 2.0),
        star_radius_m=2.0,
        planet_radius_m=1.0,
        response=response,
        wavelength_unit="um",
        flux_unit="W_m2_um",
    )
    print(f"Manufactured planet/star contrast: {contrast:.12g}")
    print("Expected: 0.05. No native atmosphere model was calculated.")


if __name__ == "__main__":
    main()
