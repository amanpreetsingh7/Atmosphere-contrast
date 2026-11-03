"""Test the file wrapper with manufactured data, never atmosphere grids."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import example_from_files as example


class FileExampleTests(unittest.TestCase):
    """Exercise the documented file workflow and its input error handling."""

    def arguments(self, output=None):
        """Return absolute fixture paths for any working directory."""
        root = Path(__file__).resolve().parent / "examples" / "manufactured"
        values = [
            "--star",
            str(root / "star.json"),
            "--planet",
            str(root / "planet.json"),
            "--response-config",
            str(root / "response.json"),
            "--star-radius-m",
            "2",
            "--planet-radius-m",
            "1",
        ]
        if output is not None:
            values.extend(["--output", str(output)])
        return values

    def test_print_and_save_expected_contrast(self):
        """Check loading, analytic contrast, and saved/printed agreement."""
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "results" / "contrast.json"
            captured = io.StringIO()
            with contextlib.redirect_stdout(captured):
                self.assertEqual(example.main(self.arguments(output)), 0)
            report = json.loads(captured.getvalue())
            self.assertAlmostEqual(report["planet_star_contrast"], 0.05)
            self.assertEqual(
                report, json.loads(output.read_text(encoding="utf-8"))
            )
            self.assertEqual(report["band_limits_um"], [9.0, 14.0])

    def test_reject_invalid_radius(self):
        """Report a negative radius as a command-line input error."""
        arguments = self.arguments()
        arguments[arguments.index("--planet-radius-m") + 1] = "-1"
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                example.main(arguments)
        self.assertEqual(error.exception.code, 2)

    def test_reject_non_object_configuration(self):
        """Require mappings instead of silently accepting a JSON list."""
        with tempfile.TemporaryDirectory() as temporary:
            filename = Path(temporary) / "bad.json"
            filename.write_text("[]", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "JSON object"):
                example.read_mapping(filename)


if __name__ == "__main__":
    unittest.main()
