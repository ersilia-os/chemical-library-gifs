import os
import tempfile
import pytest

from chemgifs.main import run

SAMPLE_CSV = os.path.join(os.path.dirname(__file__), "sample.csv")


def test_run_produces_gif():
    with tempfile.TemporaryDirectory() as tmp:
        output_gif = os.path.join(tmp, "output.gif")
        run(
            input_csv=SAMPLE_CSV,
            output_file=output_gif,
            color_name="white",
            size=256,
            duration_ms=200,
        )
        assert os.path.exists(output_gif)
        assert os.path.getsize(output_gif) > 0
