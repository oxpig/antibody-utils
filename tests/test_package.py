import re

import antibody_utils


def test_version():
    assert re.fullmatch(r"\d+\.\d+\.\d+", antibody_utils.__version__)
