import re
from epymetheus import epymetheus


def test_version_exists():
    assert hasattr(epymetheus, "__version__")


def test_version_is_semver():
    assert re.match(r"^\d+\.\d+\.\d+", epymetheus.__version__)
