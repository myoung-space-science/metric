import pytest

from metric._reference import _UNITS


@pytest.fixture
def named_units():
    """All unit names and symbols."""
    keys = 'symbol', 'name'
    return tuple({k: unit[k] for k in keys} for unit in _UNITS)

