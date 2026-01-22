"""
All independently defined variables for the `metric` namespace.
"""

import collections.abc
import typing

import aliasedkeys


T = typing.TypeVar('T')


_PREFIXES = [
    {'symbol': 'Y', 'name': 'yotta', 'factor': 1e+24},
    {'symbol': 'Z', 'name': 'zetta', 'factor': 1e+21},
    {'symbol': 'E', 'name': 'exa', 'factor': 1e+18},
    {'symbol': 'P', 'name': 'peta', 'factor': 1e+15},
    {'symbol': 'T', 'name': 'tera', 'factor': 1e+12},
    {'symbol': 'G', 'name': 'giga', 'factor': 1e+9},
    {'symbol': 'M', 'name': 'mega', 'factor': 1e+6},
    {'symbol': 'k', 'name': 'kilo', 'factor': 1e+3},
    {'symbol': 'h', 'name': 'hecto', 'factor': 1e+2},
    {'symbol': 'da', 'name': 'deca', 'factor': 1e+1},
    {'symbol': '', 'name': '', 'factor': 1e0},
    {'symbol': ' ', 'name': None, 'factor': 0.0},
    {'symbol': 'd', 'name': 'deci', 'factor': 1e-1},
    {'symbol': 'c', 'name': 'centi', 'factor': 1e-2},
    {'symbol': 'm', 'name': 'milli', 'factor': 1e-3},
    {'symbol': 'μ', 'name': 'micro', 'factor': 1e-6},
    {'symbol': 'n', 'name': 'nano', 'factor': 1e-9},
    {'symbol': 'p', 'name': 'pico', 'factor': 1e-12},
    {'symbol': 'f', 'name': 'femto', 'factor': 1e-15},
    {'symbol': 'a', 'name': 'atto', 'factor': 1e-18},
    {'symbol': 'z', 'name': 'zepto', 'factor': 1e-21},
    {'symbol': 'y', 'name': 'yocto', 'factor': 1e-24},
]


UNITY = {'#', '1'}
"""Strings that represent dimensionless units."""


SYSTEMS = {'mks', 'cgs'}
"""The metric systems known to this module."""


_UNITS = [
    {
        'symbol': 'm',
        'quantity': 'length',
        'name': 'meter',
    },
    {
        'symbol': 'au',
        'quantity': 'length',
        'name': 'astronomical unit',
    },
    {
        'symbol': 'Rs',
        'quantity': 'length',
        'name': 'solar radius',
        'plural': 'solar radii',
    },
    {
        'symbol': 'g',
        'quantity': 'mass',
        'name': 'gram',
    },
    {
        'symbol': 'nuc',
        'quantity': 'mass number',
        'name': 'nucleon',
    },
    {
        'symbol': 'amu',
        'quantity': 'mass number',
        'name': 'atomic mass unit',
    },
    {
        'symbol': 's',
        'quantity': 'time',
        'name': 'second',
    },
    {
        'symbol': 'min',
        'quantity': 'time',
        'name': 'minute',
    },
    {
        'symbol': 'h',
        'quantity': 'time',
        'name': 'hour',
    },
    {
        'symbol': 'd',
        'quantity': 'time',
        'name': 'day',
    },
    {
        'symbol': 'A',
        'quantity': 'current',
        'name': 'ampere',
    },
    {
        'symbol': 'K',
        'quantity': 'temperature',
        'name': 'kelvin',
    },
    {
        'symbol': 'mol',
        'quantity': 'amount',
        'name': 'mole',
    },
    {
        'symbol': '#',
        'quantity': 'number',
        'name': 'count',
    },
    {
        'symbol': 'cd',
        'quantity': 'luminous intensity',
        'name': 'candela',
    },
    {
        'symbol': 'rad',
        'quantity': 'plane angle',
        'name': 'radian',
    },
    {
        'symbol': 'deg',
        'quantity': 'plane angle',
        'name': 'degree',
    },
    {
        'symbol': 'sr',
        'quantity': 'solid angle',
        'name': 'steradian',
    },
    {
        'symbol': 'Hz',
        'quantity': 'frequency',
        'name': 'hertz',
    },
    {
        'symbol': 'J',
        'quantity': 'energy',
        'name': 'joule',
    },
    {
        'symbol': 'erg',
        'quantity': 'energy',
        'name': 'erg',
    },
    {
        'symbol': 'eV',
        'quantity': 'energy',
        'name': 'electronvolt',
    },
    {
        'symbol': 'N',
        'quantity': 'force',
        'name': 'newton',
    },
    {
        'symbol': 'dyn',
        'quantity': 'force',
        'name': 'dyne',
    },
    {
        'symbol': 'Pa',
        'quantity': 'pressure',
        'name': 'pascal',
    },
    {
        'symbol': 'W',
        'quantity': 'power',
        'name': 'watt',
    },
    {
        'symbol': 'C',
        'quantity': 'charge',
        'name': 'coulomb',
    },
    {
        'symbol': 'statC',
        'quantity': 'charge',
        'name': 'statcoulomb',
    },
    {
        'symbol': 'statA',
        'quantity': 'current',
        'name': 'statampere',
    },
    {
        'symbol': 'statV',
        'quantity': 'potential',
        'name': 'statvolt',
    },
    {
        'symbol': 'e',
        'quantity': 'charge',
        'name': 'fundamental charge',
    },
    {
        'symbol': 'V',
        'quantity': 'potential',
        'name': 'volt',
    },
    {
        'symbol': 'Ω',
        'quantity': 'resistance',
        'name': 'ohm',
    },
    {
        'symbol': 'S',
        'quantity': 'conductance',
        'name': 'seimens',
    },
    {
        'symbol': 'F',
        'quantity': 'capacitance',
        'name': 'farad',
    },
    {
        'symbol': 'Wb',
        'quantity': 'magnetic flux',
        'name': 'weber',
    },
    {
        'symbol': 'Mx',
        'quantity': 'magnetic flux',
        'name': 'maxwell',
    },
    {
        'symbol': 'Oe',
        'quantity': 'magnetic intensity',
        'name': 'Oersted',
    },
    {
        'symbol': 'H',
        'quantity': 'inductance',
        'name': 'henry',
    },
    {
        'symbol': 'T',
        'quantity': 'magnetic induction',
        'name': 'tesla',
    },
    {
        'symbol': 'G',
        'quantity': 'magnetic induction',
        'name': 'gauss',
    },
    {
        'symbol': 'lm',
        'quantity': 'luminous flux',
        'name': 'lumen',
    },
    {
        'symbol': 'lx',
        'quantity': 'illuminance',
        'name': 'lux',
        'plural': 'luxes',
    },
    {
        'symbol': 'Bq',
        'quantity': 'radioactivity',
        'name': 'becquerel',
    },
    {
        'symbol': 'Ci',
        'quantity': 'radioactivity',
        'name': 'Curie',
    },
    {
        'symbol': 'Gy',
        'quantity': 'dosage',
        'name': 'gray',
    },
    {
        'symbol': 'P',
        'quantity': 'viscosity',
        'name': 'poise',
    },
    {
        'symbol': '1',
        'quantity': 'identity',
        'name': 'unitless',
        'plural': 'unitless',
    },
]


def _build_named_units(
    prefixes: typing.List[typing.Dict[T, typing.Any]],
    units: typing.List[typing.Dict[T, str]],
) -> aliasedkeys.Mapping[T, typing.Dict[T, typing.Any]]:
    """Helper for building the collection of named units."""
    # Find the null metric prefix (0).
    null_prefix = next(
        prefix for prefix in prefixes if prefix['factor'] == 0.0
    )
    # Find the metric prefix corresponding to unity (1).
    base_prefix = next(
        prefix for prefix in prefixes if prefix['factor'] == 1.0
    )
    # Find the unitless quantity.
    unitless = next(
        unit for unit in units if unit['quantity'] == 'identity'
    )
    # Build mapping from (prefix, unit) combinations, as appropriate.
    mapped = {}
    for unit in units:
        if unit is unitless:
            # Prevent meaningless unitless expressions (e.g., cenitunitless).
            mapped[(unitless['symbol'], unitless['name'])] = {
                'prefix': base_prefix, 'base': unitless
            }
        else:
            for prefix in prefixes:
                if prefix is base_prefix:
                    # Add the unscaled unit.
                    mapped[(unit['symbol'], unit['name'], unit['plural'])] = {
                        'prefix': base_prefix, 'base': unit
                    }
                elif prefix is not null_prefix:
                    # Add the full (prefix, unit) combination.
                    key = [
                        f"{prefix['symbol']}{unit['symbol']}",
                        f"{prefix['name']}{unit['name']}",
                        f"{prefix['name']}{unit['plural']}",
                    ]
                    if prefix['symbol'] == 'μ':
                        key += [f"u{unit['symbol']}"]
                    mapped[tuple(key)] = {'prefix': prefix, 'base': unit}
    return aliasedkeys.Mapping(mapped)


def _add_plurals(
    collection: list[dict[T, str]],
    name_key: T='name',
    plural_key: T='plural',
) -> list[dict[T, str]]:
    """Insert trivial plurals."""
    # NOTE: This implementation takes advantage of the fact that an existing
    # value of `entry[plural_key]` will overwrite the inserted value during
    # unpacking. It is equivalent to
    # - result = []
    # - for entry in collection:
    # -     if plural_key not in entry:
    # -         entry[plural_key] = f"{entry[name_key]}s"
    # -     result.append(entry)
    # - return result
    return [
        {plural_key: f"{entry[name_key]}s", **entry}
        for entry in collection
    ]


NAMED_UNITS = _build_named_units(_PREFIXES, _add_plurals(_UNITS))


# A note about angles: Kalinin (2019) "On the status of plane and solid in the
# International System of Units (SI)" makes a very compelling argument that the
# plane- and solid-angle units should not be '1'; instead:
# - the plane-angle unit should be 'radian' ('rad'),
# - the solid-angle unit should be 'steradian' ('sr'), and
# - the plane angle should be a base quantity,
# - 'radian' should be considered a base unit,
# - 'steradian' should be considered a derived unit, with 'sr = rad^2'.

# References and notes on quantities, dimensions, and units:
# - https://en.wikipedia.org/wiki/International_System_of_Quantities#Base_quantities
# - https://www.nist.gov/pml/weights-and-measures/metric-si/si-units
# - This module uses 'H' to represent the dimension of temperature because the
#   SI character, $\Theta$, is not an ASCII character. I also considered 'O'
#   because of its similarity to $\Theta$, but it looks too much like '0'
#   (zero), which, ironically, looks more like $\Theta$ than 'O'.
# - This module adds an identity quantity with '1' as its dimension and unit

BASE_QUANTITIES = [
    {'name': 'identity', 'dimension': '1', 'unit': '1'},
    {'name': 'amount', 'dimension': 'N', 'unit': 'mol'},
    {'name': 'current', 'dimension': 'I', 'unit': 'A'},
    {'name': 'length', 'dimension': 'L', 'unit': 'm'},
    {'name': 'luminous intensity', 'dimension': 'J', 'unit': 'cd'},
    {'name': 'mass', 'dimension': 'M', 'unit': 'kg'},
    {'name': 'temperature', 'dimension': 'H', 'unit': 'K'},
    {'name': 'time', 'dimension': 'T', 'unit': 's'},
]


QUANTITIES = {
    'acceleration': 'velocity / time',
    'amount': {
        'dimensions': {
            'mks': 'N',
            'cgs': 'N',
        },
        'units': {
            'mks': 'mol',
            'cgs': 'mol',
        },
    },
    'area': 'length^2',
    'capacitance': {
        'dimensions': {
            'mks': '(T^2 * I)^2 / (M * L^2)',
            'cgs': 'L',
        },
        'units': {
            'mks': 'F',
            'cgs': 'cm',
        },
    },
    'charge': {
        'dimensions': {
            'mks': 'I * T',
            'cgs': '(M^1/2 * L^3/2) / T',
        },
        'units': {
            'mks': 'C',
            'cgs': 'statC',
        },
    },
    'charge density': 'charge / volume',
    'conductance': {
        'dimensions': {
            'mks': '(T^3 * I^2) / (M * L^2)',
            'cgs': 'L / T',
        },
        'units': {
            'mks': 'S',
            'cgs': 'cm / s',
        },
    },
    'conductivity': 'conductance / length',
    'current': {
        'dimensions': {
            'mks': 'I',
            'cgs': '(M^1/2 * L^3/2) / T^2',
        },
        'units': {
            'mks': 'A',
            'cgs': 'statA',
        },
    },
    'current density': 'current / area',
    'displacement': {
        'dimensions': {
            'mks': 'I * T / L^2',
            'cgs': 'M^1/2 / (L^1/2 * T)',
        },
        'units': {
            'mks': 'C / m^2',
            'cgs': 'statC / m^2',
        },
    },
    'dosage': {
        'dimensions': {
            'mks': 'L^2 / T^2',
            'cgs': 'L^2 / T^2',
        },
        'units': {
            'mks': 'Gy',
            'cgs': 'erg / g', # Historically, 'rad', but that's in use.
        },
    },
    'electric charge': 'charge',
    'electric field': 'potential / length',
    'electromotance': 'potential',
    'energy': {
        'dimensions': {
            'mks': '(M * L^2) / T^2',
            'cgs': '(M * L^2) / T^2',
        },
        'units': {
            'mks': 'J',
            'cgs': 'erg',
        },
    },
    'energy density': 'energy / volume',
    'fluence': 'particle fluence',
    'flux': 'particle flux',
    'force': {
        'dimensions': {
            'mks': '(M * L) / T^2',
            'cgs': '(M * L) / T^2',
        },
        'units': {
            'mks': 'N',
            'cgs': 'dyn',
        },
    },
    'frequency': {
        'dimensions': {
            'mks': '1 / T',
            'cgs': '1 / T',
        },
        'units': {
            'mks': 'Hz',
            'cgs': 'Hz',
        },
    },
    'identity': {
        'dimensions': {
            'mks': '1',
            'cgs': '1',
        },
        'units': {
            'mks': '1',
            'cgs': '1',
        },
    },
    'illumunance': { # See note about radian (Kalinin 2019).
        'dimensions': {
            'mks': 'J / L^2',
            'cgs': 'J / L^2',
        },
        'units': {
            'mks': 'cd * sr / m^2',
            'cgs': 'cd * sr / cm^2',
        },
    },
    'impedance': {
        'dimensions': {
            'mks': '(M * L^2) / (T^3 * I)',
            'cgs': 'T / L',
        },
        'units': {
            'mks': 'ohm',
            'cgs': 's / cm',
        },
    },
    'inductance': {
        'dimensions': {
            'mks': '(M * L^2) / (I * T)^2',
            'cgs': 'T^2 / L',
        },
        'units': {
            'mks': 'H',
            'cgs': 's^2 / cm',
        },
    },
    'induction': 'magnetic induction',
    'integral flux': 'flux * energy',
    'length': {
        'dimensions': {
            'mks': 'L',
            'cgs': 'L',
        },
        'units': {
            'mks': 'm',
            'cgs': 'cm',
        },
    },
    'luminous flux': { # See note about radian (Kalinin 2019).
        'dimensions': {
            'mks': 'J',
            'cgs': 'J',
        },
        'units': {
            'mks': 'cd * sr',
            'cgs': 'cd * sr',
        },
    },
    'luminous intensity': {
        'dimensions': {
            'mks': 'J',
            'cgs': 'J',
        },
        'units': {
            'mks': 'cd',
            'cgs': 'cd',
        },
    },
    'magnetic field': 'magnetic induction',
    'magnetic flux': {
        'dimensions': {
            'mks': '(M * L^2) / (T^2 * I)',
            'cgs': '(M^1/2 * L^3/2) / T',
        },
        'units': {
            'mks': 'Wb',
            'cgs': 'Mx',
        },
    },
    'magnetic induction': {
        'dimensions': {
            'mks': 'M / (T^2 * I)',
            'cgs': 'M^1/2 / (L^1/2 * T)',
        },
        'units': {
            'mks': 'T',
            'cgs': 'G',
        },
    },
    'magnetic intensity': {
        'dimensions': {
            'mks': 'I / L',
            'cgs': 'M^1/2 / (L^1/2 * T)',
        },
        'units': {
            'mks': 'A / m',
            'cgs': 'Oe',
        },
    },
    'magnetic moment': {
        'dimensions': {
            'mks': 'I * L^2',
            'cgs': '(M^1/2 * L^5/2) / T',
        },
        'units': {
            'mks': 'A * m^2',
            'cgs': 'Oe * cm^3',
        },
    },
    'magnetization': 'magnetic intensity',
    'magnetomotance': 'current',
    'mass': {
        'dimensions': {
            'mks': 'M',
            'cgs': 'M',
        },
        'units': {
            'mks': 'kg',
            'cgs': 'g',
        },
    },
    'mass density': 'mass / volume',
    'mass number': 'number',
    'momentum': {
        'dimensions': {
            'mks': '(M * L) / T',
            'cgs': '(M * L) / T',
        },
        'units': {
            'mks': 'kg * m / s',
            'cgs': 'g * cm / s',
        },
    },
    'momentum density': 'momentum / volume',
    'number': 'identity',
    'number density': '1 / volume',
    'particle distribution': '1 / (length * velocity)^3',
    'particle fluence': 'number / (area * solid_angle * energy / mass_number)',
    'particle flux': 'fluence / time',
    'permeability': {
        'dimensions': {
            'mks': '(M * L) / (I * T)^2',
            'cgs': '1',
        },
        'units': {
            'mks': 'H / m',
            'cgs': '1',
        },
    },
    'permittivity': {
        'dimensions': {
            'mks': 'T^4 * I / (M * L^3)',
            'cgs': '1',
        },
        'units': {
            'mks': 'F / m',
            'cgs': '1',
        },
    },
    'plane angle': { # See note about radian (Kalinin 2019).
        'dimensions': {
            'mks': '1',
            'cgs': '1',
        },
        'units': {
            'mks': 'rad',
            'cgs': 'rad',
        },
    },
    'polarization': 'charge / area',
    'potential': {
        'dimensions': {
            'mks': '(M * L^2) / (T^3 * I)',
            'cgs': '(M^1/2 * L^1/2) / T',
        },
        'units': {
            'mks': 'V',
            'cgs': 'statV',
        },
    },
    'power': {
        'dimensions': {
            'mks': 'M * L^2 / T^3',
            'cgs': 'M * L^2 / T^3',
        },
        'units': {
            'mks': 'W',
            'cgs': 'erg / s',
        },
    },
    'power density': 'power / volume',
    'pressure': {
        'dimensions': {
            'mks': 'M / (L * T^2)',
            'cgs': 'M / (L * T^2)',
        },
        'units': {
            'mks': 'Pa',
            'cgs': 'dyn / cm^2', # also barye (Ba)?
        },
    },
    'radioactivity': {
        'dimensions': {
            'mks': '1 / T',
            'cgs': '1 / T',
        },
        'units': {
            'mks': 'Bq',
            'cgs': 'Ci',
        },
    },
    'rate': 'frequency',
    'ratio': 'identity',
    'reluctance': {
        'dimensions': {
            'mks': '(I * T)^2 / (M * L^2)',
            'cgs': '1 / L',
        },
        'units': {
            'mks': 'A / Wb',
            'cgs': '1 / cm',
        },
    },
    'resistance': 'impedance',
    'resistivity': 'resistance * length',
    'temperature': {
        'dimensions': {
            'mks': 'H',
            'cgs': 'H',
        },
        'units': {
            'mks': 'K',
            'cgs': 'K',
        },
    },
    'thermal conductivity': 'power / (length * temperature)',
    'time': {
        'dimensions': {
            'mks': 'T',
            'cgs': 'T',
        },
        'units': {
            'mks': 's',
            'cgs': 's',
        },
    },
    'solid angle': { # See note about radian (Kalinin 2019).
        'dimensions': {
            'mks': '1',
            'cgs': '1',
        },
        'units': {
            'mks': 'sr',
            'cgs': 'sr',
        },
    },
    'speed': 'velocity',
    'vector potential': {
        'dimensions': {
            'mks': '(M * L) / (T^2 * I)',
            'cgs': '(M^1/2 * L^1/2) / T',
        },
        'units': {
            'mks': 'Wb / m',
            'cgs': 'G * cm',
        },
    },
    'velocity': 'length / time',
    'viscosity': {
        'dimensions': {
            'mks': 'M / (L * T)',
            'cgs': 'M / (L * T)',
        },
        'units': {
            'mks': 'kg / (m * s)',
            'cgs': 'P',
        },
    },
    'volume': 'length^3',
    'vorticity': 'frequency',
    'wavenumber': '1 / length',
    'work': 'energy',
}


C = 2.99792458e10
"""The speed of light in cm/s."""
PI = 3.141592653589793
"""The ratio of a circle's circumference to its diameter."""


_CONVERSIONS = {
    ('Rs', 'm'): 6.96e8,
    ('F', 'cm'): C**2 * 1e-9,
    ('C', 'statC'): 10*C,
    ('e', 'C'): 1.6022e-19,
    ('S', 'cm / s'): C**2 * 1e-5,
    ('A', 'statA'): 10*C,
    # ('C / m^2', 'statC / m^2'): 4*PI * C * 1e-3,
    ('Gy', 'erg / g'): 1e4,
    ('J', 'erg'): 1e7,
    ('eV', 'J'): 1.6022e-19,
    ('N', 'dyn'): 1e5,
    ('ohm', 's / cm'): 1e5 / C**2,
    ('H', 's^2 / cm'): 1e5 / C**2,
    # ('m', 'cm'): 1e2,
    ('au', 'm'): 1.495978707e11,
    ('Wb', 'Mx'): 1e8,
    ('T', 'G'): 1e4,
    ('A / m', 'Oe'): 4*PI * 1e-3,
    # ('A * m^2', 'Oe * cm^3'): 1e-3,
    # ('kg', 'g'): 1e3,
    ('nuc', 'kg'): 1.6605e-27,
    ('amu', 'kg'): 1.6605e-27,
    # ('kg * m / s', 'g * cm / s'): 1e5,
    ('H / m', '1'): 1e7 / 4*PI,
    ('F / m', '1'): 36*PI * 1e9,
    ('rad', 'deg'): 180 / PI,
    ('V', 'statV'): 1e6 / C,
    ('W', 'erg / s'): 1e7,
    ('Pa', 'dyn / cm^2'): 1e1,
    ('Bq', 'Ci'): 1.0 / 3.7e10,
    ('A / Wb', '1 / cm'): 4*PI * 1e-9,
    ('s', 'min'): 1.0 / 60.0,
    ('s', 'h'): 1.0 / 3600.0,
    ('s', 'd'): 1.0 / 86400.0,
    ('Wb / m', 'G * cm'): 1e6,
    ('kg / (m * s)', 'P'): 1e1,
}


class Graph(collections.abc.Collection):
    """A collection of weighted edges."""

    def __init__(self, base: typing.Mapping=None) -> None:
        """Initialize an instance.

        Parameters
        ----------
        base : mapping, optional
            The mapping of weighted connections from which to initialize this
            instance. Items in the mapping must have the form `((start, end),
            weight)`.
        """
        self.connections: typing.Dict[typing.Tuple[str, str], float] = {}
        """The forward and reverse links in this graph."""
        base = base or {}
        for (start, end), weight in base.items():
            self.add_connection(start, end, weight)

    def __contains__(self, connection: typing.Tuple[str, str]):
        """True if `connection` is available."""
        return connection in self.connections

    def __len__(self) -> int:
        """The number of connections. Called for len(self)."""
        return len(self.connections)

    def __iter__(self):
        """Iterate over connections. Called for iter(self)."""
        return iter(self.connections)

    @property
    def nodes(self):
        """The distinct nodes in this graph."""
        return {n for connection in self.connections for n in connection}

    def get_adjacencies(self, node: str):
        """Retrieve the connections to this node.
        
        Parameters
        ----------
        node : string
            The key corresponding to the target node.

        Returns
        -------
        `~dict`
            A dictionary whose keys represent the nodes connected to `node` and
            whose values represent the corresponding edge weight. An empty
            dictionary represents a node with no connections.
        """
        return {
            end: v for (start, end), v in self.connections.items()
            if start == node
        } if node in self.nodes else {}

    def get_weight(self, start: str, end: str):
        """Retrieve the weight of this link, if possible."""
        if (start, end) in self.connections:
            return self.connections[(start, end)]
        raise KeyError(
            f"No connection between {start!r} and {end!r}"
        ) from None

    def add_connection(
        self,
        start: str,
        end: str,
        weight: float=None,
    ) -> None:
        """Add a connection (with optional weight) to the graph."""
        forward = ((start, end), weight)
        inverse = 1 / weight if weight else weight
        reverse = ((end, start), inverse)
        for edge, value in (forward, reverse):
            if edge not in self.connections:
                self.connections[edge] = value

    def __repr__(self) -> str:
        """An unambiguous representation of this object."""
        return f"{self.__class__.__qualname__}({self})"

    def __str__(self) -> str:
        """A simplified representation of this object."""
        return '\n'.join(
            f"({n0} -> {n1}): {wt}"
            for (n0, n1), wt in self.connections.items()
        )


CONVERSIONS = Graph(_CONVERSIONS)

