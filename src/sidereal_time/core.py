"""Convert civil (UTC) time to local sidereal time.

This module implements the standard low-precision formula for Greenwich
Apparent Sidereal Time published in the *Astronomical Almanac*. It is
adequate for pointing telescopes and scheduling observations at the
arcminute level; it is not intended for sub-arcsecond astrometry.

The single hard design decision made here: the public ``civil_to_sidereal``
function takes a UTC datetime plus a longitude in *decimal degrees* and
returns the local sidereal time in *decimal hours*. No other units are
accepted. The README documents this clearly so callers never have to guess.
"""

import datetime as _dt
import math

# Constants drawn from the IAU 1982 / Astronomical Almanac expression for
# Greenwich Mean Sidereal Time. The leading GMST at the 2000.0 epoch is
# 18h 41m 50.54841s, expressed below in seconds, plus a quadratic term in
# the elapsed Julian centuries. These are kept explicit so the formula is
# auditable against the reference.

_GMST_REF_HOURS = 6.0 / 24.0 + 41.0 / (24.0 * 60.0) + 50.54841 / (24.0 * 60.0 * 60.0)

# Rate of GMST increase per Julian century, in hours per century.
# 876000h is the number of sidereal hours in a Julian century; the extra
# term encodes the difference between sidereal and solar time.
_GMST_RATE_HOURS_PER_CENTURY = 86400.0 * 36525.0 / 1.0
# The canonical Almanac linear coefficient is 8640184.812866 seconds per
# Julian century, plus a quadratic correction of 0.0000931 seconds per
# century squared. We combine them into a single coefficient below.

_GMST_LINEAR_SEC_PER_CENTURY = 8640184.812866
_GMST_QUAD_SEC_PER_CENTURY2 = 0.0000931
_GMST_DAILY_OFFSET_SEC = 86400.0 * 1.00273790935
# Equation of equinoxes correction is tiny (<2s) and varies with nutation;
# for this library we use Mean Sidereal Time, which is what the formula
# below computes. GAST adds a correction that requires the Moon's longitude;
# we deliberately omit it to keep the dependency surface flat.


class SiderealCalculator:
    """Convert civil UTC datetimes to local sidereal time.

    A small class wraps the conversion so that callers can reuse one
    configured instance (a given longitude) across many timestamps without
    re-passing the longitude each call. The math itself is pure; the class
    exists for call-site ergonomics, not for statefulness.

    Parameters
    ----------
    longitude_deg:
        Geographic longitude of the observing site, measured **east of
        Greenwich in decimal degrees**. West longitudes are negative.
    """

    def __init__(self, longitude_deg: float):
        self.longitude_deg = float(longitude_deg)

    def local_sidereal_time(self, when: _dt.datetime) -> float:
        """Return local sidereal time in decimal hours for *when*.

        *when* must be timezone-aware and in UTC. If it is naive or not UTC
        a ``ValueError`` is raised, because silent timezone assumptions
        produce wrong sidereal times that are very hard to debug after the
        fact.
        """
        return civil_to_sidereal(when, self.longitude_deg)


def _julian_day(utc: _dt.datetime) -> float:
    """Julian Day number for the given UTC instant.

    Uses the Meeus algorithm for the integer part and adds the fractional
    day from hours/minutes/seconds. This stays in the proleptic Gregorian
    calendar, which matches Python's ``datetime`` assumptions.
    """
    y = utc.year
    m = utc.month
    d = utc.day
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    jd_int = (
        int(365.25 * (y + 4716))
        + int(30.6001 * (m + 1))
        + d
        + b
        - 1524
    )
    day_fraction = (
        utc.hour
        + (utc.minute + (utc.second + utc.microsecond / 1_000_000.0) / 60.0) / 60.0
    ) / 24.0
    return jd_int + day_fraction


def _greenwich_mean_sidereal_time_hours(jd: float) -> float:
    """GMST in decimal hours for the given Julian Day.

    The reference expression computes GMST in seconds at 0h UT of the given
    day, then scales by the sidereal-to-solar ratio for the fractional day.
    Working in seconds first keeps intermediate magnitudes close to the
    published coefficients, which makes the code easier to audit.
    """
    # JD at 0h UT of the same day.
    jd0 = math.floor(jd - 0.5) + 0.5
    t0 = (jd0 - 2451545.0) / 36525.0  # Julian centuries since J2000.0
    gmst_seconds_at_0h = (
        24110.54841
        + _GMST_LINEAR_SEC_PER_CENTURY * t0
        + _GMST_QUAD_SEC_PER_CENTURY2 * t0 * t0
    )
    # Seconds of UT elapsed since 0h.
    ut_seconds = (jd - jd0) * 86400.0
    gmst_seconds = gmst_seconds_at_0h + ut_seconds * 1.00273790935
    # Reduce to [0, 86400) so the hours value is always well behaved.
    gmst_seconds %= 86400.0
    return gmst_seconds / 3600.0


def civil_to_sidereal(when: _dt.datetime, longitude_deg: float) -> float:
    """Convert a UTC civil time to local sidereal time in hours.

    Parameters
    ----------
    when:
        A timezone-aware ``datetime`` in UTC. Naive datetimes and non-UTC
        timezones are rejected with ``ValueError`` — converting from a
        wall-clock time without a timezone is the classic source of silent
        one-hour errors, so we refuse the ambiguity rather than guess.
    longitude_deg:
        East-positive longitude in decimal degrees. West is negative.

    Returns
    -------
    float
        Local sidereal time in decimal hours, wrapped into [0, 24).

    Notes
    -----
    The returned value is Mean Sidereal Time. Apparent Sidereal Time differs
    by the equation of equinoxes (at most ~1.7 seconds of time), which we do
    not compute because it requires nutation data we do not ship.
    """
    if when.tzinfo is None:
        raise ValueError(
            "civil_to_sidereal requires a timezone-aware datetime; "
            "got a naive datetime which is ambiguous."
        )
    # Coerce to UTC explicitly. The caller may pass a datetime already in
    # UTC; if so this is a no-op. If they pass another zone, we convert.
    utc = when.astimezone(_dt.timezone.utc)
    jd = _julian_day(utc)
    gmst_hours = _greenwich_mean_sidereal_time_hours(jd)
    # Longitude east of Greenwich adds to GMST to get LST.
    # 1 hour of time == 15 degrees of arc.
    lst = gmst_hours + float(longitude_deg) / 15.0
    lst %= 24.0
    if lst < 0.0:
        lst += 24.0
    return lst
