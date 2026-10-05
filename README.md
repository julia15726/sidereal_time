# Sidereal Time

Convert civil (UTC) time to local sidereal time for a given east-positive longitude. Standard library only, no dependencies.

## Usage

```python
import datetime
from sidereal_time import SiderealCalculator, civil_to_sidereal

when = datetime.datetime(2024, 6, 21, 0, 0, 0, tzinfo=datetime.timezone.utc)

# Functional form:
lst = civil_to_sidereal(when, longitude_deg=-70.5)  # decimal hours

# Or a reusable calculator pinned to one site:
site = SiderealCalculator(longitude_deg=-70.5)
lst = site.local_sidereal_time(when)
```

`civil_to_sidereal(when, longitude_deg)` returns a `float` in decimal hours on `[0, 24)`. `SiderealCalculator(longitude_deg).local_sidereal_time(when)` does the same with the longitude baked in.

## Why this exists

Pointing a telescope or planning an observation needs the local sidereal time: the right ascension currently on the meridian. Computing it well to sub-arcsecond accuracy requires nutation and Earth-orientation-parameter corrections that pull in nontrivial data. For most amateur and many professional scheduling needs, the IAU 1982 mean-sidereal-time formula is plenty accurate (arcminute level) and depends on nothing but the date.

This library implements that formula and nothing else. The trade-off is simplicity over the last word in precision: we compute **Mean** Sidereal Time, not Apparent. The difference is the equation of equinoxes, at most ~1.7 seconds of time (~0.07 degrees of arc). If you need that, use a heavier library.

## Units and the edge you will hit

Longitude is **decimal degrees, east-positive** (west is negative). The return is **decimal hours**. Mixing these up is the single most common mistake: 15 degrees equals one hour, so a sign error or a units swap produces a silently wrong answer that still looks plausible. The functions accept only these units; there is no flag to switch to radians or H/M/S.

The datetime passed in must be **timezone-aware**. Naive datetimes are rejected with `ValueError`, because a wall-clock time with no zone is ambiguous and a wrong guess is invisible. The library coerces non-UTC zones to UTC internally, so passing `+02:00` is fine and equivalent to passing the corresponding UTC instant.
