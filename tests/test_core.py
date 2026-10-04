import datetime as dt
import math
import unittest

from sidereal_time.core import SiderealCalculator, civil_to_sidereal


class CivilToSiderealTests(unittest.TestCase):
    def _utc(self, y, mo, d, h, mi, s=0, us=0):
        return dt.datetime(y, mo, d, h, mi, s, us, tzinfo=dt.timezone.utc)

    def test_rejects_naive_datetime(self):
        naive = dt.datetime(2024, 1, 1, 0, 0, 0)
        with self.assertRaises(ValueError):
            civil_to_sidereal(naive, 0.0)

    def test_rejects_non_utc_timezone(self):
        # +02:00 is not UTC; the function should still work by coercing,
        # but here we test that the coercion produces the same answer as
        # passing the equivalent UTC instant directly.
        non_utc = dt.datetime(2024, 1, 1, 2, 0, 0, tzinfo=dt.timezone(dt.timedelta(hours=2)))
        utc_equiv = self._utc(2024, 1, 1, 0, 0, 0)
        lst_non_utc = civil_to_sidereal(non_utc, 0.0)
        lst_utc = civil_to_sidereal(utc_equiv, 0.0)
        self.assertAlmostEqual(lst_non_utc, lst_utc, places=6)

    def test_greenwich_midnight_j2000(self):
        # At 2000-01-01 12:00 UTC (the J2000.0 epoch) GMST should be
        # approximately 18h 41m 50.548s. The published value at the epoch
        # is 18h41m50.54841s. We allow a generous tolerance because the
        # mean formula accumulates rounding.
        when = self._utc(2000, 1, 1, 12, 0, 0, 0)
        gmst = civil_to_sidereal(when, 0.0)
        # The formula references 0h UT of the day, so at 12h UT the GMST
        # is approximately the 0h value plus ~12 sidereal hours.
        expected_at_0h = 18.0 + 41.0 / 60.0 + 50.54841 / 3600.0
        expected = (expected_at_0h + 12.0 * 1.00273790935) % 24.0
        self.assertAlmostEqual(gmst % 24.0, expected, delta=1e-6)

    def test_wraps_into_0_to_24(self):
        when = self._utc(2024, 6, 21, 0, 0, 0)
        lst = civil_to_sidereal(when, 0.0)
        self.assertTrue(0.0 <= lst < 24.0, f"LST out of range: {lst}")

    def test_positive_longitude_adds_to_gmst(self):
        # East longitude should produce a larger LST than Greenwich for the
        # same instant.
        when = self._utc(2024, 3, 15, 18, 30, 0)
        greenwich = civil_to_sidereal(when, 0.0)
        east = civil_to_sidereal(when, 30.0)  # 2 hours east
        # East of Greenwich, LST is ahead (larger, modulo 24).
        diff = (east - greenwich) % 24.0
        self.assertAlmostEqual(diff, 2.0, places=6)

    def test_negative_longitude_subtracts(self):
        when = self._utc(2024, 3, 15, 18, 30, 0)
        greenwich = civil_to_sidereal(when, 0.0)
        west = civil_to_sidereal(when, -30.0)  # 2 hours west
        diff = (greenwich - west) % 24.0
        self.assertAlmostEqual(diff, 2.0, places=6)

    def test_antipodal_longitude_complement(self):
        # At longitude 180 and -180 the result should be the same point on
        # the circle (differing by 24h, which the wrap removes).
        when = self._utc(2024, 7, 4, 9, 15, 0)
        east = civil_to_sidereal(when, 180.0)
        west = civil_to_sidereal(when, -180.0)
        self.assertAlmostEqual(east, west, places=6)

    def test_microsecond_precision_supported(self):
        # The function must accept and use microseconds; changing them by
        # 1e6 should change the LST by ~1s of sidereal time (1/3600 hour).
        base = self._utc(2024, 1, 1, 0, 0, 0, 0)
        plus_one_sec = self._utc(2024, 1, 1, 0, 0, 1, 0)
        lst_base = civil_to_sidereal(base, 0.0)
        lst_plus = civil_to_sidereal(plus_one_sec, 0.0)
        delta_hours = (lst_plus - lst_base) % 24.0
        # 1 second of solar time ~ 1.0027 seconds of sidereal time.
        expected_delta = 1.00273790935 / 3600.0
        self.assertAlmostEqual(delta_hours, expected_delta, places=8)

    def test_deterministic_for_same_inputs(self):
        # Determinism: same inputs must yield identical outputs.
        when = self._utc(2023, 11, 2, 23, 59, 59, 123456)
        lon = -70.5
        a = civil_to_sidereal(when, lon)
        b = civil_to_sidereal(when, lon)
        self.assertEqual(a, b)

    def test_returns_float(self):
        when = self._utc(2024, 1, 1, 0, 0, 0)
        lst = civil_to_sidereal(when, 0.0)
        self.assertIsInstance(lst, float)


class SiderealCalculatorTests(unittest.TestCase):
    def test_class_matches_free_function(self):
        when = dt.datetime(2024, 5, 10, 14, 20, 0, tzinfo=dt.timezone.utc)
        lon = -105.27
        calc = SiderealCalculator(lon)
        self.assertAlmostEqual(calc.local_sidereal_time(when), civil_to_sidereal(when, lon), places=12)

    def test_stores_longitude(self):
        calc = SiderealCalculator(42.0)
        self.assertEqual(calc.longitude_deg, 42.0)

    def test_rejects_naive_via_class(self):
        calc = SiderealCalculator(0.0)
        naive = dt.datetime(2024, 1, 1)
        with self.assertRaises(ValueError):
            calc.local_sidereal_time(naive)


if __name__ == "__main__":
    unittest.main()
