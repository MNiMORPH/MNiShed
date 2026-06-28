"""
Generate the autumn-senescence figure for the documentation.

Contrasts the two ``senescence_method`` forms of the phenology coefficient
(:meth:`mnished.Buckets.phenology_Kc`): a fixed calendar window (``doy``) versus
day-length-driven brown-down (``photoperiod``). The point of the figure is
latitude transfer.

* *Left:* the autumn brown-down at a single latitude (~47.5°N). The two forms
  give similar autumn shapes; the photoperiod form is anchored to a critical day
  length *below the equinox*, which is what makes it transfer correctly (right).
* *Right:* the brown-down midpoint date versus latitude. A fixed calendar window
  is, by construction, the same date everywhere (flat line). The photoperiod cue
  moves with latitude in the physically correct direction — higher latitudes
  senesce earlier — because a sub-equinox day length is reached sooner where the
  day shortens faster. The slope is modest (~0.6 day per degree): pure photoperiod
  captures part of the autumn gradient, not all of it (temperature also matters;
  see the NDVI complement, issue #26), but unlike a fixed date it needs no
  per-basin re-dating.

Day length is computed from a standard astronomical declination model, and the
senescence ramps use the exact expressions in ``phenology_Kc``. No model run is
needed — the figure illustrates the Kc mechanism itself.

Usage (from the repository root)::

    python docs/figures/plot_phenology_senescence.py

Writes ``docs/source/_static/phenology_senescence.png``.
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after Agg backend)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
OUT = os.path.join(ROOT, "docs", "source", "_static", "phenology_senescence.png")

# Phenology defaults (match mnished phenology_params)
CRIT_HR, SPAN_HR = 11.0, 1.8           # senescence_photoperiod__hr / _span__hr
START_DOY, END_DOY = 260, 305          # senescence_start_doy / _end_doy
SOLSTICE_DOY = 172                     # post-solstice gate
CROW_WING_LAT = 47.5


def day_length(latitude_deg, doy):
    """Astronomical day length [hr] at a latitude for a day-of-year array."""
    decl = np.radians(23.44) * np.sin(2 * np.pi * (doy - 81) / 365.0)
    lat = np.radians(latitude_deg)
    h = np.clip(-np.tan(lat) * np.tan(decl), -1.0, 1.0)
    return (24.0 / np.pi) * np.arccos(h)


def senesce_photoperiod(latitude_deg, doy):
    """Day-length-driven brown-down fraction (0 = full canopy, 1 = dormant)."""
    N = day_length(latitude_deg, doy)
    s = np.clip((CRIT_HR - N) / SPAN_HR, 0.0, 1.0)
    return np.where(doy >= SOLSTICE_DOY, s, 0.0)


def senesce_doy(doy):
    """Fixed-calendar brown-down fraction (latitude-independent)."""
    return np.clip((doy - START_DOY) / (END_DOY - START_DOY), 0.0, 1.0)


def midpoint_doy(frac, doy):
    """Day-of-year where a monotonic brown-down fraction crosses 0.5."""
    return float(np.interp(0.5, frac, doy))


def main():
    doy = np.arange(182, 366)          # Jul -> Dec, the brown-down half-year
    months = [182, 213, 244, 274, 305, 335, 365]
    mlabels = ["Jul", "Aug", "Sep", "Oct", "Nov", "Dec", ""]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.6, 3.4),
                                   gridspec_kw=dict(wspace=0.28))

    # --- Left: the two brown-down forms at the Crow Wing latitude ----------
    axL.plot(doy, 1 - senesce_doy(doy), color="#a63", lw=2.2,
             label="calendar (doy)")
    axL.plot(doy, 1 - senesce_photoperiod(CROW_WING_LAT, doy), color="#2a7",
             lw=2.2, ls="--", label="photoperiod")
    axL.set_title(f"Autumn brown-down at {CROW_WING_LAT:.0f}°N", fontsize=9.5)
    axL.set_ylabel("canopy fraction remaining\n(1 − senescence)")
    axL.set_xticks(months)
    axL.set_xticklabels(mlabels)
    axL.set_xlim(182, 365)
    axL.set_ylim(-0.03, 1.05)
    axL.set_xlabel("month")
    axL.legend(fontsize=8, loc="lower left", framealpha=0.9)

    # --- Right: brown-down midpoint date vs latitude -----------------------
    lats = np.linspace(35, 55, 41)
    mid_pp = np.array([midpoint_doy(senesce_photoperiod(la, doy), doy)
                       for la in lats])
    mid_doy = midpoint_doy(senesce_doy(doy), doy)
    axR.plot(mid_pp, lats, color="#2a7", lw=2.4, label="photoperiod")
    axR.axvline(mid_doy, color="#a63", lw=2.0, ls=":",
                label="calendar (fixed everywhere)")
    axR.axhline(CROW_WING_LAT, color="gray", lw=0.8, ls="--")
    axR.text(mid_pp[0] - 1, CROW_WING_LAT + 0.5, "Crow Wing", fontsize=7.5,
             color="gray", ha="right")
    axR.set_title("Brown-down date tracks latitude\n(photoperiod only)",
                  fontsize=9.5)
    axR.set_xlabel("senescence midpoint (day of year)")
    axR.set_ylabel("latitude (°N)")
    axR.legend(fontsize=8, loc="lower right", framealpha=0.9)

    fig.savefig(OUT, dpi=130, bbox_inches="tight")

    slope = (mid_pp[-1] - mid_pp[0]) / (lats[-1] - lats[0])
    print(f"wrote {os.path.relpath(OUT, ROOT)}")
    print(f"  photoperiod midpoint: {mid_pp[0]:.0f} d @35°N -> "
          f"{mid_pp[-1]:.0f} d @55°N  (slope {slope:.2f} day/°, north earlier)")
    print(f"  calendar midpoint: {mid_doy:.0f} d (all latitudes)")


if __name__ == "__main__":
    main()
