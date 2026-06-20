# Cannon River — forward run

Runs MNiShed for the Cannon River catchment (Minnesota, USA) with a
fixed parameter set and prints goodness-of-fit diagnostics.

**Catchment:** Cannon River near Red Wing, MN (USGS 05355200; 3800 km²)
**Period:** 1992–1995 (daily)

## Requirements

- `mnished` installed (`pip install mnished`)

## Usage

```bash
python run_forward.py
```

Prints KGE\_logKGE\_logFDC, KGE\_logFDC, AIC, and BFI, then saves
`forward_run.png`.

## Files

| File | Description |
|------|-------------|
| `cannon_cfg.yml` | Model configuration (reservoirs, snowmelt, modules) |
| `CannonTestInput.csv` | Daily forcing and observed discharge |
| `run_forward.py` | Run script |

## Data attribution

Precipitation and temperature forcing are derived from the Livneh et al.
(2015) gridded meteorology product.  The catchment-average daily values were
extracted from the relevant grid cells using the
[LivnehPierce hydro-extractor](https://github.com/MNiMORPH/LivnehPierce-hydro-extractor)
(Pierce et al., 2021).

Observed streamflow is from the USGS National Water Information System
(gauge 05355200).

**References**

- Livneh, B., et al. (2015). A spatially comprehensive, hydrometeorological
  data set for Mexico, the U.S., and Southern Canada 1950–2013. *Scientific
  Data*, 2, 150042. https://doi.org/10.1038/sdata.2015.42
- Pierce, A., et al. (2021). LivnehPierce hydro-extractor.
  https://github.com/MNiMORPH/LivnehPierce-hydro-extractor

## Parameter estimation

To calibrate the parameters in `cannon_cfg.yml` see the companion example
in `../cannon_inverse/`.
