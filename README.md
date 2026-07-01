# Sampling-Toolbox
Toolbox used to generate Latin Hypercube Samples for Space Weather Indicies

The workflow:
- Loads observational datasets (solar wind, geomagnetic indices, atmospheric winds)
- Fits statistical distributions (single and bimodal models)
- Generates samples using LHS
- Outputs results in NetCDF format

---

## 📦 Features

- Latin Hypercube Sampling for efficient parameter space coverage
- Single and bimodal distribution fitting to 50 years of observational data for each input parameter
- Supports multiple space weather variables:
  - IMF components: `bx`, `by`, `bz`
  - Solar wind: `swvel`, `swden`
  - Geomagnetic indices: `al`
  - Atmospheric zonal wind: `un`
  - Solar flux: `f107`
- NetCDF output for easy scientific integration
- Simple runner script for non-technical users

---

