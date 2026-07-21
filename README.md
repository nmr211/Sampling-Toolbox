# Before Starting
An introduction to the Sampling-Toolbox is provided in the paper.

Paper Name + DOI

We recomend reading this (short) paper before getting started.

# Sampling-Toolbox

Toolbox used to generate Latin Hypercube Samples (LHS) for Space Weather Indices.

The workflow:
- Loads observational datasets for solar wind, geomagnetic, magnetospheric, atmospheric winds
- Fits statistical distributions (single and bimodal models) to individual datasets
- Generates LHS based on selected input parameters
- Output file generates samples of chosen input parameters in NetCDF format (`lhs_samples.nc`)

---

## 📦 Features

- Latin Hypercube Sampling for efficient parameter space coverage
- Single and bimodal distribution fitting to 50 years of observational data for each input parameter
- Supports multiple space weather variables:
  - IMF components: `bx`, `by`, `bz`
  - Solar wind: `swvel`, `swden`
  - Geomagnetic indices: `al`
  - Atmospheric winds: `un`
  - Solar flux: `f107`
- NetCDF output for easy scientific integration
- Simple runner script for non-technical users

---

## 🚀 Getting Started

### 1. Download and Unzip the Data

Download the toolbox files and unzip the provided datasets called Indices.zip, AL_data.zip, and Zonal_Wind_Data.zip.

After extracting the data, ensure the following data folders/files are available in the toolbox directory:

- `AL_data` — geomagnetic index data
- `Indices` — space weather index datasets
- `Zonal wind` — atmospheric zonal wind observations

These datasets contain the observational records used to fit the statistical distributions for each input parameter.

---

### 2. Configure the Sampling Options

The sampling process is controlled through the script. 
Edit `run_toolbox.py` to define the number of samples, variables, and output filename.

```python
SAMPLES = 10000

VARIABLES = [
    "bx", "by", "bz",
    "swvel", "swden",
    "al", "un", "f107"
]

OUTPUT_FILE = "lhs_samples.nc"
```

### 3. Run Toolbox

Once the datasets have been downloaded, extracted, and the configuration in `run_toolbox.py` has been updated, the toolbox can be executed.

Run the toolbox from the command line using:

```bash
python run_toolbox.py
```

The script will then generate the requested number of samples.

After completion, the output file will be created:

```text
lhs_samples.nc
```

This file contains the generated samples for the selected variables and can be used as input for further scientific analysis or modelling applications.
