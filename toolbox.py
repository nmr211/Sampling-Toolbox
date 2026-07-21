# -*- coding: utf-8 -*-

import argparse
import netCDF4 as nc
import numpy as np
import h5py
from scipy import stats
from scipy.optimize import minimize
from scipy.stats import qmc
from netCDF4 import Dataset


# =============================================================================
# Functions (Bimodal likelihoods)
# =============================================================================

def bimodal_johnsonsu_nll(params, data):
    w = params[0]
    a1, b1, loc1, scale1 = params[1:5]
    a2, b2, loc2, scale2 = params[5:9]

    if w <= 0 or w >= 1:
        return np.inf
    if scale1 <= 0 or scale2 <= 0:
        return np.inf
    if b1 <= 0 or b2 <= 0:
        return np.inf

    pdf = (
        w * stats.johnsonsu.pdf(data, a1, b1, loc1, scale1)
        + (1 - w) * stats.johnsonsu.pdf(data, a2, b2, loc2, scale2)
    )

    return -np.sum(np.log(pdf + 1e-12))


def bimodal_gennorm_nll(params, data):
    w = params[0]
    b1, loc1, scale1 = params[1:4]
    b2, loc2, scale2 = params[4:7]

    if w <= 0 or w >= 1:
        return np.inf
    if b1 <= 0 or b2 <= 0:
        return np.inf
    if scale1 <= 0 or scale2 <= 0:
        return np.inf

    pdf = (
        w * stats.gennorm.pdf(data, b1, loc1, scale1)
        + (1 - w) * stats.gennorm.pdf(data, b2, loc2, scale2)
    )

    return -np.sum(np.log(pdf + 1e-12))


def bimodal_skewnorm_nll(params, data):
    w = params[0]
    a1, loc1, scale1 = params[1:4]
    a2, loc2, scale2 = params[4:7]

    if w <= 0 or w >= 1:
        return np.inf
    if scale1 <= 0 or scale2 <= 0:
        return np.inf

    pdf = (
        w * stats.skewnorm.pdf(data, a1, loc1, scale1)
        + (1 - w) * stats.skewnorm.pdf(data, a2, loc2, scale2)
    )

    return -np.sum(np.log(pdf + 1e-12))


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Latin Hypercube Sampling for space weather variables"
    )

    parser.add_argument("--samples", type=int, default=10000,
                        help="Number of samples to generate")

    parser.add_argument("--variables", nargs="+",
                        default=["bx", "by", "bz", "swvel", "swden", "al", "un", "f107"],
                        help="Variables to sample")

    parser.add_argument("--output", type=str, default="lhs_samples.nc",
                        help="Output NetCDF file")

    return parser.parse_args()


# =============================================================================
# Data loading
# =============================================================================

def load_data():

    with h5py.File('Indices/geo500101g.002.hdf5', 'r') as GPI_dataset:
        data_group = GPI_dataset['Data']['Table Layout']
        f107_data = data_group['f10.7'][:]
        f107 = f107_data[np.isfinite(f107_data)] * 1e22  # Convert to SFU

    with h5py.File('Indices/imf631127g.002.hdf5', 'r') as IMF_dataset:
        IMF_table = IMF_dataset['Data']['Table Layout']

        bx = IMF_table['bxgsm'][:][np.isfinite(IMF_table['bxgsm'][:])] * 1e9
        by = IMF_table['bygsm'][:][np.isfinite(IMF_table['bygsm'][:])] * 1e9
        bz = IMF_table['bzgsm'][:][np.isfinite(IMF_table['bzgsm'][:])] * 1e9
        swden = IMF_table['swden'][:][np.isfinite(IMF_table['swden'][:])]
        swvel = IMF_table['swspd'][:][np.isfinite(IMF_table['swspd'][:])]
        

    with h5py.File('AL_data.nc', 'r') as f:
        al = f['AL'][:]

    with nc.Dataset('Zonal_Wind_Data.nc', 'r') as f:
        un = f['un'][:]

    data = {
        "bx": bx[np.isfinite(bx)],
        "by": by[np.isfinite(by)],
        "bz": bz[np.isfinite(bz)],
        "swden": swden[np.isfinite(swden)],
        "swvel": swvel[np.isfinite(swvel)],
        "al": al[np.isfinite(al)],
        "un": un[np.isfinite(un)],
        "f107": f107[np.isfinite(f107)],
    }

    return data


# =============================================================================
# Distribution fitting
# =============================================================================

def fit_distributions(data):
    fits = {}

    fits["bx"] = stats.t.fit(data["bx"])
    fits["by"] = stats.t.fit(data["by"])
    fits["bz"] = stats.johnsonsu.fit(data["bz"])

    fits["swvel"] = stats.skewnorm.fit(data["swvel"])
    fits["swden"] = stats.johnsonsu.fit(data["swden"])

    return fits


def fit_bimodal(data):
    results = {}

    p0_al = [0.5, 1, 1, np.percentile(data["al"], 25), np.std(data["al"]),
             1, 1, np.percentile(data["al"], 75), np.std(data["al"])]

    results["al"] = minimize(bimodal_johnsonsu_nll, p0_al,
                             args=(data["al"],), method="L-BFGS-B")

    p0_un = [0.5, 2, np.percentile(data["un"], 25), np.std(data["un"]),
             2, np.percentile(data["un"], 75), np.std(data["un"])]

    results["un"] = minimize(bimodal_gennorm_nll, p0_un,
                             args=(data["un"],), method="L-BFGS-B")

    p0_f107 = [0.5, 0, np.percentile(data["f107"], 25), np.std(data["f107"]),
               0, np.percentile(data["f107"], 75), np.std(data["f107"])]

    results["f107"] = minimize(bimodal_skewnorm_nll, p0_f107,
                               args=(data["f107"],), method="L-BFGS-B")

    return results


# =============================================================================
# Sampling helpers
# =============================================================================

def clip(x, data):
    return np.clip(x, np.min(data), np.max(data))


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    data = load_data()

    fits = fit_distributions(data)
    bimodal = fit_bimodal(data)

    n = args.samples
    selected = set(args.variables)

    sampler = qmc.LatinHypercube(d=1)
    u = sampler.random(n=n).flatten()

    samples = {}

    # -------------------------
    # Continuous distributions
    # -------------------------

    if "bx" in selected:
        samples["bx"] = clip(stats.t.ppf(u, *fits["bx"]), data["bx"])

    if "by" in selected:
        samples["by"] = clip(stats.t.ppf(u, *fits["by"]), data["by"])

    if "bz" in selected:
        samples["bz"] = clip(stats.johnsonsu.ppf(u, *fits["bz"]), data["bz"])

    if "swvel" in selected:
        samples["swvel"] = clip(stats.skewnorm.ppf(u, *fits["swvel"]), data["swvel"])

    if "swden" in selected:
        samples["swden"] = clip(stats.johnsonsu.ppf(u, *fits["swden"]), data["swden"])

    # -------------------------
    # Bimodal distributions
    # -------------------------

    if "al" in selected:
        p = bimodal["al"].x
        w = p[0]
        a1, b1, loc1, scale1 = p[1:5]
        a2, b2, loc2, scale2 = p[5:9]

        u2 = np.random.rand(n)
        raw = np.where(
            u2 < w,
            stats.johnsonsu.rvs(a1, b1, loc1, scale1, size=n),
            stats.johnsonsu.rvs(a2, b2, loc2, scale2, size=n),
        )
        samples["al"] = clip(raw, data["al"])

    if "un" in selected:
        p = bimodal["un"].x
        w = p[0]
        b1, loc1, scale1 = p[1:4]
        b2, loc2, scale2 = p[4:7]

        u2 = np.random.rand(n)
        raw = np.where(
            u2 < w,
            stats.gennorm.rvs(b1, loc1, scale1, size=n),
            stats.gennorm.rvs(b2, loc2, scale2, size=n),
        )
        samples["un"] = clip(raw, data["un"])

    if "f107" in selected:
        p = bimodal["f107"].x
        w = p[0]
        a1, loc1, scale1 = p[1:4]
        a2, loc2, scale2 = p[4:7]

        u2 = np.random.rand(n)
        raw = np.where(
            u2 < w,
            stats.skewnorm.rvs(a1, loc1, scale1, size=n),
            stats.skewnorm.rvs(a2, loc2, scale2, size=n),
        )
        samples["f107"] = clip(raw, data["f107"])

    # -------------------------
    # Save output
    # -------------------------

    with Dataset(args.output, "w", format="NETCDF4") as ncfile:
        ncfile.createDimension("sample", n)

        for name, arr in samples.items():
            var = ncfile.createVariable(name, "f8", ("sample",))
            var[:] = arr
            var.description = f"LHS sampled variable: {name}"

    print("Sampling complete.")
    print({k: v.shape for k, v in samples.items()})
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()