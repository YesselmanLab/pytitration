import re
import os
import json
import logging
import sys
from typing import List, Dict
import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
import pickle

from rna_map.mutation_histogram import (
    MutationHistogram,
    convert_dreem_mut_histos_to_mutation_histogram,
)
from dreem.bit_vector import MutationHistogram as DREEMMutationHistogram

# logging #####################################################################

APP_LOGGER_NAME = "PYTITRATION"


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=False, file_name=None):
    """
    Set up the logger for the app
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    # pylint: disable=C0103
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(module_name):
    """
    Get the logger for the module
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)


log = get_logger("CLI")

# titration functions ##########################################################


def normalize_data_full(data):
    if np.min(data) == np.max(data):
        return data
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def normalize_data(data):
    if np.min(data) == np.max(data):
        return data
    return (data) / (np.max(data))


def normalized_hill_equation(conc, K, n):
    """
    Assumes a range from 0 to 1
    :param conc: concentration of titration agent either mg2+ or a ligand
    :param K: dissociation constant
    :param n: hill coefficient
    """
    return ((conc / K) ** n) / (1 + (conc / K) ** n)


def fit_bootstrap(p0, x, y, function, n_runs=100, n_sigma=1.0):
    """
    Uses bootstrap method to estimate the 1 sigma confidence interval of the parameters
    for a fit to data (x,y) with function function(x, params).

    1 sigma corresponds to 68.3% confidence interval
    2 sigma corresponds to 95.44% confidence interval

    :param p0: initial guess for parameters
    :param x: x - independent values
    :param y: y - dependent values, what are you trying to fit to
    :param function: function to fit to, should be a python function
    :param n_runs: number of bootstrap runs - default 100
    :param n_sigma: number of sigma to use for confidence interval - default 1.0
    """
    errfunc = lambda p, x, y: function(x, p[0], p[1]) - y
    # Fit first time
    pfit, perr = optimize.leastsq(errfunc, p0, args=(x, y), full_output=0)
    # Get the stdev of the residuals
    residuals = errfunc(pfit, x, y)
    sigma_res = np.std(residuals)
    sigma_err_total = np.sqrt(sigma_res**2)
    # 100 random data sets are generated and fitted
    ps = []
    for i in range(n_runs):
        random_delta = np.random.normal(0.0, sigma_err_total, len(y))
        random_y = y + random_delta
        random_fit, _ = optimize.leastsq(errfunc, p0, args=(x, random_y), full_output=0)
        ps.append(random_fit)
    ps = np.array(ps)
    mean_pfit = np.mean(ps, 0)
    err_pfit = n_sigma * np.std(ps, 0)
    return mean_pfit, err_pfit


# plotting functions ###########################################################


def generate_titration_plot(df):
    pstart = [1, 1]
    norm_data = -normalize_data(np.array(df["avg"])) + 1
    pfit, perr = fit_bootstrap(pstart, df["conc"], norm_data, normalized_hill_equation)
    # print(g[["mg_conc", "gaaa_avg"]])
    plt.scatter(df["conc"], norm_data)
    xs, ys = [], []
    for j in np.arange(0, df["conc"].max(), 0.25):
        y = normalized_hill_equation(j, pfit[0], pfit[1])
        xs.append(j)
        ys.append(y)
    plt.plot(xs, ys)
    plt.ylim(-0.05, 1.1)
    plt.ylabel("Normalized Reactivity")
    plt.xlabel("Concentration (mM)")


# helper functions #############################################################


def str_to_range(x):
    """
    Convert a string representation of a range of numbers to a list of integers.

    Given a string representation of a range of numbers, this function returns a
    list of integers corresponding to the numbers in the range. The string can
    contain single numbers separated by commas, and ranges of numbers separated by
    a hyphen.

    :param x: A string representation of a range of numbers.
    :return: A list of integers corresponding to the numbers in the range.
    """
    return sum(
        (
            i if len(i) == 1 else list(range(i[0], i[1] + 1))
            for i in (
                [int(j) for j in i if j] for i in re.findall(r"(\d+),?(?:-(\d+))?", x)
            )
        ),
        [],
    )


def get_pop_avg(mut_histo: Dict) -> List[float]:
    """
    Returns the population average of the histogram
    :param inc_del: if True, include deletions in the average
    """
    nuc_coords = list(range(mut_histo["start"], mut_histo["end"] + 1))
    pop_avg = []
    for pos in nuc_coords:
        try:
            mut_frac = mut_histo["mut_bases"][pos] / mut_histo["info_bases"][pos]
        except:
            mut_frac = 0.0
        pop_avg.append(round(mut_frac, 5))
    return pop_avg


def sub_select_data(data: List, selection: List) -> List:
    """
    Select a subset of data based on a list of indices.

    Given a list of data, this function returns a subset of the data based on a
    list of indices.

    :param data: A list of data.
    :param selection: A list of indices.
    :return: A subset of the data.
    """
    return [data[i] for i in selection]


@click.command()
@click.option(
    "-c", "--csv", help="Input csv file with dir and conc columns", required=True
)
@click.option("-n", "--name", help="the construct name to use", required=True)
@click.option(
    "-r",
    "--range",
    "nuc_range",
    help="the nucleotides to use to average",
    required=True,
)
def cli(csv, name, nuc_range):
    """
    A simple program that takes chemical mapping data with different titrating conditions
    and computes the Kd and hill coefficient for the titration.
    """
    setup_applevel_logger()
    log.info("Starting pytitration")
    log.info("CSV: %s", csv)
    log.info("Name: %s", name)
    log.info("Range: %s", nuc_range)
    df = pd.read_csv(csv)
    # make sure this is 1 indexed not 0 indexed
    num_range = [i + 1 for i in str_to_range(nuc_range)]
    data = []
    if "dir" not in df.columns:
        log.error("No 'dir' column in CSV file.")
        exit(1)
    if "conc" not in df.columns:
        log.error("No 'conc' column in CSV file.")
        exit(1)
    for i, row in df.iterrows():
        pickle_path = row["dir"] + "/output/BitVector_Files/mutation_histos.p"
        if not os.path.isfile(pickle_path):
            raise ValueError("No Pickel file found at: " + pickle_path)
        with open(pickle_path, "rb") as f:
            mut_histos = pickle.load(f)
        type_val = type(list(mut_histos.values())[0])
        if type_val == DREEMMutationHistogram:
            mut_histos = convert_dreem_mut_histos_to_mutation_histogram(mut_histos)
        if name not in mut_histos:
            raise ValueError("No data found for: " + name)
        cur_mut_histo = mut_histos[name].get_dict()
        pop_avg = get_pop_avg(cur_mut_histo)
        avg = np.average(sub_select_data(pop_avg, num_range))
        data.append([row["conc"], avg])
    df = pd.DataFrame(data, columns=["conc", "avg"])
    log.info("writing titration data to titration_data.csv")
    df.to_csv("titration_data.csv", index=False)
    pstart = [1, 1]
    norm_data = -normalize_data(np.array(df["avg"])) + 1
    pfit, perr = fit_bootstrap(pstart, df["conc"], norm_data, normalized_hill_equation)
    log.info("Kd = %f", pfit[0])
    log.info("n = %f", pfit[1])
    log.info("saving plot to titration_plot.png")
    generate_titration_plot(df)
    plt.savefig("titration_plot.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    cli()
