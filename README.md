# pytitration

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

a short script to compute kd and hill coefficients from chemical mapping data

## Install

To install pytitration 

```shell
python -m pip install git+https://github.com/Yesselmanlab/pytitration
```


# How to run 

```shell
pytitration --help
Usage: pytitration [OPTIONS]

  A simple program that takes chemical mapping data with different titrating
  conditions and computes the Kd and hill coefficient for the titration.

Options:
  -c, --csv TEXT    Input csv file with dir and conc columns  [required]
  -n, --name TEXT   the construct name to use  [required]
  -r, --range TEXT  the nucleotides to use to average  [required]
  --help            Show this message and exit.
```

```shell 
pytitration --csv test/resources/test.csv --name atp_ttr_10 --range "110-111"              
PYTITRATION.CLI - INFO - Starting pytitration
PYTITRATION.CLI - INFO - CSV: test/resources/test.csv
PYTITRATION.CLI - INFO - Name: atp_ttr_10
PYTITRATION.CLI - INFO - Range: 110-111
PYTITRATION.CLI - INFO - writing titration data to titration_data.csv
PYTITRATION.CLI - INFO - Kd = 6.854161
PYTITRATION.CLI - INFO - n = 0.637045
PYTITRATION.CLI - INFO - saving plot to titration_plot.png
```