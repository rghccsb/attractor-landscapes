# Attractor Landscapes

Code and data accompanying the paper [Attractor Landscapes as a Model Selection Criterion in Data Poor Environments](https://doi.org/10.1101/2021.11.09.466986).

## Terms of Use

Prototype tools and models developed specifically under the scope of this project are freely available to the community for **non-commercial use and research purposes only**. These are available under a set of standard licensing terms protecting ownership by the Center and legislating not-for-profit use responsibilities and liabilities.

## How to Use

The models are provided in the `models` directory, with the network defined in `models/20200415_cv_MCM_Guimera.json` and the 19 model parameters defined in `models/20200415_cv_MCM_Guimera.solutions.out`. Additionally, the sqlite database `models/20200415_cv_MCM_Guimera.db` has these models as well as some of the attractors.

To run the Python code, a conda environment file is provided at `environment.yml`. You can create the environment with the dependencies by running `conda env create --file environment.yml`. After the environment is created, you activate it by running `conda activate attractor-landsacpes`. To discover the point attractors (attractors that have 1 state) run `python scripts/point_attractor_search.py`. To discover the cyclic attractors (attractors that have more than 1 state) run `python scripts/cyclic_attractor_search.py`. To discover the attractors with different drugs applied, run `python scripts/drug_attractor_search.py`. And finally, to discover the performance of the cyclic attractors, run `python scripts/cyclic_attractor_performance.py`.

To recreate the figures, you can run `python scripts_figures/generate_attractor_table.py` and `Rscript scripts_figures/significance_test.R`. The results are stored in the `figures` directory.
