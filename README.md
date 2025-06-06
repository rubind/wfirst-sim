# If you are looking for the distance-modulus covariance matrices from https://arxiv.org/abs/2506.04327, they are here:

https://github.com/rubind/roman_fisher_results

# wfirst-sim

A few notes about how to run things:

You will need two GitHubs:

https://github.com/rubind/wfirst-sim

and

https://github.com/rubind/wfirst-sim-data

after cloning:

```export WFIRST_SIM_DATA=/path/to/wfirst-sim-data```

```export WFIRST_SIM=/path/to/wfirst-sim```

The simulation package:

wfirst-sim/scripts/stan_cosmo/STEP1_simple_survey.py

Which takes a parameter file and the name of an output pickle file:

```python STEP1_simple_survey.py paramfile.csv survey.pickle```

And for making some summary plots:

```python STEP1A_plot_survey survey.pickle```

Then wfirst-sim/scripts/stan_cosmo/STEP2_Analytic_Fisher.py runs with the pickle file and makes a distance-modulus covariance matrix (comb_mat.fits).

```python STEP2_Analytic_Fisher.py survey.pickle```

wfirst-sim/scripts/stan_cosmo/FoM.py takes that comb_mat.fits and makes FoM values. We’ve been using FoM_0.26.

```python FoM.py */comb_mat*.fits```

wfirst-sim/scripts/stan_cosmo/STEP1_survey_grid.py generates sets of surveys (it's setup for our Slurm queue) so that's good to look at as well.
