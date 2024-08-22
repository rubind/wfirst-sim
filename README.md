# wfirst-sim

A few notes about how to run things:

Two githubs:

https://github.com/rubind/wfirst-sim
https://github.com/rubind/wfirst-sim-data

```export WFIRST_SIM_DATA=/path/to/wfirst-sim-data```
```export WFIRST_SIM=/path/to/wfirst-sim```

The simulation package:

wfirst-sim/scripts/stan_cosmoSTEP1_simple_survey.py

Which takes a parameter file and the name of an output pickle file:

```python STEP1_simple_survey.py paramfile.csv survey.pickle```

Then STEP2_Analytic_Fisher.py runs with the pickle file and makes a distance-modulus covariance matrix (comb_mat.fits).

```python STEP2_Analytic_Fisher.py survey.pickle```

FoM.py takes that comb_mat.fits and makes FoM values. We’ve been using FoM_0.26.

```python FoM.py */comb_mat*.fits```

STEP1_survey_grid.py generates sets of surveys (it's setup for our Slurm queue) so that's good to look at as well.
