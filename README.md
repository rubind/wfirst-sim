# wfirst-sim

A few notes about how to run things:

Two githubs:

https://github.com/rubind/wfirst-sim
https://github.com/rubind/wfirst-sim-data

export WFIRST_SIM_DATA=/path/to/wfirst-sim-data
export WFIRST_SIM=/path/to/wfirst-sim

The simulation package:

wfirst-sim/scripts/stan_cosmoSTEP1_simple_survey.py

Which takes a parameter file and the name of an output pickle file.

Then STEP2_Analytic_Fisher.py runs with the pickle file and makes a distance-modulus covariance matrix (comb_mat.fits).

FoM.py takes that comb_mat.fits and makes FoM values. We’ve been using FoM_0.26.
