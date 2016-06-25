Developer test scripts
======================

This folder contains scripts that are used for various testing purposes
during code development. None of the scripts actually needs to work, they
may actually be outdated, but could still be useful for further testing
if needed.

`benchmark_cta_analysis.py`
>  This script peforms a benchmark of the various CTA analysis types
>  (unbinned, binned, stacked), and if matplotlib is installed, creates
>  a plot of the benchmark results.

`check_models.py`
>  Checks a variety of models using an unbinned analysis.

`make_pull_at_sensitivity_limit.py`
>  This script creates pull distributions as function of energy
>  for a test source that is at the 5 sigma sensitivity limit for
>  50h of effective observing time. It takes as input file the
>  output file generated using cspull.
>  If the Python processing module is installed the script can make
>  use of multiple cores and runs in parallel.

`make_survey.py`
>  This script simulates a CTA galactic plane survey.

`simulate_events.py`
>  Simulates events and compares the simulation to the model.

