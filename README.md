# Estimate LW
Replace dhsvm forcing longwave with estimates from empirical approaches

In order to run this program you can run the following command:

python est_lw.py [dhsvm forcing] [elevation] [optional:output location]

dhsvm forcing: dhsvm forcing to correct longwave for
elevation: the elevation in meters of the forcing location
The output location can optionally be specified. If omitted, the original forcing will be overwritten.

By default, the above will run with all-sky option 7 and clear-sky option 13.
If a different method is preferred, they can be specified by running est_lw() instead of th default option, get_lw().

