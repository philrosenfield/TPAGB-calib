# Calibrating TP-AGB Models #

Attempting to refactor original TP-AGB calibration code to work on
HST data that isn't ANGST.


## Contributors ##
 * Phil Rosenfield


## Initial Preparations for Each New Data Set##
1. Prepare the data and run MATCH using HybridMC to get the SFH
2. Prepare template galaxy input for trilegal (if simulation mass isn't known, it can be passed to vary_the_sfh)
3. Make a function or two that define data color cuts, offsets, and select RGB and TPAGB stars

## Running ##
4: Make nsfhs of galaxy_input files and trilegal AMR files
5: Run trilegal for each, keeping only a subset of the data based on #3
6: Run the analysis routines
