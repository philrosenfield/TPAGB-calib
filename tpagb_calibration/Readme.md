# Calibrating TP-AGB Models #

Utilities to assess TP-AGB models using star formation histories and synthetic stellar populations

## Contributors ##
 * Phil Rosenfield

## Initial Preparations for Each New Data Set##
* Reduce the data
* Run MATCH using HybridMC to get the SFH (see plotting.interactive_match_cmdlimits and calcsfh_parallel)
* Prepare template galaxy input for TRILEGAL (if simulation mass isn't known, it can be passed to sfhs.vary_sfh)

## Running ##
* Run prepare_data.py in the match hmc directory to create input file for vary_sfh
* Run sfhs.vary_sfh and execute output script in trilegal directory
* Add ast test and compress files using analysis.add_asts
* Make lfs, cmds, and do stats tests with analysis.analyze

