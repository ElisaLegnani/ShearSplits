#!/bin/bash           

export RUNNAME=test

export BASE=/global/cfs/cdirs/des/elisa/cosmosis-standard-library
source /dvs_ro/cfs/projectdirs/des/zuntz/cosmosis-global/setup-cosmosis3

export DATAFILE=/global/u1/e/elegnani/ShearSplits/cosmosis_runs/output_dir/DESshear.fits
export VALUESINCLUDE=/global/u1/e/elegnani/ShearSplits/cosmosis_runs/cosmosis_config/values_cosmo.ini
export OUT_DIR=/global/u1/e/elegnani/ShearSplits/cosmosis_runs/output_dir/

cosmosis /global/u1/e/elegnani/ShearSplits/cosmosis_runs/cosmosis_config/params.ini