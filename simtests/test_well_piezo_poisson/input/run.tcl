###############################################
# test script for piezo / poisson calculation #
###############################################

set options(-verbose) 1
set options(-conf)    "conf.tcl"
set options(-model)   "kp6x6WZ"
set options(-grid)    "InGaN.asc.gz"
set options(-piezo)   1
set options(-scale_coords) 0.6
set options(-strain)  1
set options(-plot)    "strain:piezo_ascii:potential_ascii"

source "../../../input/tow.tcl"
