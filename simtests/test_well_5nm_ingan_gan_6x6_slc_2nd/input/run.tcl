set options(-verbose)   1
set options(-slc)       "densities.dat"
set options(-nocache)   1
set options(-plot)      "bands:transitions:probabilities_ascii:optical"
set options(-model)     "kp6x6WZ"
set options(-grid)      "InGaN.asc.gz"
set options(-log)       "output.log"
set options(-krange)    "0.0001:1.5"
set options(-nk)        32
set options(-deg)       2
source "../../../input/tow.tcl"
