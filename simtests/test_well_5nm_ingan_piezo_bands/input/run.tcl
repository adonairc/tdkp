set options(-verbose)   1
set options(-strain)    1
set options(-piezo)     1
set options(-nocache)   1
set options(-plot)      "piezo:bandedges:bands:transitions:strain:potential"
set options(-model)     "kp6x6WZ"
set options(-grid)      "InGaN.asc.gz"
set options(-log)       "output.log"
set options(-krange)    "0.0001:1.5"
set options(-nk)        32
set options(-conf)      "conf.tcl"
set options(-bands)     1

source "../../../input/tow.tcl"
