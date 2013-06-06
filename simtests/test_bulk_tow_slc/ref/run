# ------------------------------------------------
# run.tcl
#
# for testing tow.tcl
# ------------------------------------------------

if { [info exists env(OMP_NUM_THREADS)] } {
  set wtr $env(OMP_NUM_THREADS)
  set env(OMP_NUM_THREADS) 1
}

set models [list "kp4x4" "kp6x6" "kp8x8"]
foreach m $models {
    if { [info exists options] } {
        unset options
    }
    set options(-verbose)   1
    set options(-slc)       "densities.dat"
    set options(-bulk)      "GaAs"
    set options(-verbose)   0
    set options(-nocache)   1
    set options(-conf)      "conf.tcl"
    set options(-noquit)    1
    set options(-plot)      "bands:transitions"
    set options(-model)     $m
    set options(-log)       [format "output_%s.log" $m]
    set options(-krange)    "0.0001:1.5"
    set options(-nk)        32
    # create config file
    set f [open "conf.tcl" "w"]
    puts $f "set vb_bands_ascii_file \"vb_disp_$m.dat\""
    puts $f "set cb_bands_ascii_file \"cb_disp_$m.dat\""
    puts $f "set optics_ascii_file   \"optics_$m.dat\""
    puts $f "set transitions_ascii_file \"transitions_$m.dat\""
    close $f
    exec cat conf_add.tcl >> conf.tcl
    source "../../../input/tow.tcl"

}


if { [info exists wtr] } {
  set env(OMP_NUM_THREADS) $wtr
}
