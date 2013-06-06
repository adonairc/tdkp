#!/bin/bash

if [ $1 != "" ]; then
  if [ ! -d $1 ]; then
    mkdir $1 $1/input $1/run $1/ref $1/misc
    
    echo "{" > $1/testdef
    echo "  file = ," >> $1/testdef
    echo "  tag  = ," >> $1/testdef
    echo "  col  = ," >> $1/testdef
    echo "  rel  = ," >> $1/testdef
    echo "  abs  =  " >> $1/testdef
    echo "}" >> $1/testdef
  else
    echo "error, $1 already exists"
  fi
fi
