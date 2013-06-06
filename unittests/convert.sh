#!/bin/bash

for ff in `find . -name "*.grd*"`
do
	echo $ff | sed -e 's/\.grd.*/.asc/'
done


