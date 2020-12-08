#!/bin/bash

FILES=./Si:H_Extended_XYZ_structures/Iterative_Training_round7/*.xyz
output=./GAP_fitting/Training_Data/usingGAP_iterationStructures/Iterative_Training_round7.xyz
rm $output
for i in $FILES
do
	eval cat $i >> "${output}"
done
