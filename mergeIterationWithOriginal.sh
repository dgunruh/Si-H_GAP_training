#!/bin/bash

FILES=./Si:H_Extended_XYZ_structures/selectedStructures_round6/*.xyz
output=./GAP_fitting/Training_Data/trainingGAP_iterationStructures/additionalStructures_6.xyz
rm $output
for i in $FILES
do
	eval cat $i >> "${output}"
done

cp ./GAP_fitting/Training_Data/round5Iteration.xyz ./GAP_fitting/Training_Data/round6Iteration.xyz
cat ${output} >> ./GAP_fitting/Training_Data/round6Iteration.xyz
