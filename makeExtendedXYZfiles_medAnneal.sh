#!/bin/bash
round=7
lbound=31
ubound=35


label=amorph
config=amorph
for ((i=$lbound; i <= $ubound; i++)); do
	for ((j=1; j<=6; j++)); do
		input=${label}_heating${i}_${j}.out
		output=${label}_heating${i}_${j}.xyz
		python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
	done
done
