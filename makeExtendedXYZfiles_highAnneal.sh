#!/bin/bash
round=5
lbound=16
ubound=20


label=liquid
config=liquid
for ((i=$lbound; i <= $ubound; i++)); do
	for ((j=1; j<=3; j++)); do
		input=${label}_highAnneal${i}_${j}.out
		output=${label}_highAnneal${i}_${j}.xyz
		python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
	done
done
