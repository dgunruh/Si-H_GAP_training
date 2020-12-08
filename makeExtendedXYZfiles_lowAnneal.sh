#!/bin/bash
round=3
lbound=11
ubound=15


label=amorph
config=amorph
for ((i=$lbound; i <= $ubound; i++)); do
	for ((j=1; j<=5; j++)); do
		input=${label}_lowAnneal${i}_${j}.out
		output=${label}_lowAnneal${i}_${j}.xyz
		python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
	done
done
