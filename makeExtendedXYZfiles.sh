#!/bin/bash
round=2
lbound=6
ubound=10


label=amorph
config=amorph
for ((i=$lbound; i <= $ubound; i++)); do
	input=${label}_optimize$i.out
	output=${label}_optimize$i.xyz
	python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
done

label=interstitial
config=interstitial
for ((i=$lbound; i <= $ubound; i++)); do
	input=${label}_optimize$i.out
	output=${label}_optimize$i.xyz
	python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
done

label=vacancy
config=vacancy
for ((i=$lbound; i <= $ubound; i++)); do
	input=${label}_optimize$i.out
	output=${label}_optimize$i.xyz
	python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
done

label=divacancy
config=divacancy
for ((i=$lbound; i <= $ubound; i++)); do
	input=${label}_optimize$i.out
	output=${label}_optimize$i.xyz
	python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
done

label=liquid
config=liquid
for ((i=$lbound; i <= $ubound; i++)); do
	input=${label}_optimize$i.out
	output=${label}_optimize$i.xyz
	python3 Espresso-to-ExtendedXYZ.py $round $input $output $config
done

