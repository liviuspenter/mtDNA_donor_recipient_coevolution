#!/bin/bash

for exp in 1 5 10 50 100 500 1000
do
	sbatch ./mixing.souporcell.with.reference.sh $exp
done

