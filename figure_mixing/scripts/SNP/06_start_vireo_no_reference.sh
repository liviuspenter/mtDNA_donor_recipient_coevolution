#!/bin/bash

for exp in 1 5 10 50 100 500 1000
do
	sbatch ./mixing.vireo.no.reference.sh $exp 
done

