#!/bin/bash

M=5000
N=5000
LIMIT=25000
while [[ $M -le $LIMIT ]]; do
    echo "Execution Times for $M x $N" >> ./results/sizes5000.txt
    ./a.out $M $N 0 > problem.dat
    ./gpusub problem.dat $M $N 10 10 0 >> ./results/results5000.txt
    ((M+=5000))
done
