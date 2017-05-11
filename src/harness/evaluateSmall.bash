#!/bin/bash

M=500
N=100
NLIM=500
LIMIT=20500
make clean
g++ makeMatrix.cpp 
make
while [[ $N -le $NLIM ]]; do
    while [[ $M -le $LIMIT ]]; do
	echo "Execution Times for $M x $N" >> ./results/sizesSmall.txt
	./a.out $M $N 0 > problem.dat
	./gpusub problem.dat $M $N 10 10 0 >> ./results/resultsSmall.txt
	((M+=5000))
    done
    echo " " >> ./results/resultsSmall.txt
    M=500
    ((N+=100))
done
