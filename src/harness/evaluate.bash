#!/bin/bash

M=1000
N=100
NLIM=500
LIMIT=20000
make clean
g++ makeMatrix.cpp 
make
while [[ $M -le $LIMIT ]]; do
    while [[ $N -le $NLIM ]]; do
	echo "Execution Times for $M x $N" >> results.txt
	./a.out $M $N 0 > problem.dat
	./gpusub problem.dat $M $N 10 10 0 >> results.txt
	((N+=100))
    done
    N=100
    ((M+=1000))
done
