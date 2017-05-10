---
layout: default
title: GPU Parallel Best Subset Regression
permalink: /
---
# Project Summary

The LSQ best subset regression algorithm can be used to identify accuracte linear models with a small number of variables. 
The performance of this algorithm is dominated by the cost of performing one inherently sequential kernel repeatedly. 
It is possible to parallelize this algorithm to be suitable for GPU computing by eliminating dependencies by repeating certain calculations. 

# Techincal Challenge

The original version of this algorithm was not designed to be run on a parallel computing platform. 
The LSQ algorithm iteratively reads in a row of data at a time, and updates a previously computed QR factorization result. 
This method sequentially updates other values using a Givens rotation approach. 
The challenge with this method is that the Givens rotation has a lot of sequential dependencies and does not lend itself well to being parallelized. 
In particular, there was one part of the algorithm that I beleived could be run completely independent of each other, but as I was implementing my parallel version, I realized that there was another dependency that I had overlooked. 
Another challenge is that it is not obvious how to compare the results from other state-of-the-art algorithms with the results from the LSQ algorithm.
Since the LSQ algorithm iteratively updates the QR factorization result, the error in residuals between this iterative approach and other methods are not comparable.
As such, I have not been able to compare my results against other high performance algorithms, and have only been able to compare the performance against a sequential single-core LSQ implementation. 

# Preliminary Results

I am still in the process of generating speedup plots. Currently, my initial implementation is accurate but about 4x slower than the sequential CPU code for fairly large problems. 
I am hoping to fix a few more bugs, and then improve my performance a bit more before the final deadline on Friday.

# Summary 

The final results on Friday are expected to show that it is possible to parallelize this sequential LSQ algorithm. 
Hopefully, the results will show that for smaller problems, the CPU sequential version of the LSQ algorithm will be faster than the GPU version becasue of the overhead of executing a GPU code and transferring data to the GPU from the CPU.
However, as the problem size increases, the GPU algorithm should begin to perform better than the CPU algorithm as the problem increases in size. 
 

# Updated Schedule

| Date            | Goals | Status |
|---|---|---|
| April 15, 2017  | Set up a test harness to compare the performance of different best subset regression implementations with some example problems. | Complete |
| April 23, 2017  | Complete an initial CUDA implementation of the LSQ algorithm, and confirm its accuracy by using the example problems created in the first week. | C++ version complete |
| April 26, 2017  | Extend harness to use other state-of-the-art linear least squares solvers. | Complete, but has issues |
| April 30, 2017  | Identify and experiment with different parts of the algorithm that can be parallelized. If time permits implement a hybrid CPU-GPU algorithm. | In progress |
| May 3, 2017     | Carry out experiments with each of the best subset regression implementations on well-conditioned problems. | Complete |
| May 7, 2017     | Carry out experiments with each of the best subset regression implementations on ill-conditioned problems. | Delayed for lack of time |
| May 10, 2017    | Resolve last minute issues and generate speedup graphs that will be used to show the effectiveness of the parallel GPU algorithm. | In progress |


## Please use the links above to access the project [Proposal](https://bsauk.github.io/GPU_LSQ/proposal) and [Checkpoint](https://bsauk.github.io/GPU_LSQ/checkpoint) reports.