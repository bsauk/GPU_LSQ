---
layout: default
title: Checkpoint
permalink: /checkpoint/
---
# Progress

For the last few weeks, my effort has been directed at translating the code I want to parallelize from Fortran to C++.
I went through the trouble of converting the LSQ algorithm as well as methods for performing best subset regression to C++ from Fortran for a few reasons.
First, I wanted to understand the algorithm, and thought that translating it from one language to another would help me to understand it.
Additionally, I wanted to translate it to C++ so that the code would be easier to convert to CUDA code when I begin parallelizing it.
The C++ code will serve as a comparison to the CUDA code to determine how much faster the parallel implementation is than the sequential version of the algorithm.
I recently finished up the code translation and now have an operating best subset regression algorithm implemented in C++.

I have also set up a test harness that I will be able to use to compare the performance of the different best subset regression implementations.
This is a C++ file that will measure the execution time of each of the algorithms and compare the performance results.
Additionally, I have the ability to compare the results from the parallel algorithm with the sequential version to ensure that accuracy is maintained.
I also have created a few test problems that I will use to determine the speedup for using the parallel algorithm on different sized problems.
Currently the matrices I have generated are well-conditioned matrices, but if time permits, I will also test some more ill-conditioned matrices.

# Revised Project Schedule

I am currently behind in the schedule that I created during the project proposal.
When I begina, I had forgotten that in addition to translating the LSQ algorithm into C++, I would also have to create a best subset regression framework in C++.
This required me to translate more Fortran code into C++, and took longer than I had previously anticipated.
However, since the code is all in C++, it should not be too much hassle to implement at least a sequential CUDA version so I will be able to spend the next few weeks parallelizing my code.
I still believe that I will be able to deliver on most of my goals, however, I am not sure if I will be able to achieve my hope to achieve goals.
The only goal I am concerned about, involves comparing the performance of the best subset regression algorithm against a few different state-of-the-art linear least squares solvers.
This goal is more challenging than I had expected, because the output from the LSQ algorithm is different than the output of other state-of-the-art least square solvers so I will have to generate other data structures to compare their performance directly on the best subset regression problem.
While this goal is still achievable, it will depend on how long it takes me to parallelize my LSQ algorithm. 

## Parallelism Competition Goals

1. Implement the LSQ.f90 algorithm in a parallel framework in CUDA.
2. Hopefully, compare the performance of solving the best subset regression problem using a few different state-of-the-art linear least squares solvers.
3. Show speedup graphs that compare the performance of this GPU implementation with other implementations.
4. Attempt to develop a hybrid LSQ algorithm.
5. Compare the accuracy of the GPU LSQ implementation against the accuracy of other methods to solve QR factorization. 

# Parallelism Competition

My goals for the parallelism compettion are still the same, I am going to show graphs that demonstrate the benefit of using this parallel algorithm for best subset regression.
The graph will compare the performance of this algorithm against the performance of other state-of-the-art algorithms.
Where the speedups will be shown over a wide range of problems to see how the performance of each algorithm scales with problem size.

# Concerns

At this point, I need to spend some time looking through the LSQ algorithm more to identify opportunities for parallelism.
The code seems inherently sequential, but I believe there are some locations that can be run in parallel.
My concern currently is that there will not be enough parallelism at those points to generate a significant speedup.
However, if that is the case, I will be able to still explore my other axis of parallelism, speeduping the best subset regression algorith, and not just speeding up each instance of the linear least squares problem solver.
Overall, there may be some difficulties moving forward, but I believe I should be able to deliver on my goals by the parallelism competition.

# Updated Schedule

| Date            | Goals | Status |
|---|---|---|
| April 15, 2017  | Set up a test harness to compare the performance of different best subset regression implementations with some example problems. | Complete |
| April 23, 2017  | Complete an initial CUDA implementation of the LSQ algorithm, and confirm its accuracy by using the example problems created in the first week. | C++ version complete |
| April 26, 2017  | Extend harness to use other state-of-the-art linear least squares solvers. |  |
| April 30, 2017  | Identify and experiment with different parts of the algorithm that can be parallelized. If time permits implement a hybrid CPU-GPU algorithm. |  |
| May 3, 2017     | Carry out experiments with each of the best subset regression implementations on well-conditioned problems. |  |
| May 7, 2017     | Carry out experiments with each of the best subset regression implementations on ill-conditioned problems. |  |
| May 10, 2017    | Resolve last minute issues and generate speedup graphs that will be used to show the effectiveness of the parallel GPU algorithm. |  |
