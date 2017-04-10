# Summary

I will implement a GPU parallel unconstrained linear least squares algorithm based off of Alan Miller's algorithm designed for Fortran 90 [1]. I will compare the performance of
solving best subset regression problems against the original Fortran implementation. The performance of the algorithm will be compared against the original implementation and
the performance of a more exhaustive way to solve the best subset regression problem.

# Background

Unconstrained linear least squares is used to develop linear models from data. There are many ways to solve the linear least squares problem that trade-off computational
complexity for numerical accuracy. In particular, there are many dense linear algebra libraries that can solve these problems with QR factorization such as: LAPACK [2] for
single core CPUs, PLASMA [3] for multi-core CPUs, and MAGMA [4] for hybrid CPU-GPU environments. However, these methods are not able to produce accurate solutions for very
ill-conditioned matrices. 

In particular, I am interested in solving linear least squares problems for solving the best subset regression problem. Best subset regression is a method used to determine the
most accurate model given a large set of data. When generating a linear model, not all measured variables may have a direct impact on the value of some dependent variable. In
those cases, the best subset regression problem attempts to find the best subset of variables that can be used to predict an output response [5]. There are multiple approaches
that can be taken to solve the best subset regression problem. The first involves solving a linear least squares problem for every subset of variables from the original input
data. This approach is time consuming, and relies on having high performance linear least squares problem solvers. The downside to this approach is that as the number of
variables in increases, it becomes intractable to solve every linear least squares problem for all subsets of variables. Also if more information is added or removed, the
problem needs to be completely resolved again without using the information obtained from solving a similar problem. Another approach is to solve a single instance of the linear
least squares problem with a method such as QR factorization, and then to add or remove rows or variables and then update or downdate the QR factorization solution [6]. LAPACK,
PLASMA, and MAGMA all can be used in the former approach to quickly solve a linear least squares problem, and then find the best solution out of all of the potential solutions.
Alternatively, LSQ.f90 relies on finding a solution to the linear least squares problem, and then updating or downdating the QR factorization solution to test the accuracy of
all subset models.

There are a few different ways to parallelize the best subset regression problem. By using Alan Miller's approach, it should be possible to parallelize multiple indpendent QR
factorization updates or downdates simultaneously. The solution to each linear least squares problem can also be parallelized to improve the performance of the most critical
part of each iteration in the algorithm.  

# The Challenge

1. From examining the LSQ.f90 algorithm, there is a large portion of the code which is not currently parallelizable and would require a different strategy than what is currently
implemented to speedup the solution to the best subset regression problem. I will also have to identify data dependencies that exist in the algorithm. 

2. The LSQ.f90 algorithm has potential for divergent execution, depending on the input data, and I think that the communication to computation ratio will be problem dependent,
and poor for matrices with a much larger number of rows than columns. 

3. There will be a trade-off between how finely grained the least squares algorithm can be and how many independent linear least squares problems can be executed simultaneously
on a NVIDIA GPU. I may have to use some tuning rules or heuristics to identify how much I should parallelize each independent linear least squares problem.

# Resources

I will be starting with the LSQ.f90 algorithm created by Alan Miller [1]. First, I will implement this code in CUDA and identify opportunties for independent calculations that
can be run simultaneously in parallel. I will start development on the latedays clusters so that I have access to NVIDIA K40 GPUs to test the performance of my parallel
algorithm after implementing it. I am planning on developing this algorithm for a single GPU implementation, but if time permits, it may be interesting to pursue either a hybrid
CPU-GPU or a multi-GPU implementation if there is enough parallelism to exploit. 

# Goals and Deliverables
## Plan to Achieve

1. Implement the LSQ.f90 algorithm in a parallel framework in CUDA. Also experiment with different update or downdate sizes.

2. Compare the performance of solving the best subset selection regression using a few different state-of-the-art linear least squares solvers.

3. I will show speedup graphs that compare the performance of this GPU implementation with other implementations that I compare against to see how well this implementation does.
By comparing the solvers over a range of different sized problems, I will demonstrate how the algorithm scales with problem size. 

## Hope to Achieve

1. Implement a hybrid CPU-GPU LSQ algorithm that is able to increase performance by utilizing both types of processing units efficiently. 

2. Compare the accuracy of the GPU LSQ implementation against the accuracy of solving a large number of linear least squares problems. This is especially important in the case
where a matrix is ill-conditioned and cannot be solved accurately with other techniques.

# Platform Choice

The reason why I have chosen to use a GPU and to use CUDA is that every QR factorization update or downdate should be a simple independent process that can be handled by an
individual thread or thread block. With a GPU I should be able to update or downdate a large number of solutions simultaneously and determine the best subset regression in a
reasonable amount of time.

# Schedule

| Date            | Goals |
|---|---|
| April 15, 2017  | Set up a test harness to compare the performance of different best subset regression implementations with some example problems. |
| April 23, 2017  | Complete an initial CUDA implementation of the LSQ algorithm, and confirm its accuracy by using the example problems created in the first week. |   
| April 30, 2017  | Identify and experiment with different parts of the algorithm that can be parallelized. If time permits implement a hybrid CPU-GPU algorithm. |    
| May 7, 2017     | Carry out experiments with each of the best subset regression implementations on the well-conditioned example problems. If time permits, also experiment with
ill-conditioned problems. |
| May 10, 2017    | Resolve last minute issues and generate speedup graphs that will be used to show the effectiveness of the parallel GPU algorithm. |


# References

[1] A. Miller. LSQ.f90, Current as of 10, April, 2017. jblevins.org/mirror/amiller/.

[2] E. Anderson, Z. Bai, C. Bischof, L. S. Blackford, J. Demmel, J. Dongarra, J. Du Croz, S. Hammarling, A. Greenbaum, A. McKenney, and D. Sorensen. LAPACK users’ guide (third ed.). Society for Industrial and Applied Mathematics, Philadelphia, PA, USA, 1999.

[3] B. Hadri, H. Ltaief, E. Agullo, and J. Dongarra. Tile QR factorization with parallel panel processing for multicore architectures. In Parallel & Distributed Processing (IPDPS), 2010 IEEE International Symposium on, pages 1–10, 2010.

[4] S. Tomov, J. Dongarra, and M. Baboulin. Towards dense linear algebra for hybrid GPU accelerated manycore systems. Parallel Computing, 36:232–240, 2010.

[5] http://support.minitab.com/en-us/minitab/17/topic-library/modeling-statistics/regression-and-correlation/basics/basics-of-best-subsets-regression/

[6] R. Andrew, N. Dingle. Implementing QR factorization updating algorithms on GPUs. Parallel Computing, 40:161-172, 2014. 