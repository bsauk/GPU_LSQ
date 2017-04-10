# **GPU Parallel Unconstrained Linear Least Squares Algorithm**

Benjamin Sauk

## Summary

I will implement a GPU parallel unconstrained linear least squares algorithm based off of Alan Miller's algorithm designed for Fortran 90. I will compare the performance of solving best subset regression against the original Fortran implementation. I will also compare the performance of the GPU implementation against the performance of a more exhaustive way to solve the best subset regression problem.

## Background

Unconstrained linear least squares is used to develop linear models from data. There are many ways to solve the linear least squares problem that tradeoff computational complexity for numerical accuracy. In particular, there are many dense linear algebra libraries that can solve these problems with QR factorization such as: LAPACK for single core CPUs, PLASMA for multi-core CPUs, and MAGMA for hybrid CPU-GPU environments. However, these methods are not able to produce accurate solutions for very ill-conditioned matrices. The unconstrained linear least squares problem is used in a variety of disciplines such as signal processing, machine learning, and best subset regression. 

In particular, I am interested in solving linear least squares problems for best subset regression. Best subset regression is a method to determine the most accurate model given a large set of data. The problem in this context involves determining what subset of variables generates a the best linear model that matches the input data. This problem requires the solution of a large number of linear least squares problems and comparing the results to identify the best subset of variables.

There are multiple approaches that can be taken to solve the best subset regression problem. The first involves solving a linear least squares problem for every subset of variables from the original input data. This approach is time consuming, and relies on the performance of solving each instance of the linear least squares problem. The downsides to this approach are that as the number of variables in a problem increases, it becomes intractable to solve that many linear least squares problems, and that if more information is added, the problem needs to be completely resolved again. Another approach is to solve a single instance of the linear least squares problem with a method such as QR factorization, and then to add or remove rows or variables and then update or downdate the QR factorization solution with a new set of variables. LAPACK, PLASMA, and MAGMA all can be used in the former approach to quickly solve a linear least squares problem, and then find the best solution out of all of the subsets. Alternatively, LSQ.f90 relies on finding a solution to the linear least squares problem, and then updating or downdating the solution by adding or removing data from the matrix.

In the problem of best subset regression, there exists two different axis of parallelism that can be exploited. First, since a large number of linear least squares problems will need to be solved, it is possible to solve independent problems simultaneously. Also, there is the ability to parallelize each solution of the linear least squares problem by using Alan Miller's LSQ algorithm in places where there are independent loop operations. In the case of the LSQ.f90 algorithm, it should be possible to parallelize multiple update or downdate iterations simultaneously. 

## The Challenge

1. From examining the LSQ.f90 algorithm, there is a large portion of the code which is not currently parallelizable and would require a different strategy than what is currently implemented to speedup the solution to the best subset regression problem. I will also have to identify data dependencies that exist in the algorithm. 

2. The LSQ.f90 algorithm has potential for divergent execution, depending on the input data, and I think that the communication to computation ratio will be problem dependent. 

3. There will be a trade-off between how finely grained the least squares algorithm can be and how many independent linear least squares problems can be executed simultaneously on a NVIDIA GPU. I may have to use some type tuning rules or heuristics to identify how much I should parallelize each independent linear least squares problem.

## Resources

I will be starting with the LSQ.f90 algorithm created by Alan Miller. I will begin by implementing this code in CUDA and identify opportunties for independent calculations that can be run simultaneously in parallel. I will also start development on the GHC machines or the latedays clusters so that I have access to NVIDIA K40 GPUs to test the performance of my parallel algorithm after implementing it. I was currently planning on developing this algorithm for a single GPU implementation, but if time permits, it may be interesting to pursue a multi-GPU implementation if there is enough parallelism to exploit. 

## Goals and Deliverables
### Plan to Achieve

1. Implement the LSQ.f90 algorithm in a parallel framework in CUDA.

2. Compare the performance of solving best subset selection regression using a few different state-of-the-art linear least squares solvers. This way I will be able to compare the performance of using the update or downdate techniques in the LSQ.f90 algorithm against the performance of quickly solving each linear least squares problem in a best subset regression problem.

3. I will show speedup graphs that compare the performance of this GPU implementation with other implementations that I compare against to see how well this implementation does. I will show the performance of each solver over a range of different sized matrices so that I can observe how well the algorithm scales with problem size.

### Hope to Achieve

1. Implement a hybrid CPU-GPU LSQ algorithm that is able to increase performance by utilizing both types of processing units efficiently. 

2. Compare the accuracy of the GPU LSQ implementation against the accuracy of solving linear a large number of linear least squares problems. This is especially important in the case where a matrix is ill-conditioned and cannot be solved accurately with other techniques.

## Platform Choice

The reason why I have chosen to use a GPU and to use CUDA is every QR factorization update or downdate should be a simple independent process that can be handled by an individual thread or thread block. With a GPU I am able to update or downdate a lot of solutions simultaneously and determine the best subset regression in a hopefully reasonable amount of time.

## Schedule

| Date            | Goals |
|---|---|
| April 15, 2017  | Set up a script to compare the performance of different best subset regression implementations with some example problems. |
| April 23, 2017  | Complete an initial CUDA implementation of the LSQ algorithm, and confirm its accuracy by using the example problems created in the first week. |   
| April 30, 2017  | Identify and experiment with different parts of the algorithm that can be parallelized. If time permits implement a hybrid CPU-GPU algorithm. |    
| May 7, 2017     | Carry out experiments with each of the best subset regression implementations on the example well-conditioned example problems. If time permits, also experiment with ill-conditioned problems. |
| May 10, 2017    | Clean up last minute issues and generate speedup graphs that will be used to show the effectiveness of the parallel GPU algorithm. |

<!---
Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/bsauk/lsq.io/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.

--->