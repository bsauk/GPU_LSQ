# **GPU Parallel Unconstrained Linear Least Squares Algorithm**

Benjamin Sauk

## Summary

I will implement a GPU parallel unconstrained linear least squares algorithm based off of Alan Miller's algorithm designed for Fortran 90. I will compare the performance of common linear least square algorithms against the performance of this algorithm for a range of different sized problems. 

## Background

Unconstrained linear least squares is used to develop linear models from data. There are many ways to solve the linear least squares problem that tradeoff computational complexity for numerical accuracy. In particular, there are many dense linear algebra libraries that can solve these problems with QR factorization such as: LAPACK for single core CPUs, PLASMA for multi-core CPUs, and MAGMA for hybrid CPU-GPU environments. However, these methods are not able to produce accurate solutions for very ill-conditioned matrices. The unconstrained linear least squares problem is used in a variety of disciplines such as signal processing, machine learning, and best subset regression.

In particular, I am interested in solving linear least squares problems for best subset regression. Best subset regression is a method to determine the most accurate model given a large set of data. The problem in this context involves determining what subset of variables generates a linear model that best matches the input data. This problem requires the solution of a large number of linear least squares problems and comparing the results to identify the best subset of variables.

In the problem of best subset regression, there exists two different axis of parallelism that can be exploited. First, since a large number of linear least squares problems will need to be solved, it is possible to solve independent problems simultaneously. Also, there should be the ability to parallelize each solution of the linear least squares problem by using Alan Miller's LSQ algorithm in places where there are independent loop operations.

## The Challenge



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