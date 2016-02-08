# Simulating 1+1 dimensional directed stationary Eden (FPP) model

This simulates and visualizes certain statistics for the origin-rooted tree in the directed 1+1 dimensional directed stationary Eden (First Passage Percolation) model from the paper [Stationary Eden model on groups](http://arxiv.org/abs/1410.4944) co-authored with Eviatar B. Procaccia. The simulation is optimized in the sense that it is building the tree as a Markov process adapted to a certain filtration. Furthermore, it is performing only the computation over the edges in the filtration that lie in the boundary of the currently generated cluster. The complexity for a single simulation is O(Height^2) where Height is the (random) height of the tree (cutoff at MAX height), so the worst case complexity is O(MAX^2). Computing the average complexity is equivalent to resolving the open mathematical problem of computing the power decay exponent in the tail distribution of the Height. Assuming the conjectured value of -2/3 for this exponent, the average complexity would be O(MAX^(4/3)).


## Note:
* The simulation produces outputs into a file one line for each simulation. Each line consisting of four integers, height, maximal left displacement (min_i l(i)), maximal right displacement (max_i r(i)) and the maximal width (max_i(r(i)-l(i))) of the tree (here i-th level of the tree consists of nodes {(j,i): l(i) <= j <= r(i)}. File results_20kx50k contains the output of 50 thousand simulations each with height capped at 20 thousand using exponential distribution. File results_50kx100k-uniform contains the output of 100 thousand simulations each with height capped at 50 thousand using uniform distribution.
* In C files NUM denotes the number of simulations, MAX is the height cutoff, set DUB to be twice the value of NUM.
* SFPP.c simulates the (more natural) model with the exponentially distributed edge weights, while SFPP_uniform.c simulates the model with the uniformly (0,1) distributed edge weights (faster). The code is virtually identical, but has been separated into two files in the view of the potential (academic) audience. For other distributions use the same code with the inverse cumulative distribution function.
* Using gcc compile SFPP.c with the -lm option.

## Dependencies:
* Python 2.7; Pylab (all for plotting)


## Authors:
Tonci Antunovic