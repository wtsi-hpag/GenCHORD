# Deforest

The Deforester is a tool for inferring the average copy-number of large contiguous blocks of sequences within a genome from the raw coverage data, using Bayesian methods to infer the underlying stochastic parameters and hence extract a moving average, which is constrained to integer multiples of a fundamental frequency.

This tool was designed for analysing Chromothripsis in conjunction with the [steppingStone tool](https://github.com/wtsi-hpag/steppingStone), but has applications outside of that usecase. 

The tool has been verified to work on data from both long and short read platforms.