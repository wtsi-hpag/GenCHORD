.. theory

#####################
Theory of Operation
#####################

When a genome is sequenced and aligned against a reference, a number of portions of the DNA can vary drastically vary from the reference - this can be due to normal variations in DNA between individuals, or due to substantial variations such as caused by genetic diseases or cancer. 

One potential structural variation can be that a region of DNA is either duplicated (potentially multiple times) or deleted from one or both copies of a chromosome, resulting in Copy Number Variations. When the size of the region which has been duplicated or deleted is shorter than the length of the reads used to perform the sequencing, this difference can be resolved by the alignment tool. If, however, the anomalous region is *larger* than the average read size, then the tools will instead register this as a variation in coverage, with only the (small number of) reads which cross the boundary between copies being evidence to the contrary.

The image below shows how an increase in coverage is evidence of a duplicated region with respect to the reference.

.. tikz:: A demonstration of coverage-aliased copy numbers
	
	\node[anchor=east] at (0,0.05) {True genome:};
	\fill (0,0) rectangle (2,0.1);
	\fill[blue] (2,0) rectangle (3,0.1);
	\fill[blue,pattern=north west lines, pattern color=blue] (3,0) rectangle (4,0.1);
	\fill (4,0) rectangle (6,0.1);
	\def\readheight{-1}
	\node[anchor=east] at (0,{\readheight-0.05}) {Reads:};
	\foreach[count=\i] \x in {0.3,1,1.4,1.8,2.2,2.5,2.8,3.3,3.7,4.1,4.7,5.1,5.3,5.6}
	{
		\fill[green!50!black] (\x,\readheight) rectangle ({\x+0.2},{\readheight -0.05});
	}
	\foreach[count=\i] \x in {0.2,0.6,1.1,1.5,2,2.35,2.75,3,3.4,3.8,4.5,5.1,5.8}
	{
		\fill[green!50!black] (\x,{\readheight-0.1}) rectangle ({\x+0.2},{\readheight -0.15});
	}
	\foreach[count=\i] \x in {0.3,1.2,1.5,1.8,2.1,2.4,2.9,3.3,4.1,5.3,5.5}
	{
		\fill[green!50!black] (\x,{\readheight+0.1}) rectangle ({\x+0.2},{\readheight +0.05});
	}
	\def\refheight{-2}
	\draw[->,red] (2,\readheight-0.15)--(2,\refheight);
	\draw[->,red] (4,\readheight-0.15)--(3,\refheight);
	\node[anchor=west,red] at (4,{(\refheight+\readheight-0.15)/2}) {Cannot distinguish};
	\node[anchor=east] at (0,{\refheight -0.05}) {Reference genome:};
	\fill (0,\refheight) rectangle (2,{\refheight -0.1});
	\fill[blue] (2,\refheight) rectangle (3,{\refheight -0.1});
	\fill (3,\refheight) rectangle (5,{\refheight -0.1});
	\def\areadheight{-3}
	\node[anchor=east] at (0,{\areadheight-0.05}) {Aligned Reads:};
	\foreach[count=\i] \x in {0.3,1,1.4,1.8,2.2,2.5,2.8,3.3,3.7,4.1,4.7,5.1,5.3,5.6}
	{
		\def\xoffset{0}
		\def\yoffset{0}
		\newdimen\pos
		\pos = \x cm
		\ifdim \pos > 3cm
			\def\xoffset{1}
			\ifdim \pos < 4cm
				\def\yoffset{0.3}
			\fi
		\fi
		\fill[red!50!black] ({\x-\xoffset},{\areadheight -\yoffset-0.1}) rectangle ({\x-\xoffset+0.2},{\areadheight -\yoffset-0.15});
	}
	\foreach[count=\i] \x in {0.2,0.6,1.1,1.5,2,2.35,2.75,3,3.4,3.8,4.5,5.1,5.8}
	{
		\def\xoffset{0}
		\def\yoffset{0}
		\newdimen\pos
		\pos = \x cm
		\ifdim \pos > 3cm
			\def\xoffset{1}
			\ifdim \pos < 4cm
				\def\yoffset{0.3}
			\fi
		\fi
		\fill[red!50!black] ({\x-\xoffset},{\areadheight -\yoffset}) rectangle ({\x-\xoffset+0.2},{\areadheight -\yoffset-0.05});
	}
	\foreach[count=\i] \x in {0.3,1.2,1.5,1.8,2.1,2.4,2.9,3.3,4.1,5.3,5.5}
	{
		\def\xoffset{0}
		\def\yoffset{0}
		\newdimen\pos
		\pos = \x cm
		\ifdim \pos > 3cm
			\def\xoffset{1}
			\ifdim \pos < 4cm
				\def\yoffset{0.3}
			\fi
		\fi
		\fill[red!50!black] ({\x-\xoffset},{\areadheight -\yoffset+0.1}) rectangle ({\x-\xoffset+0.2},{\areadheight -\yoffset+0.05});
	}
	\def\plotheight{-6}
	\node[anchor=east] at (0,{\plotheight-0.05+0.75}) {Coverage Data:};
	\draw[->](0.2,\plotheight)--(5,\plotheight);
	\node[anchor=north] at(2.5,\plotheight) {\tiny Index};
	\node[anchor=south,rotate=90] at(0.2,{\plotheight+0.75}) {\tiny Coverage};
	\draw[->](0.2,\plotheight)--++(0,1.5);
	\draw [red] plot [smooth] coordinates {(0.2,{\plotheight+0.5}) (0.4,{\plotheight+0.7}) (1,{\plotheight+0.4}) (1.3,{\plotheight+0.55}) (1.6,{\plotheight+0.7}) (1.7,{\plotheight+0.4})
		(2,{\plotheight+1.3})  (2.1,{\plotheight+1.1})  (2.4,{\plotheight+1.4})  (2.6,{\plotheight+1.45}) (2.8,{\plotheight+1.3})	(3,{\plotheight+1.3})  
		(3.1,{\plotheight+0.5}) (3.3,{\plotheight+0.6}) (3.5,{\plotheight+0.4}) (3.8,{\plotheight+0.7}) (4,{\plotheight+0.3}) (4.5,{\plotheight+0.6})  (5,{\plotheight+0.5})    
		};

********************
The Challenge
********************

	A given coverage file is, however, extremely noisy, even when using high-fidelity long read platforms. Whilst the human eye can pick out certain large scale patterns on a chromosomal scale, the exact criteria for what can be considered a genuine CNV and what is merely noise can be somewhat difficult to distinguish. 

	One potential solution would be to pass a basic smoothing filter over the coverage data, thereby "denoising" it. 

	However, a smoothing operation is in fact antithetical to what we are attempting to find, which is **sharp transition edges** between regions of different multiplicity. A smoothing kernel, by definition, obscures those edges and can entirely eliminate smaller scale transitions should the smoothing length scale be larger than the length scale over which multiplicity can vary.

	The Deforester tool is a method by which these transition edges can be located, without smoothing out all of the information of interest.

********************
Underlying Theory
********************

	The Deforester tool (so named because it reduces the "forest" of coverage data into a singular "tree") works on a series of successive Bayesian Inferences. Bayes' Theorem says that, given data :math:`D`, a hypothesis :math:`H_1` is better than the hypothesis :math:`H_2` if it meets the criteria:

	.. math::

		\frac{p(D | H_1) \text{prior}(H_1)}{p(D|H_2) \text{prior}(H_2)} > 1

	For convenience, we often work in logarithmic space, in which case the relationship is more easily expressed as:

	.. math::

		\mathfrak{p}(D|H_1) + \log(\text{prior}(H_1) > \mathfrak{p}(D|H_2) + \log(\text{prior}(H_2)

	Where :math:`\mathfrak{p}` is the log-probability of observing the data if the hypothesis were true, and :math:`\text{prior}(H)` is the *prior probability* we have that :math:`H` is true, before we looked at the data.

======================
Algorithmic Overview
======================

	As a brief summary of the algorithm used to compute the `tree` from the raw data:

	- First, compute the most likely value of the hyperparameters of the model (:math:`\gamma, \sigma, \nu`)
	- Iterate through the genome, one chromosome at a time
		- Iterate through each chromosome, one base at a time
		- At each base :math:`i`, compute the most likely harmonic for the range :math:`[i,i+L]`
			- If the harmonic is the previous base, increment :math:`i \to i+1` and repeat
			- If the harmonic is different:
				- Perform a Bayesian test to ensure validity
				- Refined the search window to locate the exact edge location
				- Increment :math:`i \to i +L` and repeat


===================
Statistical Model
===================

	We must therefore formulate a statistical model and aggregate the knowledge that we have about our system. The following are the basic principles upon which Deforester operates:

	- The genome is sampled to a given average depth across all chromosomes, termed the *fundamental frequency*, and denoted :math:`\nu`. We assume that we know this in advance (in practice, we compute the most likely value).
	- Every section of the genome has a given *multiplicity* or *resonance*, denoted :math:`q`, denoting the copy-number.
		- A region with an unresolved CNV is sampled to a depth of :math:`q\nu`.
		- Most portions of the human genome should have :math:`q = 2`, due to diploidy
	- The sampling process is noisy
		- For each sampling depth, the process is distributed according to a Poisson distribution with mean :math:`\lambda`
		- Variations in sampling mean that :math:`\lambda` varies across the genome. We assume this is random and Gamma-distributed such that the mean :math:`\bar{\lambda} = \nu`, with standard deviation :math:`\sigma`
		- There are then errors due to misalignment etc (this causes regions with true-coverage of zero to have small segments aligned to them, for example). We assume this error is Gaussian with standard deviation :math:`\gamma`
		- The noise is distributed according to the standard iid assumptions
	- Deviations in :math:`q` must be on a large scale (else they would have been resolved) - :math:`q` should vary only over regions greater than a lengthscale :math:`L`.
	- If index :math:`i` has :math:`q_i`, then we need a significant amount of evidence that :math:`q_{i+1} \neq q_i`.


	The result of the above is that the probability that the component at base :math:`i` was generated due to the :math:`q^\text{th}` resonance is given by:

	.. math::

		p(q_i = q | k_i, \nu, \sigma,\gamma) &= \sum_{k=0}^\infty N(k;k_i,\gamma) \int_0^\infty \text{Poisson}(k;\lambda) \Gamma(\lambda; q\nu, \sigma) \mathrm d \lambda
		\\
		&= \sum_{k=0}^\infty N(k;k_i,\gamma) \times NB\left(k; \frac{q^2 \nu^2}{\sigma^2}, \frac{q\nu}{q \nu + \sigma^2}\right)

	Where we have used the notation that :math:`N(k;\mu,\sigma)` is the usual normal distribution probability given mean :math:`\mu` and variance :math:`\sigma^2`, and :math:`\text{Poisson}(k;\lambda)` is the usual Poisson probability of observing an integer event :math:`k` given an average rate :math:`\lambda` and :math:`\Gamma(\lambda,\mu,\sigma)` is the Gamma distribution (expressed in terms of its first and second moments, instead of the standard notation). In the second line we used the property that the Gamma distribution is the conjugate prior of the Poisson distribution, with the result being the Negative Binomial distribution, :math:`NB(k;\mu,\sigma)`.

	This is a function which we can precompute for a suitable range of :math:`k` and :math:`q` to be reused across the entire genome. 

	**The quantity of interest for us is the value of *q*, which corresponds to the most likely multiplicity of any base, since *q* is constrained to be an integer: our goal is to find the most likely value of *q* for each base in the dataset**

	Since we are using the standard *iid* assumptions, the probability that a *sequence* of bases between :math:`i` and :math:`j` was generated on the :math:`q^\text{th}` harmonic (i.e. in a region with multiplicity :math:`q`), is:

	.. math::

		\mathfrak{p}(q_{i \to j} = q | \{k_n\}, \nu, \sigma,\gamma) = \sum_{n=i}^j \log\left(p(q_i = q | k_i, \nu, \sigma,\gamma)\right)

	At each base, we therefore compute :math:`\mathfrak{p}(i \to i + L)` for :math:`0 \leq q \leq q_\text{max}`, and find the value which is largest. The inclusion of :math:`L` ensures that we only identify transitions which are larger than some pre-specified minimum size, which helps remove random noise since we know from the specification of the problem that we should only see variations larger than the read size.

	If the maximum-likelihood value of *q* is *the same* than the current `running` value of *q* (i.e. the assigned *q* to the base :math:`i-1`), then we increment the counter by 1, and perform the test on the set :math:`i+1 \to i+L+1`: we did not find a transition. 

	If the maximum value of *q* is *different* than the `running` value, then we have identified a *transition candidate*. To see if this candidate is statistically valid, we perform a Bayesian test. The candidate is accepted if:

	.. math::

		\mathfrak{p}(q_{i \to i+L} = q_\text{new} ) > \mathfrak{p}(q_{i \to i+L} = q_\text{old} ) - \log(\text{Prior}(q_\text{old},q_\text{new}))

	The prior controls how much evidence we need before we change our minds (our prior belief being that subsequent *q* should be the same). We set our prior to be such that :math:`\log(\text{Prior}(q_\text{old},q_\text{new})) = \alpha`, i.e. a constant. This value is usually around -4, meaning we need the transition model to be 100 times more likely than the constant model.

	We then increment the counter by :math:`L` (since our requirement is that the model does not show variations smaller than this scale), and continue the search.

---------------------
Refined Search
---------------------

		Despite initially seeming like a good idea, the restriction that we move an entire block of length :math:`L` can lead to poor fits. Consider the following example sequence:


		``10  9  10  11  8  1000  999  1004  7  11  13  10``

		It seems clear that there is a transition between an average of 10, to an average of 1000, back to an average of 10. With :math:`\nu=5`, we would expect something along the lines of:
		
		``2  2  2  2  2  200  200  200  2  2  2  2``

		If, however, we analysed this sequence using the restrictions above with :math:`L=3` and , then we would find something like the following:

		``2 2 2 60 60 60 120 120 120 2 2 2``

		This is because we move the search window one base at a time (since we do not know when a transition might start -- and in reality :math:`L\gg 1`, so making jumps of length :math:`L` would wildly miss the beginning of transitions), and so we accidentally smooth over the transition edge, making the inference a disappointing compromise which meets neither set of data.

		Once a transition has been accepted, we therefore perform an additional check, moving the start and end of the window around by a number of bases (up to some pre-specified fraction of :math:`L`), searching for the exact location of the most likely beginning of the transition, and, given this more accurate model of where the transition begins, if there is a different value of :math:`q` which should be assigned.

		The end result is that it is indeed possible for the model to generate jumps from :math:`q_1 \to q_2 \to q_1` where the total length of the :math:`q_2` sequence is less than :math:`L`, seemingly changing the definition of :math:`L` from the `minimum transition size`: however, the model must first be *willing* to move a block of width :math:`L` for us to begin this more accurate search.