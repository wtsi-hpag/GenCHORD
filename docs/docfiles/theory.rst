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
	\draw[->](0,\plotheight)--(5,\plotheight);
	\draw[->](0,\plotheight)--++(0,1.5);
	\draw [red] plot [smooth] coordinates {(0,{\plotheight+0.5}) (0.4,{\plotheight+0.7}) (1,{\plotheight+0.4}) (1.3,{\plotheight+0.55}) (1.6,{\plotheight+0.7}) (1.7,{\plotheight+0.4})
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

	Where :math:`\matfrak{p}` is the log-probability of observing the data if the hypothesis were true, and :math:`\text{prior}(H)` is the *prior probability* we have that :math:`H` is true, before we looked at the data.