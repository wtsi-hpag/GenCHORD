from matplotlib import pyplot as pt
import math
import scipy
import numpy as np



def Poisson(ks,mu):
	if mu > 0:
		return np.exp(ks * math.log(mu) -mu - scipy.special.loggamma(ks+1))
	else:
		return (ks==0)
def NB(ks,mu,var):
	sigma = math.sqrt(var)
	if mu > 0:
		r = (mu/sigma)**2
		p = mu/(mu + sigma**2)

		loggams = scipy.special.loggamma(ks+r) - (scipy.special.loggamma(r) + scipy.special.loggamma(ks+1))
		return np.exp(loggams) * np.power(1-p,ks) * p**r
	else:
		return (ks==0)

def Full(ks,mu,sigma,lm,gamma):
	p1 = (1-gamma) * NB(ks,mu,sigma)

	p2 = gamma/lm * (ks <= lm) #* scipy.special.gammainc(ks+1,lm)

	return p1 + p2

mu = 30
i = 0
fig,axs = pt.subplots(2,sharex=True)
for mu in [30,0]:
	sigma = 10

	gamma = 1e-3
	N = 100
	ks = np.linspace(0,N,N+1)

	poiss = Poisson(ks,mu)

	axs[i].plot(ks,poiss,label="Poisson")

	c = 0
	styles = ['solid','dotted','dashed','dashdot']
	for sigma in [1,10,100]:
		nb = NB(ks,mu,sigma)
		leg = "NB, variance = " + str(sigma)
		if c== 0:
			ls, = axs[i].plot(ks,nb,label=leg)
		else:
			axs[i].plot(ks,nb,color=ls.get_color(),linestyle=styles[c],label=leg)
		c+=1

	c = 0
	sigma = 10
	for lm in [80,100,1000]:
		full = Full(ks,mu,sigma,lm,gamma)
		leg = "EB, variance=10, cutoff =" + str(lm)
		if c== 0:
			ls, = axs[i].plot(ks,full,label=leg)
		else:
			axs[i].plot(ks,full,color=ls.get_color(),linestyle=styles[c],label=leg)
		c+=1


	if i == 1:
		axs[i].legend()
	print("Poiss:",np.sum(poiss))
	print("NB:",np.sum(nb))
	# print("Full:",np.sum(full))
	axs[i].set_yscale('log')
	axs[i].set_ylim([1e-10,1])
	axs[i].set_title("Mean = " + str(mu))
	i+=1


fig.text(0.5, 0.04, "Coverage, k", ha='center')
fig.text(0.02, 0.5, "Probability of Occurence", va='center', rotation='vertical')
pt.draw()
pt.pause(0.02)
input("Enter to exit")
