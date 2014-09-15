import numpy as np
import scipy.special as ss
import pylab

import pymixture


np.random.seed(2)
train_size=1000000
s = pymixture.Mixture(train_size=train_size)
lpi = -25
figsize = (8 * 1.5, 6 * 1.5)

prefix = "mix-ic-"; ic = s.first_proposal(); L = np.log(2.38**2 / 64); skip = 0
#prefix = "mix-mle-"; ic = s.positions; L = -13.8; skip = 0
#prefix = "mix-100-"; ic = s.first_proposal(); L = np.log(2.38**2 / 64); skip = 100

fs = 20

pylab.ion()
#xt = [6, 60, 600, 6000, 60000, 600000]
xt = [100, 1000, 10000, 100000, 1000000]
xl = [str(i) for i in xt]

colors = ['r', 'orange', 'y', 'g', 'c', 'b', 'purple', 'magenta']
colors = colors * 1
nc = len(colors)
log_p = np.zeros((nc, train_size))
log_q = np.zeros((nc, train_size))

theta = ic
for j in range(skip):
	proposal = s.next_proposal(theta, L)
	lp = s.log_posterior(theta)
	lq = s.log_posterior(proposal)
	print j	
	if (lq - lp) > np.log(np.random.random()):
		accept = 1.0
		print 'accept'
		theta = proposal.copy()
	else:
		accept = 0.0
		print 'reject'
	L += 1.0 / (j + 1)**0.5 * (accept - 0.234)

for j in range(nc):
	proposal = s.next_proposal(theta, L)
	
	log_p[j] = s.log_contributions(theta)
	log_q[j] = s.log_contributions(proposal)
	print j
	
	if (log_q[j].sum() - log_p[j].sum()) > np.log(np.random.random()):
		accept = 1.0
		print 'accept'
		theta = proposal.copy()
	else:
		accept = 0.0
		print 'reject'
	
	L += 1.0 / (j + 1)**0.5 * (accept - 0.234)

mhr = log_p.sum(axis=1) - log_q.sum(axis=1)
ind = mhr.argsort()
mhr = mhr[ind]
log_p = log_p[ind]
log_q = log_q[ind]
pos = np.arange(1, train_size+1)

pylab.figure(2, figsize=figsize)
pylab.clf()

for j in range(nc):
	r = log_q[j] - log_p[j]
	sum_r = np.cumsum(r)
	estimate = sum_r / np.arange(1, train_size + 1) * train_size
	#error = exact[j] - estimate

	var_r = np.cumsum(r**2) / pos - (sum_r / pos)**2
	ebar = (train_size - pos) * (1. / np.sqrt(pos) + 1. / np.sqrt(train_size - pos)) * np.sqrt(var_r)
	ebar[-1] = 10**-12

	# rho = cov / var_p / var_q = (exp_pq - exp_p * exp_q) / var_p / var_q
	var_p = np.cumsum(log_p[j]**2) / pos - (np.cumsum(log_p[j]) / pos)**2
	var_q = np.cumsum(log_q[j]**2) / pos - (np.cumsum(log_q[j]) / pos)**2
	cov_pq = 0.99 * np.sqrt(var_p * var_q)
	var_r_approx = var_p + var_q - 2 * cov_pq
	ebar_approx = (train_size - pos) * (1. / np.sqrt(pos) + 1. / np.sqrt(train_size - pos)) * np.sqrt(var_r_approx)
	ebar_approx[-1] = 10**-12
	
	#pylab.figure(1)
	#pylab.loglog(np.arange(1, train_size), np.abs(error[:-1]), 'b', alpha=0.2, linewidth=2)

	ind = np.cast[int](np.logspace(0, np.log10(train_size - 1), 1000))
	ind = list(set(ind)) + [0]
	ind.sort()
	
	pylab.figure(2)
	pylab.subplot(2,1,1)
	pylab.semilogx(pos[ind], estimate[ind], color=colors[j], linewidth=2)
	#pylab.semilogx(pos, estimate + ebar, color=colors[j], alpha=0.2, linewidth=1)
	#pylab.semilogx(pos, estimate - ebar, color=colors[j], alpha=0.2, linewidth=1)
	
	pylab.fill_between(pos[ind], estimate[ind] + ebar[ind], estimate[ind] - ebar[ind], color=colors[j], alpha=0.2)

	pylab.subplot(2,1,2)
	branch_1 = 0.5 + 0.5 * ss.erf((estimate - np.log(np.random.random())) / np.sqrt(2) / ebar)
	#branch_2 = 0.5 + 0.5 * ss.erf((estimate) / np.sqrt(2) / ebar)
	pylab.semilogx(pos[ind], branch_1[ind], '-', color=colors[j], linewidth=2)
	#pylab.semilogx(pos, branch_2, ':', color=colors[j], linewidth=2)


	
"""
pylab.figure(3)
pos = np.arange(1, train_size + 1)

sum_p = np.cumsum(log_p[j]) 
var_p = (np.cumsum(log_p[j]**2) - sum_p**2 / pos) / pos
estimate_p = sum_p / pos * train_size
error_p = np.sqrt(var_p * train_size)

sum_q = np.cumsum(log_q[j]) 
var_q = (np.cumsum(log_q[j]**2) - sum_q**2 / pos) / pos
estimate_q = sum_q / pos * train_size
error_q = np.sqrt(var_q * train_size)

r = log_p[j] - log_q[j]
sum_r = np.cumsum(r)
var_r = (np.cumsum(r**2) - sum_r**2 / pos) / pos
estimate_r = sum_r / pos * train_size

#float(N - n) * (1. / np.sqrt(n) + 1. / np.sqrt(N - n)) * s_n
#error_r = np.sqrt(var_r * train_size)

error_r = (train_size - pos) * (1. / np.sqrt(pos) + 1. / np.sqrt(train_size - pos)) * np.sqrt(var_r)

pylab.semilogx(pos, estimate_r, 'b')
#pylab.semilogx(pos, estimate_p - estimate_q + error_p + error_q, 'g')
#pylab.semilogx(pos, estimate_p - estimate_q - error_p - error_q, 'g')
pylab.semilogx(pos, estimate_r + error_r, 'r')
pylab.semilogx(pos, estimate_r - error_r, 'r')
"""

"""
pylab.figure(1)
pylab.xticks([5, 50, 500, 5000, 50000], ['5', '50', '500', '5000', '50000'], fontsize=fs)
pylab.yticks(fontsize=fs)
a = list(pylab.axis())
a[1] = train_size
pylab.axis(a)
pylab.title('Absolute error', fontsize=fs)	
pylab.xlabel('N', fontsize=fs)
pylab.draw()
"""

pylab.figure(2)
pylab.subplot(2,1,1)
pylab.semilogx([1, train_size], [1, 1], 'k:', linewidth=2)
pylab.xticks(xt, ['$\\mathdefault{10^%d}$' % np.log10(tt) for tt in xt], fontsize=fs)
(a, b) = pylab.yticks(fontsize=fs)
if (prefix == "mix-ic-"):
	yt = [-400000,  -200000,        0,   200000,  400000]
	ytn = ['$\\mathdefault{-4 \\times 10^{5}}$', '$\\mathdefault{-2 \\times 10^{5}}$', '0', 
		   '$\\mathdefault{2 \\times 10^{5}}$', '$\\mathdefault{4 \\times 10^{5}}$']
	pylab.yticks(yt, ytn)

#pylab.yticks(a[::2], fontsize=fs)
#pylab.yticks(a[1::2], fontsize=fs)
a = list(pylab.axis())
a[0] = xt[0]
a[1] = xt[-1]
if (prefix == 'mix-ic-'):
	a[2] = -500000
	a[3] = 500000

if (prefix == 'mix-mle-'):
	a[2] = -1000
	a[3] = 1000

pylab.axis(a)
pylab.ylabel("log(MH ratio) estimate", fontsize=fs+2)	
#pylab.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pylab.xlabel('', fontsize=fs)
pylab.title('Model of the MH ratio and acceptance probability\n', fontsize=fs+2)

pylab.subplot(2,1,2)
pylab.xticks(xt, xl, fontsize=fs)
pylab.yticks(fontsize=fs)
pylab.semilogx([2, train_size], [0.5, 0.5], 'k:', linewidth=2)
a = list(pylab.axis())
eps = 0.01
a[0] = xt[0]
a[1] = xt[-1]
a[2] = 0 - eps
a[3] = 1 + eps
pylab.axis(a)
pylab.xlabel('subsample size $(m)$', fontsize=fs+2)
pylab.ylabel('predictor $(\\psi_{\\rho}^{(m)})$', fontsize=fs+2)
pylab.draw()



pylab.savefig('../figs/%straces.pdf' % prefix)
#pylab.savefig('../figs/%straces-thick.pdf' % prefix)

"""
log_p = np.zeros(ndim)
perm = np.random.permutation(ndim)
#perm = np.arange(ndim)
for i in range(1, ndim + 1):
	M = np.dot(s.data[:, perm[:i]], theta[perm[:i], :]) / i * ndim
	log_p[i-1] = (M[np.arange(train_size), s.response] - np.log(np.exp(M).sum(axis=1))).sum()
	print i, log_p[i-1]

pylab.figure(1)
pylab.plot(log_p)
#pylab.axis(a)
pylab.draw()
"""

pylab.figure(4, figsize=figsize)
pylab.clf()
p = np.cumsum(log_p[0]) / pos * pos[-1]
q = np.cumsum(log_q[0]) / pos * pos[-1]
r = log_p[0] - log_q[0]
var_r = (np.cumsum(r**2) - sum_r**2 / pos) / pos
ebar = np.sqrt(pos[-1] * (pos[-1] - pos) / pos * var_r) 
ebar = (train_size - pos) * (1. / np.sqrt(pos) + 1. / np.sqrt(train_size - pos)) * np.sqrt(var_r)
pylab.subplot(3, 1, 1)
pylab.title('Example trajectories of estimates and error\n', fontsize=fs)
pylab.semilogx(pos[ind][1:], p[ind][1:], 'b', linewidth=2)
pylab.semilogx(pos[ind][1:], q[ind][1:], 'r:', linewidth=3)
pylab.ylabel('log(post) estimate', fontsize=fs)
(a, b) = pylab.xticks()
pylab.xticks(a[1:], fontsize=fs)
a = list(pylab.axis())
a[0] = xt[0]
a[1] = train_size
if (prefix == 'mix-ic-'):
	a[2] = -6.6*10**6
	a[3] = -5.4*10**6
	pylab.legend(["current state, $\\theta$", "proposal, $\\theta'$"], loc='upper right', fontsize=fs-2)

if (prefix == 'mix-mle-'):
	a[2] = -3.5*10**6
	a[3] = -3.8*10**6
	pylab.legend(["current state, $\\theta$", "proposal, $\\theta'$"], loc='lower right', fontsize=fs-2)

pylab.axis(a)
(a, b) = pylab.yticks()
if (prefix == 'mix-ic-'):
	a = [-6.5*10**6, -6*10**6, -5.5*10**6]

if (prefix == 'mix-mle-'):
	a = [-3.8*10**6, -3.7*10**6, -3.6*10**6, -3.5*10**6]

b = ['$\\mathdefault{%1.1f \\times 10^6}$' % (u/10**6) for u in a]
pylab.yticks(a, b, fontsize=fs-6)

pylab.subplot(3, 1, 2)
pylab.semilogx(pos[ind][1:], (q-p)[ind][1:], 'b', linewidth=2)
jj = np.invert(np.isnan(ebar[ind]))
pylab.fill_between(pos[ind][jj], (q-p)[ind][jj] + 2*ebar[ind][jj], (q-p)[ind][jj] - 2*ebar[ind][jj], color='b', alpha=0.2)
pylab.fill_between(pos[ind][jj], (q-p)[ind][jj] + ebar[ind][jj], (q-p)[ind][jj] - ebar[ind][jj], color='b', alpha=0.2)
pylab.semilogx([1, train_size], [1, 1], 'k:', linewidth=2)
(a, b) = pylab.xticks()
pylab.xticks(a[1:], fontsize=fs)
(a, b) = pylab.yticks()
if (prefix == 'mix-ic-'):
	a = [-10**5, 0, 10**5]
	b = ['$\\mathdefault{%d \\times 10^5}$' % (u/10**5) for u in a]
	b[1] = 0
	pylab.yticks(a, b, fontsize=fs-6)

if (prefix == 'mix-mle-'):
	a = [-500, 0, 500]
	pylab.yticks(a, fontsize=fs-6)

a = list(pylab.axis())
a[0] = xt[0]
a[1] = train_size
if (prefix == 'mix-ic-'):
	a[2] = -2*10**5
	a[3] = 2*10**5

if (prefix == 'mix-mle-'):
	a[2] = -800
	a[3] = 800

pylab.axis(a)
pylab.ylabel("log(MH ratio)\n", fontsize=fs-2)
pylab.legend(('estimate',), loc='upper right', fontsize=fs-2)

pylab.subplot(3, 1, 3)
pylab.loglog(pos[ind][1:-1], (np.abs(((q-p)[-1] - (q-p))[ind][1:-1])/(q-p)[-1]), 'b', linewidth=2)
pylab.loglog(pos[ind][1:-1], (np.abs(ebar[ind][1:-1]/(q-p)[-1])), 'g--', linewidth=3)
pylab.loglog(pos[ind][1:-1], (np.abs(2 * ebar[ind][1:-1]/(q-p)[-1])), 'c:', linewidth=3)
(a, b) = pylab.xticks()
pylab.xticks(a[1:], fontsize=fs)
a = list(pylab.axis())
a[0] = xt[0]
a[1] = train_size
pylab.axis(a)
(a, b) = pylab.yticks()
pylab.yticks(a[::2][1:-1], fontsize=fs-4)

pylab.xlabel('subsample size $(m)$', fontsize=fs)
pylab.ylabel("|relative error|\n(log scale)\n", fontsize=fs-2)	
pylab.legend(('actual', 'model $(1\\sigma)$', 'model $(2\\sigma)$'), loc='lower left', fontsize=fs-2, ncol=3)
pylab.savefig('../figs/%straces-pair.pdf' % prefix)

pylab.figure(5, figsize=figsize)
pylab.clf()
pylab.plot(log_p[0][ind], log_q[0][ind], 'b.')
pylab.plot([lpi, 0], [lpi, 0], 'k:', linewidth=2)
pylab.xlabel("\n$\\log$ $p(x_i | \\theta)$", fontsize=fs+2)
pylab.ylabel("$\\log$ $p(x_i | \\theta')$", fontsize=fs+2)
pylab.xticks(fontsize=fs)
pylab.yticks(fontsize=fs)
cc = np.corrcoef(log_p[0], log_q[0])[0, 1]
pylab.text(lpi * 0.95, lpi * 0.15, '$N$ = %d\ncorrelation > %2.6f' % (train_size, cc), fontsize=fs)
pylab.axis([lpi, 0, lpi, 0])
pylab.title('Likelihoods at the proposal versus the current state', fontsize=fs)
pylab.savefig('../figs/%scorrelation.pdf' % prefix)

pylab.figure(7, figsize=figsize)
pylab.clf()
pylab.title('Traces for different permutations of the data\n', fontsize=fs)
r = []
for i in range(nc):
	lp = log_p[0].copy()
	lp = lp[np.random.permutation(len(lp))]
	r += [np.cumsum(lp) / pos * pos[-1]]

ii = np.argsort([j[100] for j in r])[::-1]
for (i,c) in zip(ii, colors):
	pylab.semilogx(pos[ind][1:], r[i][ind][1:], c, linewidth=2)
	
(a, b) = pylab.xticks()
pylab.xticks(a[1:], fontsize=fs)
pylab.ylabel('log(posterior) estimate', fontsize=fs)
pylab.xlabel('subsample size $(m)$', fontsize=fs)
if (prefix == 'mix-ic-'):
	yt = [-6.6*10**6, -6.4*10**6, -6.2*10**6, -6.0*10**6, -5.8*10**6]

if (prefix == 'mix-mle-'):
	yt = [-3.9*10**6, -3.8*10**6, -3.7*10**6, -3.6*10**6, -3.5*10**6]

b = ['$\\mathdefault{%1.1f \\times 10^6}$' % (u/10**6) for u in yt]
pylab.yticks(yt, b, fontsize=fs-6)
a = list(pylab.axis())
a[0] = xt[0]
a[1] = xt[-1]
a[2] = yt[0]
a[3] = yt[-1]
pylab.axis(a)
pylab.savefig('../figs/%spermutations.pdf' % prefix)

pylab.figure(8, figsize=figsize)
pylab.clf()
u = np.random.random(10**6)
pylab.subplot(2,1,1)
pylab.hist(np.log(u), 100)
pylab.text(-15.5, 10**5, '$N$ = 1,000,000\n100 bins', fontsize=fs-2)
pylab.xticks(fontsize=fs)
(a, b) = pylab.yticks()
pylab.yticks(a[:-1:2])
pylab.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pylab.yticks(fontsize=fs)
pylab.rc('font', **{'size':fs})
pylab.ylabel('count', fontsize=fs)
pylab.title('Histograms of $\\log$ $u$, $u$ ~ Unif(0, 1)', fontsize=fs)
pylab.subplot(2,1,2)
pylab.hist(np.log(u), 100, log=True)
pylab.ylabel('count (log scale)', fontsize=fs)
pylab.xticks(fontsize=fs)
pylab.yticks(fontsize=fs)
(a, b) = pylab.yticks()
pylab.yticks(a[2:-1:2])
pylab.xlabel('$\\log$ $u$\n', fontsize=fs)
pylab.savefig('../figs/%suniform.pdf' % prefix)

pylab.figure(6, figsize=figsize)
pylab.clf()
pylab.subplot(3,1,1)
pylab.hist(log_p[0], 100)
pylab.ylabel('count', fontsize=fs-2)
pylab.xticks(fontsize=fs-2)
(a, b) = pylab.yticks()
pylab.yticks(a[::2], fontsize=fs-2)
pylab.text(-34.3, 5000, '$N$ = 50000\n100 bins', fontsize=fs-2)
pylab.text(-34.3, 40000, 'log(posterior) at current state\ns = %2.2f' % log_p[0].std(), fontsize=fs-2)
pylab.title('Histograms of per-datum contributions', fontsize=fs-2)
pylab.subplot(3,1,2)
pylab.hist(log_q[0], 100)
pylab.text(-34.3, 40000, "log(posterior) at proposal\ns = %2.2f" % log_q[0].std(), fontsize=fs-2)
pylab.ylabel('count', fontsize=fs-2)
pylab.xticks(fontsize=fs-2)
(a, b) = pylab.yticks()
pylab.yticks(a[::2], fontsize=fs-2)
pylab.subplot(3,1,3)
pylab.hist(log_q[0] - log_p[0], 100)
pylab.text(-5.7, 40000, 'log(MH ratio)\ns = %2.2f' % (log_p[0] - log_q[0]).std(), fontsize=fs-2)
pylab.ylabel('count', fontsize=fs-2)
pylab.xticks(fontsize=fs-2)
(a, b) = pylab.yticks()
pylab.yticks(a[::2], fontsize=fs-2)
pylab.xlabel('per-datum contribution', fontsize=fs-2)
pylab.savefig('../figs/%shistogram.pdf' % prefix)

"""
pylab.figure(9)
pylab.clf()
#p = np.cumsum(log_p[0]) / pos * pos[-1]
#q = np.cumsum(log_q[0]) / pos * pos[-1]
pylab.subplot(2,1,1)
pylab.hist(log_p[0], 100)
pylab.subplot(2,1,2)
pylab.hist(log_p[0] - log_q[0], 100)
"""