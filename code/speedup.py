import numpy as np
import scipy.misc as sm
import tabular as tb
import pylab

fs=20
lw=3
figsize = (8 * 1.2, 6 * 1.2)

def weak_bounds():

	pylab.ion()
	pylab.figure(1, figsize=figsize)
	pylab.clf()

	qvec = [0.001, 0.01, 0.1]

	m = np.logspace(0, 5)
	pylab.loglog(m, m, 'k:', linewidth=lw)
	
	cvec = ['b', 'c', 'g']

	"""
	for (q, c) in zip(qvec, cvec):
		K = np.log(q)/np.log(1-q)
		M = np.linspace(1, K)
		MM = np.logspace(np.log10(M[-1]), 5)
		SS = np.zeros(len(MM)) + K
		SS[MM < K] = MM[MM < K]
		SS[MM > 2**K] += np.log2(MM[MM > 2**K] - K)
		pylab.loglog(MM, SS, c + '-', linewidth=1)	
	"""

	for (q, c) in zip(qvec, cvec):
		K = np.log(q)/np.log(1-q)
		print q, K
		M = np.linspace(1, K)
		S = (1 - q - (1-q)**M) / q
		MM = np.array(M.tolist() + np.logspace(np.log10(M[-1]), 5).tolist())
		SS = np.array(S.tolist() + [S[-1]] * 50)	
		SL = SS * 0 + K
		SL[MM < K] = MM[MM < K]
		SL[MM > 2**K] += np.log2(MM[MM > 2**K] - K)
		pylab.fill_between(MM, SL, SS, color=c, alpha=0.2)

	A = 1
	for (q, c) in zip(qvec, cvec):
		p = 1 - q
		
		R = np.log(q)/np.log(1-q)

		B = p**np.arange(1, R)
		D = np.ones(len(B))
		
		B = np.concatenate([B, q * p**np.arange(0, R)])
		D = np.concatenate([D, np.arange(R)])

		for pow in range(2,88):
			B = np.concatenate([B, q**pow * p**np.arange(0, R)])
			D = np.concatenate([D, sm.comb(np.arange(R), pow)])
			assert len(B) == len(D), 'len(B) != len(D)'
			if len(B) > 10**8:
				print pow, 'breaking'
				break
			
		
		B = np.concatenate(([1], B))
		D = np.concatenate(([1], D))		
		i = B.argsort()[::-1]
		B = (D[i] * B[i]).cumsum()
		D = D[i].cumsum()
		j = np.nonzero((D >= A) & (D <= 10**5))[0]		
		#pylab.loglog(np.arange(A, 100001), A*C[np.arange(A-1, 100000)/A])
		pylab.loglog(D[j], A*B[j/A], c, linewidth=lw)
		pylab.draw()

	pylab.loglog(m, np.log2(m+1), 'purple', linewidth=lw)
		
	pylab.yscale('log')
	pylab.xscale('log')
		
		#pylab.loglog(MM, SS, c + '-', linewidth=1)	

	pylab.xlabel('Number of cores $(J)$', fontsize=fs)
	pylab.ylabel('Expected speedup $(E[S_J])$', fontsize=fs)
	pylab.title('Expected speedup with simple bounds', fontsize=fs)
	pylab.legend(['$E[S_J] = J$'] + 
				 [('$q = %1.4f' % q).strip('0') + '$' for q in qvec] +
				 ['$E[S_J] = \log_2 (J+1)$'], loc='upper left', fontsize=fs)
	pylab.xticks(fontsize=fs)
	pylab.yticks(fontsize=fs)
	pylab.axis((1, 10**4, 1, 10**4))
	pylab.savefig('../figs/expected-speedup.pdf')

def upper_bound():

	pylab.ion()
	pylab.figure(3, figsize=figsize)
	pylab.clf()
	
	q = np.logspace(-3, np.log10(0.5))
	K = np.log(q)/np.log(1-q)	
	pylab.loglog(q, K, linewidth=lw)
	pylab.xlabel('q', fontsize=fs)
	pylab.ylabel('K', fontsize=fs)
	pylab.xticks(fontsize=fs)
	pylab.yticks(fontsize=fs)
	#pylab.title('Weak upper bound on speedup', fontsize=fs)
	pylab.savefig('../figs/speedup-K-q.pdf')

def multi_way(A=64):

	pylab.ion()
	pylab.figure(2, figsize=figsize)
	pylab.clf()

	qvec = [0.001, 0.01, 0.1]

	cvec = ['b', 'c', 'g']

	m = np.logspace(0, 5)
	pylab.loglog(m, m, 'k:', linewidth=lw)

	for (q, c) in zip(qvec, cvec):
		print q
		p = 1 - q
		
		R = np.floor(np.log(q)/np.log(1-q))

		B = p**np.arange(1, R)
		D = np.ones(len(B))
		
		B = np.concatenate([B, q * p**np.arange(0, R)])
		D = np.concatenate([D, np.arange(R)])
		

		for pow in range(2,88):
			B = np.concatenate([B, q**pow * p**np.arange(0, R)])
			D = np.concatenate([D, sm.comb(np.arange(R), pow)])
			assert len(B) == len(D), 'len(B) != len(D)'
			if len(B) > 10**8:
				print pow, 'breaking'
				break
			
		
		B = np.concatenate(([1], B))
		D = np.concatenate(([1], D))
		i = B.argsort()[::-1]
		B = (D[i] * B[i]).cumsum()
		D = D[i].cumsum()
		j = np.nonzero((D <= 10**5))[0]
		#pylab.loglog(np.arange(A, 100001), A*C[np.arange(A-1, 100000)/A])
		pylab.loglog(np.concatenate(([1], A*D[j])), np.concatenate(([1], A*B[j])), c, linewidth=lw)
		pylab.draw()

	pylab.loglog(np.concatenate(([1], m*A)), np.concatenate(([1], np.log2(m+1)*A)), 'purple', linewidth=lw)

	pylab.xlabel('Number of cores', fontsize=fs)
	pylab.ylabel('Expected speedup', fontsize=fs)
	pylab.title('Expected speedup with %d-way parallelism' % A, fontsize=fs)
	pylab.legend(['$E[S_J] = J$'] + 
				 [('$q = %1.4f' % q).strip('0') + '$' for q in qvec] +
				 ['$q = 0.5$'], loc='upper left', fontsize=fs)
	pylab.xticks(fontsize=fs)
	pylab.yticks(fontsize=fs)
	pylab.axis((1, 10**4, 1, 10**4))	
	pylab.savefig('../figs/speedup-%d.pdf' % A)
	