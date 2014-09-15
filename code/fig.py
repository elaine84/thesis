import numpy as np
import pylab

sigma = 0.5
mu = 0.0
x = np.linspace(-3, 3, 100)
y = np.exp(-(x - mu)**2 / 2 / sigma**2) / np.sqrt(2 * np.pi) / sigma + 0.01
z = x * 0
z[20] = y[20]

pylab.ion()
pylab.figure(10, figsize=(12, 3))
pylab.clf()
fig = pylab.figure(10, frameon=False)
pylab.axis([-2, 2, 0, 1.2])
a = pylab.axes()
a.set_frame_on(False)
a.get_xaxis().tick_bottom()
pylab.fill_between(x[60:], 0, y[60:], linewidth=4, color='CornflowerBlue')
pylab.plot(x, y, linewidth=4, color='k')
pylab.plot([x[60], x[60]], [0, y.max()], '--', linewidth=6, color='orange')
pylab.plot([0, 0], [0, 1], 'k:', linewidth=2)
pylab.xticks([],[])
pylab.yticks([],[])
#pylab.xticks(fontsize=16)
fs=20
pylab.text(0.4, 0.85, '$\\log$ $\\gamma_\\rho$', fontsize=fs)
pylab.text(-0.05, 1, '$\\hat{\\mu}_m$', fontsize=fs)
pylab.text(-2.2, 0.35, "$\\cal{L}(\\theta$'$) - \\cal{L}(\\theta)$ $\sim$ $\\cal{N}($$\\hat{\\mu}_m,$ $\\hat{\\sigma}^2_m)$", fontsize=fs)
pylab.savefig('../figs/accept.pdf')
