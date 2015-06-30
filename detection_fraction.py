import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

nn = np.loadtxt('support_data/splines_detection_fraction_z.txt')
rf = np.loadtxt('support_data/splines_detection_fraction_z_RF.txt')
ab = np.loadtxt('support_data/splines_detection_fraction_z_AB.txt')
fl = np.loadtxt('support_data/splines_detection_fraction_z_flux.txt')
z = nn[:,0]
an = np.piecewise(z, [z<=5.96,z>5.96], [lambda x: -0.01+1.02*np.exp(-x/1.68), lambda x: 0.02])

data1 = np.loadtxt('../traindata/Swift_train_all.txt', usecols=(2,26))
data2 = np.loadtxt('../traindata/Swift_validate_all.txt', usecols=(2,26))
data = np.vstack((data1,data2))
total, bin_edges_1 = np.histogram(data[:,0], density=False, range=(0,10), bins=50)
detect, bin_edges_2 = np.histogram(data[data[:,1]==1.0,0], density=False, range=(0,10), bins=50)
td = 1.0 * detect / total
bin_centers = (bin_edges_1[:-1] + bin_edges_1[1:]) / 2.0
tderr = 1.0 * np.minimum(np.sqrt(detect),np.ones(detect.shape)) / total


fig,ax = plt.subplots(1,3)

ax[0].plot(nn[:,0],nn[:,1],'-k',lw=2,label="NN")
ax[0].plot(rf[:,0],rf[:,1],'-r',lw=2,label="RF")
ax[0].plot(ab[:,0],ab[:,1],'-b',lw=2,label="AB")
ax[0].plot(fl[:,0],fl[:,1],'-g',lw=2,label="flux")
ax[0].plot(z,an,'-c',lw=2,label="analytic")
ax[0].grid()
ax[0].legend(loc='upper right')
ax[0].set_ylim(0,1)

ax[1].plot(nn[:,0],nn[:,1],'-k',lw=2,label="NN")
ax[1].plot(rf[:,0],rf[:,1],'-r',lw=2,label="RF")
ax[1].plot(ab[:,0],ab[:,1],'-b',lw=2,label="AB")
ax[1].plot(fl[:,0],fl[:,1],'-g',lw=2,label="flux")
ax[1].plot(z,an,'-c',lw=2,label="analytic")
ax[1].grid()
ax[1].legend(loc='upper right')
ax[1].set_yscale('log')
ax[1].set_ylim(1e-4,1)

ax[2].plot(nn[:,0],nn[:,1]/rf[:,1],'-k',lw=2,label="NN vs RF")
ax[2].plot(ab[:,0],ab[:,1]/rf[:,1],'-r',lw=2,label="AB vs RF")
ax[2].legend(loc='upper left')
ax[2].grid()

plt.show()


fig,ax = plt.subplots(1)
ax.errorbar(bin_centers,td,fmt=u'-m',lw=2,label="Data")
ax.plot(nn[:,0],nn[:,1],'-k',lw=2,label="NN")
ax.plot(rf[:,0],rf[:,1],'-r',lw=2,label="RF")
ax.plot(ab[:,0],ab[:,1],'-b',lw=2,label="AB")
ax.plot(fl[:,0],fl[:,1],'-g',lw=2,label="flux")
ax.plot(z,an,'-c',lw=2,label="analytic")
ax.grid("on")
ax.legend(loc='lower left')
ax.set_yscale('log')
ax.set_ylim(1e-4,1)
ax.set_xlabel(r'$z$', fontsize=18)
ax.set_ylabel(r'$F_{\rm det}(z)$', fontsize=18)
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.25))
plt.show()
fig.savefig('../../paper/figures/detection_fraction.png', bbox_inches='tight', pad_inches=0.05, dpi=200)

