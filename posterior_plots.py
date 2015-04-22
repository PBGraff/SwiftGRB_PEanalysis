# import necessary modules
import numpy as np
import triangle

# read in data and set true values
#datafile1 = 'chains/analysis_n0100_n1200_n2100_seed6054_post_equal_weights.dat'
#datafile2 = 'chains/analysis_n0100_n1200_n2100_seed6055_post_equal_weights.dat'
#datafile3 = 'chains/analysis_n0100_n1200_n2100_seed6056_post_equal_weights.dat'
#n0_true = 1.00
#n1_true = 2.00
#n2_true = -1.00
#extents=[(0.15,1.05), (1.9,3.5), (-1.5,-0.25), (16,22)]

datafile1 = 'chains/analysis_n084_n1207_n270_seed9641_phys_live.points'
datafile2 = 'chains/analysis_n084_n1207_n270_seed9642_phys_live.points'
datafile3 = 'chains/analysis_n084_n1207_n270_seed9643_phys_live.points'
n0_true = 0.84
n1_true = 2.07
n2_true = -0.70

data1 = np.loadtxt(datafile1, usecols=(0,1,2,3))
data2 = np.loadtxt(datafile2, usecols=(0,1,2,3))
data3 = np.loadtxt(datafile3, usecols=(0,1,2,3))
nstar_true = n0_true * (1.0 + 3.6)**n1_true
print data1.shape, data2.shape, data3.shape

truths = [n0_true, n1_true, n2_true, nstar_true]

extents = []
for i in range(4):
	tmp1 = [data1[:,i].min(), data2[:,i].min(), data3[:,i].min(), truths[i]]
	tmp2 = [data1[:,i].max(), data2[:,i].max(), data3[:,i].max(), truths[i]]
	extents.append((min(tmp1)*0.95, max(tmp2)*1.05))
print extents

# plot it
try:
	figure = triangle.corner(data1, labels=[r'$n_0$',r'$n_1$',r'$n_2$',r'$n_{\ast}$'],
                         truths=truths, quantiles=[0.16, 0.5, 0.84], show_titles=True,
                         title_args={"fontsize": 12}, bins=30,
                         extents=extents)
	figure.savefig("test_posterior_default_n0.png")
except:
	pass

try:
	figure = triangle.corner(data2, labels=[r'$n_0$',r'$n_1$',r'$n_2$',r'$n_{\ast}$'],
                         truths=truths, quantiles=[0.16, 0.5, 0.84], show_titles=True,
                         title_args={"fontsize": 12}, bins=30,
                         extents=extents)
	figure.savefig("test_posterior_default_nstar.png")
except:
	pass

try:
	figure = triangle.corner(data3, labels=[r'$n_0$',r'$n_1$',r'$n_2$',r'$n_{\ast}$'],
                         truths=truths, quantiles=[0.16, 0.5, 0.84], show_titles=True,
                         title_args={"fontsize": 12}, bins=30,
                         extents=extents)
	figure.savefig("test_posterior_default_ntotal.png")
except:
	pass
