import numpy as np
import os
import triangle
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib.ticker import MultipleLocator

def PlotTriangle(fileroot,usetruths=True):
	data = np.loadtxt(fileroot+'post_equal_weights.dat', usecols=(0,1,2,3,4,5))
	if (usetruths):
		truths = np.loadtxt(fileroot+'injected_values.txt')
		figure = triangle.corner(data, labels=names, truths=truths, bins=30, quantiles=[0.05, 0.5, 0.95],
								 show_titles=True, title_args={"fontsize": 12})
	else:
		figure = triangle.corner(data, labels=names, bins=30, quantiles=[0.05, 0.5, 0.95],
								 show_titles=True, title_args={"fontsize": 12})
	figure.savefig(fileroot+'posterior_plot.png')

def FindRecoveryCDF(fileroot):
	CDFvals = np.zeros(6)
	data = np.loadtxt(fileroot+'post_equal_weights.dat', usecols=(0,1,2,3,4,5))
	truths = np.loadtxt(fileroot+'injected_values.txt')
	for i in range(6):
		CDFvals[i] = np.sum(data[:,i]<truths[i]) * 1.0 / data.shape[0]
	return CDFvals

def PP_KStest(cdf):
	D, p = scipy.stats.kstest(cdf, 'uniform')
	return p

def PP_Plot(idx,cdf):
	x = np.linspace(0.0,1.0,num=opts.number+1)
	y = np.zeros(opts.number+1)
	y[1:] = np.sort(cdf)
	fig, ax = plt.subplots(1)
	ax.plot(x,x,'--k',lw=2)
	ax.plot(y,x,'-b',lw=2)
	ax.axis([0,1,0,1])
	ax.set_xlabel(r'$p$')
	ax.set_ylabel(r'$\Pr(p)$')
	ax.set_title(names[idx]+r' KS p-value = %.4lf'%(PP_KStest(cdf)))
	ax.xaxis.set_major_locator(MultipleLocator(0.1))
	ax.yaxis.set_major_locator(MultipleLocator(0.1))
	ax.grid()
	#plt.show()
	fig.savefig(opts.outdir+'/pp_plot_'+names2[idx]+'.png')

def RunAnalysis(i):
	root = opts.outdir+'/run'+str(i)+'_'
	cmd = './Analysis --seed='+str(seeds[i+1])+' --nlive='+str(opts.nlive)+'  --varyz1 --silent --outfile='+root
	cmd += ' --n0='+str(n0_inj[i])+' --n1='+str(n1_inj[i])+' --n2='+str(n2_inj[i])+' --z1='+str(z1_inj[i])
	if (opts.resume):
		cmd += ' --resume'
	print cmd

	# run the analysis
	os.system(cmd)
	
	# plot posterior results
	PlotTriangle(root)

	# find recovered CDF for injected values
	recoveredCDFvals[i,:] = FindRecoveryCDF(root)

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-n","--number",action="store",type="int",default=100,help="Number of analyses to run")
parser.add_option("-o","--outdir",action="store",type="string",default="chains",help="Directory for analyses' output")
parser.add_option("-l","--nlive",action="store",type="int",default=1000,help="Number of live points to use")
parser.add_option("-r","--resume",action="store_true",default=False,help="Resume previous run")
(opts,args)=parser.parse_args()

seeds = np.random.random_integers(10,30000,opts.number+1)
recoveredCDFvals = np.zeros((opts.number,6))

names = [r'$n_0$', r'$n_1$', r'$n_2$', r'$z_1$', r'$n_{\ast}$', r'$n_{\rm tot}$']
names2 = ['n0', 'n1', 'n2', 'z1', 'nstar', 'ntot']

if not os.path.exists(opts.outdir):
	os.makedirs(opts.outdir)

# run prior sampling
cmd = './Analysis --seed='+str(seeds[0])+' --nlive='+str(opts.number)+' --outfile='+opts.outdir+'/prior_ --varyz1 --zeroLogLike --silent'
if (opts.resume):
	cmd = cmd + ' --resume'
print cmd
os.system(cmd)

# retrieve prior samples for injecting
prior_samples = np.loadtxt(opts.outdir+'/prior_post_equal_weights.dat')
n0_inj = prior_samples[:,0]
n1_inj = prior_samples[:,1]
n2_inj = prior_samples[:,2]
z1_inj = prior_samples[:,3]

# loop to perform all analyses
for i in range(opts.number):
	RunAnalysis(i)

# make P-P plots
for i in range(6):
	PP_Plot(i,recoveredCDFvals[:,i])
