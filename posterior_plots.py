# import necessary modules
import numpy as np
import triangle
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

def Rz(z,n):
    n0 = n[0]
    n1 = n[1]
    n2 = n[2]
    z1 = n[3]
    R = np.piecewise(z, [z<=z1,z>z1], [lambda x: n0*np.power((1+x),n1),
        lambda x: n0*np.power((1+x),n2)*np.power((1+z1),n1-n2)])
    return R

def plotRzPost(datafile, N, maxlike, zdatafile, outfile, annotate = False, title = None):
    data = np.loadtxt(datafile, usecols=(0,1,2,3))
    zdata = np.loadtxt(zdatafile)
    samples = np.random.choice(range(data.shape[0]),size=N,replace=False)
    factor = 4.0 * np.pi / 6.0 * 0.8 * dV / (1.0 + z)
    fig, ax = plt.subplots(2,1, figsize=(6,8))
    fig.subplots_adjust(hspace=0.05)
    hist_n, hist_bin, hist_patch = ax[1].hist(zdata, bins=Nbins, range=[0,10])
    dz = 10.0 / Nbins
    ax[1].cla()
    for i in range(N):
        n = data[samples[i],:]
        if i==0:
            ax[0].plot(z, Rz(z,n), '-b', lw=1, alpha=0.1, label='posterior sample')
        else:
            ax[0].plot(z, Rz(z,n), '-b', lw=1, alpha=0.1)
        ax[1].plot(z, Rz(z,n) * factor * df, '-b', lw=1, alpha=0.1)
    ax[0].plot(z, Rz(z,maxlike), '-k', lw=2, label='max likelihood')
    ax[1].plot(z, Rz(z,maxlike) * factor, '--k', lw=2, label=r'$N_{\mathrm{int}}$')
    ax[1].plot(z, Rz(z,maxlike) * factor * df, '-k', lw=2, label=r'$N_{\mathrm{exp}}$')
    ax[1].bar(hist_bin[:-1], hist_n/dz, width=dz, color='red', alpha=0.5, label=r'$N_{\mathrm{obs}}$')
    #ax[0].set_xlabel('Redshift')
    ax[1].set_xlabel('Redshift')
    ax[0].set_ylabel(r'$R_{\mathrm{GRB}}(z)$')
    ax[1].set_ylabel(r'$N(z)/dz$')
    ax[0].grid('on')
    ax[1].grid('on')
    ax[0].legend(loc='upper left', prop={'size':10})
    ax[1].legend(loc='upper right')
    ax[0].set_xlim([0,10])
    ax[1].set_xlim([0,10])
    ax[0].set_ylim([0,40])
    #ax[1].set_ylim([0,40])
    ax[1].set_ylim([0,(int(np.max(Rz(z,maxlike) * factor))/10 + 1)*10])
    ax[0].xaxis.set_major_locator(MultipleLocator(1))
    ax[1].xaxis.set_major_locator(MultipleLocator(1))
    ax[0].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[1].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[0].yaxis.set_major_locator(MultipleLocator(5))
    #ax[1].yaxis.set_major_locator(MultipleLocator(5))
    ax[1].yaxis.set_major_locator(MultipleLocator(20))
    ax[0].yaxis.set_minor_locator(MultipleLocator(1))
    #ax[1].yaxis.set_minor_locator(MultipleLocator(1))
    ax[1].yaxis.set_minor_locator(MultipleLocator(5))
    ax[0].set_xticklabels([])
    if annotate:
        ax[0].annotate(r'$n_0$', xy = (0.05,0.5), xytext = (0.25,7), xycoords = 'data', textcoords = 'data',
            arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        ax[0].annotate(r'$n_1$', xy = (2.5,6), xytext = (1.5,12), xycoords = 'data', textcoords = 'data',
            arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        ax[0].annotate(r'$z_1$', xy = (5.2,18), xytext = (4.25,22), xycoords = 'data', textcoords = 'data',
            arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        ax[0].annotate(r'$n_2$', xy = (6.25,7), xytext = (7.25,16), xycoords = 'data', textcoords = 'data',
            arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    if title:
        ax[0].set_title(title)
    fig.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=200)
    #plt.show()


levels = 1.0 - np.exp(-0.5 * np.linspace(1.0, 3.0, num=3) ** 2)
maxlike = {'NN': [.416,1.875,-.483,3.418,3421], 'RF': [.48,1.7,-5.934,6.857,4455],
    'AB': [.489,1.681,-5.95,6.682,4392]}
#maxpost = {'NN': [.487,1.677,-.516,8.855,5756], 'RF': [.316,2.132,-.677,5.049,5552],
#    'AB': [.756,1.224,-1.155,8.196,3974]}
for model in ['NN','RF','AB']:
    datafile = 'chains/realdata_varyz1_'+model+'_post_equal_weights.dat'
    data = np.loadtxt(datafile, usecols=(0,1,2,3,5))
    figure = triangle.corner(data, labels=[r'$n_0$',r'$n_1$',r'$n_2$', r'$z_1$', r'$N_{\rm tot}$'],
                             bins=50, quantiles=[0.05, 0.5, 0.95], show_titles=True, levels=levels,
                             title_args={"fontsize": 14}, verbose=False, smooth1d=1, smooth=1,
                             range=[(0,1.6),(0.7,3.2),(-6,0),(1.,10),(1500,10000)],
                             maxlike=maxlike[model], label_kwargs={"fontsize": 20})
    figure.savefig('../../paper/figures/realdata_varyz1_posterior_'+model+'.png')
    plt.close(figure)

bestfit = [0.42, 2.07, -0.7, 3.6, 4570.315127]
levels = 1.0 - np.exp(-0.5 * np.linspace(1.0, 3.0, num=3) ** 2)
maxlike = {'NN': [0.465,1.918,-5.961,5.617,4924], 'RF': [.408,2.074,-5.962,5.219,4937],
    'AB': [.385,2.163,-5.971,4.888,4818]}
for model in ['NN','RF','AB']:
    datafile = 'chains/bestfit_varyz1_'+model+'_post_equal_weights.dat'
    data = np.loadtxt(datafile, usecols=(0,1,2,3,5))
    figure = triangle.corner(data, labels=[r'$n_0$',r'$n_1$',r'$n_2$', r'$z_1$', r'$N_{\rm tot}$'],
                             bins=50, truths=bestfit, quantiles=[0.05, 0.5, 0.95], show_titles=True,
                             title_args={"fontsize": 14}, verbose=False, levels=levels, smooth1d=1, smooth=1,
                             range=[(0,1.6),(0.7,4.0),(-6,0),(1.,10),(2500,12000)],
                             maxlike=maxlike[model], label_kwargs={"fontsize": 20})
    figure.savefig('../../paper/figures/bestfit_varyz1_posterior_'+model+'.png')
    plt.close(figure)

"""
detfrac_data = np.loadtxt('support_data/splines_detection_fraction_z_RF.txt')
detfrac = interp1d(detfrac_data[:,0], detfrac_data[:,1], kind='linear')
Ez_data = np.loadtxt('support_data/splines_Ez.txt', usecols=(0,3))
Ez = interp1d(Ez_data[:,0], Ez_data[:,1], kind='linear')

N = 200
Nbins = 20
z = np.linspace(0,10,num=2001)
df = gaussian_filter1d(detfrac(z), 10)
dV = 75.28129176631614 * Ez(z)

plotRzPost('chains/realdata_varyz1_RF_post_equal_weights.dat', N, [.48,1.7,-5.934,6.857],
    'support_data/FynboGRB_lum_z_Zonly.txt',
    '../../paper/figures/redshift_distribution_posterior_RF_realdata.png',
    title = 'Real Data')
plotRzPost('chains/bestfit_varyz1_RF_post_equal_weights.dat', N, [.408,2.074,-5.962,5.219],
    'chains/bestfit_varyz1_RF_detectedZdata.txt',
    '../../paper/figures/redshift_distribution_posterior_RF_bestfit.png', True,
    title = 'Simulated Data')

plotRzPost('chains/realdata_varyz1_NN_post_equal_weights.dat', N, [.416,1.875,-.483,3.418],
    'support_data/FynboGRB_lum_z_Zonly.txt',
    '../../paper/figures/redshift_distribution_posterior_NN_realdata.png',
    title = 'Real Data')
plotRzPost('chains/bestfit_varyz1_NN_post_equal_weights.dat', N, [0.465,1.918,-5.961,5.617],
    'chains/bestfit_varyz1_NN_detectedZdata.txt',
    '../../paper/figures/redshift_distribution_posterior_NN_bestfit.png',
    title = 'Simulated Data')
"""
