import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy import interpolate

def chi2r_pdf(df):
    x = np.linspace(chi2.ppf(0.00000001, df),chi2.ppf(0.99999999, df), 200)
    x_r = x/df # reduced x
    y = chi2.pdf(x, df)
    y = y/np.sum(y)#*N
    
    return x_r,y

PLOT_ALL = 0
PLOT_X = 1
bins = 40

prior = ['unif$_{5\sigma}$','unif$_{4\sigma}$','unif$_{3\sigma}$','unif$_{2\sigma}$','unif$_{1\sigma}$','unif$_{0.5\sigma}$','poor','good','best'] 
color = ['firebrick','orange','green','magenta','brown','lightgreen','cyan','blue','black']
R_prior = [[],[],[],[],[],[],[[5,5],[40,10],[45,15],[90,20]],[[8,4],[35,10],[40,20],[80,20]],[[10,5],[30,10],[50,15],[70,20]]]
prior_number = [10,6,5,4,7,8,11,12,13]

xx = np.linspace(0,300,1000)
dxx = xx[4]-xx[3]
xx_max = [30,80,110,160]

if PLOT_ALL:
    prior_range = range(len(prior))
else:
    prior_range = [0,6,7,8]
    prior[0] = 'unif.'

M = [400,50]
M_sum = np.sum(M)
K = 14
lw = 2.0
contrast = ['SAXS','SANS']
ymax = [6,2.5]
xlim = [[0.75,1.3],[0.3,1.8]]
if PLOT_X:
    fig,ax = plt.subplots(4,2,figsize=(10,12))
    xtext2,ytext2,letter = -0.14,0.96,['A','C','E','G','B','D','F','H']
    count = 0
    rows,cols = [0,1,2,3,0,1,2,3],[0,0,0,0,1,1,1,1]
    for jj in range(2):
        xtext = xlim[jj][0]+(xlim[jj][1]-xlim[jj][0])*0.05
        for i in prior_range:
            row,col = rows[count],cols[count]
            nu = ['K','N_\mathrm{g}','0']
            color = ['blue','grey','red']
            for k in range(len(nu)):
                col_data=4+k
                nn = prior_number[i]+10*jj
                filename = 'prior%d_weight0/R_prior%d_weight0_f1.0.dat' % (nn,nn)
                chi2r,Ng = np.genfromtxt(filename,skip_header=1,usecols=[col_data,7],unpack=True)
                print(len(chi2r))
                hist,bin_edges = np.histogram(chi2r,bins=bins,range=xlim[jj])
                x = np.zeros(bins)
                for j in range(bins):
                    x[j] = (bin_edges[j]+bin_edges[j+1])/2
                dx = x[4]-x[3]
                hist = hist/(np.sum(hist)*dx)
                av = np.mean(chi2r)
                label = r'$\nu = %s\ \langle\chi^2_\mathrm{r}\rangle$ = %1.2f' % (nu[k],av)
                ax[row,col].plot(x,hist,linewidth=lw,color=color[k],label=label)
                Ng_tot = np.mean(Ng)
        
            ## theoretical distribution
            df = M[jj]-Ng_tot
            xx,y = chi2r_pdf(df)
            dxx = xx[4]-xx[3]
            xx,y = chi2r_pdf(df)
            func_int = interpolate.interp1d(xx,y,kind='cubic')
            y_bin_centers = func_int(xx)
            y_norm = y_bin_centers/(np.sum(y_bin_centers)*dxx)
            ax[row,col].plot(xx,y_norm,lw=1,color='black',linestyle='--',alpha=1,zorder=10)
            ax[row,col].set_xlim(xlim[jj])
            ax[row,col].set_ylim(0,ymax[jj])
            ax[row,col].legend(frameon=False)
            ax[row,col].set_ylabel('probability density')
            ax[row,col].text(xtext,0.9*ymax[jj],'%s prior' % prior[i])
            if row < 3:
                ax[row,col].set_xticks([])
            else:
                ax[row,col].set_xlabel(r'$\chi_r^2$ (fit to %s data)' % contrast[jj])
            ax[row,col].text(xtext2,ytext2,letter[count],transform=ax[row,col].transAxes)
            count += 1
            ax[row,col].set_yticks([])
    plt.tight_layout()
    plt.savefig('FigureS2.pdf')
    plt.show()
