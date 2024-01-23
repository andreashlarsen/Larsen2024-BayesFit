import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy import interpolate

def chi2r_pdf(df):
    x = np.linspace(chi2.ppf(0.00000001, df),chi2.ppf(0.99999999, df), 200)
    x_r = x/df # reduced x
    y = chi2.pdf(x, df)
    y = y/np.sum(y)
    
    return x_r,y

PLOT_R = 1
PLOT_DATA = 1

# factor
f = 1.0

plt.rcParams.update({'font.size': 14})

bins = [50,60,80,200]
prior = ['no']
weight = ['$\Sigma(\chi^2)$','$\Sigma(\chi^2/M)$','$\Sigma(\chi^2 N_g/M)$','SANS','SAXS']
wj = ['unity,  ','$1/M_j$,   ','$N_{\mathrm{g},j}/M_j$,','SANS,  ','SAXS,  ']
color = ['firebrick','green','black','green','red']
if PLOT_R:
    count = 0
    if PLOT_DATA:
        fig,ax = plt.subplots(3,2,figsize=(13,13.5))
        rows,cols = [0,0,1,1,2,2],[0,1,0,1,0,1]
        name = ['SAXS','SANS']
        color_data = ['red','blue']
        xtext,ytext,letter = -0.14,0.96,['B','C','D','E','F','G']
        for i in range(2):
            row,col = rows[count],cols[count]
            j=i+1
            filename = 'data/Isim%d_f%1.1f.dat' % (j,f)
            q,I,dI,y = np.genfromtxt(filename,usecols=[0,1,2,3],unpack=True)
            M = len(q)
            ax[row,col].plot(q,y,color='grey',label='theoretical curve',zorder=1)
            ax[row,col].errorbar(q,I,yerr=dI,linestyle='none',marker='.',color=color_data[i],label='%s data ($M=%d$)' % (name[i],M),zorder=0)
            ax[row,col].set_yscale('log')
            ax[row,col].set_xlabel(r'$q$ [$\mathrm{\AA}^{-1}$]')
            ax[row,col].set_ylabel(r'$I(q)$ [cm$^{-1}$]')
            ax[row,col].set_ylim(1E-7,1)
            ax[row,col].legend(frameon=False)
            ax[row,col].text(xtext,ytext,letter[count],transform=ax[row,col].transAxes)
            count +=1
    else:
        fig,ax = plt.subplots(2,2,figsize=(12,9))
        rows,cols = [0,0,1,1],[0,1,0,1]
    R_true = [10,30,50,70]
    minbin,maxbin = [0,20,30,69],[30,50,80,72]
    for k in range(len(R_true)):
        hist_max = 0
        row,col = rows[count],cols[count]
        for w in range(len(weight)):
            for i in range(len(prior)):
                R = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (i,w,i,w,f),skip_header=1,usecols=[k],unpack=True)
                print(len(R))
                R_av = np.mean(R)
                R_std = np.std(R)
                R_dev = np.sqrt(np.mean((R-R_true[k])**2))
                R_min,R_max = np.amin(R),np.amax(R)
                hist,bin_edges = np.histogram(R,bins=bins[k],range=(R_min*0.95,R_max*1.05))
                x = np.zeros(bins[k])
                for j in range(bins[k]):
                    x[j] = (bin_edges[j]+bin_edges[j+1])/2
                dx = x[4]-x[3]
                hist = hist/(np.sum(hist)*dx) # normalize by integral
                label = r'%-s $\Delta R$: %1.1f $\mathrm{\AA}$' % (wj[w],R_dev)
                if w < 3:
                    ax[row,col].plot(x,hist,color=color[w],label=label)
                else:
                    ax[row,col].fill_between(x,hist,color=color[w],alpha=0.4,label=label)
                if np.amax(hist) > hist_max:
                    hist_max = np.amax(hist)
        ax[row,col].plot([R_true[k],R_true[k]],[-1,10],color='grey')
        ax[row,col].legend(frameon=False,loc='upper right',title=r'weight ($w_j$)')
        ax[row,col].set_ylabel('probability density [a.u.]')
        ax[row,col].set_yticks([])
        if k == 0:
            ax[row,col].set_xlabel(r'$R_\mathrm{c}$ [$\mathrm{\AA}$]')
        else:
            ax[row,col].set_xlabel(r'$R_{%d}$ [$\mathrm{\AA}$]' % k)
        ax[row,col].set_xlim(minbin[k],maxbin[k])
        ax[row,col].set_ylim(0,hist_max*1.05)
        ax[row,col].text(xtext,ytext,letter[count],transform=ax[row,col].transAxes)
        count += 1
    plt.tight_layout()
    if PLOT_DATA:
        plt.savefig('Figure2.pdf')
    plt.show()
