import numpy as np
import matplotlib.pyplot as plt
    
PLOT_R = 1
PLOT_DATA = 1

plt.rcParams.update({'font.size': 14})

bins = 50
prior = ['no']
weight = ['$\Sigma(\chi^2)$','$\Sigma(\chi^2/M)$','$\Sigma(\chi^2 N_g/M)$','SANS','SAXS']
color = ['pink','red','brown','black','darkblue','royalblue','skyblue']
factors = [0.1,0.2,0.5,1,2,5,10]
xtext,ytext,letter = -0.14,0.96,['A','B','C','D','E']
if PLOT_R:
    count = 0
    if PLOT_DATA:
        fig,ax = plt.subplot_mosaic([['A','A'],['B','C'],['D','E']],figsize=(12,14))
        scale = [1E6,1E5,1E4,1E3,1E2,1E1,1]
        f = [0.1,0.2,0.5,1.0,2.0,5.0,10.0]
        color_data = ['pink','red','brown','black','darkblue','royalblue','skyblue']
        description = ['underestimated errors','underestimated errors','underestimated errors','correct errors','overestimated errors','overestimated errors','overestimated errors']
        for i in [1,3,5]:
            j=i+1
            filename = 'data/Isim2_f%1.1f.dat' % f[i]
            q,I,dI,y = np.genfromtxt(filename,unpack=True)
            s=scale[i]
            M = len(q)
            if i == 50:
                ax[letter[count]].plot(q,s*y,color='grey',label='theoretical curve',zorder=100)
            else:
                ax[letter[count]].plot(q,s*y,color='grey',zorder=100)
            ax[letter[count]].errorbar(q,s*I,yerr=s*dI,linestyle='none',marker='.',color=color_data[i],label='%1.1f (%s)' % (f[i],description[i]), zorder=100-i)
        ax[letter[count]].set_yscale('log')
        ax[letter[count]].set_xlabel(r'$q$ [$\mathrm{\AA}^{-1}$]')
        ax[letter[count]].set_ylabel(r'$I(q)$ [cm$^{-1}$]')
        ax[letter[count]].set_ylim(1E-5,1E7)
        ax[letter[count]].legend(frameon=False,title='factor multiplied on error')
        ax[letter[count]].text(xtext/2,ytext,letter[count],transform=ax['A'].transAxes)
        count += 1
    else:
        fig,ax = plt.subplots(2,2,figsize=(12,9))
    R_true = [10,30,50,70]
    minbin,maxbin = [0,20,35,68.5],[30,48,70,72.5]
    rows,cols = [0,0,1,1],[0,1,0,1]
    for k in range(len(R_true)):
        hist_max = 0
        if k == 3:
            bins = 150
        w = 0
        for jf in range(len(factors)):
            f = factors[jf]
            for i in range(len(prior)):
                R = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (i,w,i,w,f),skip_header=1,usecols=[k],unpack=True)
                print(len(R))
                R_av = np.mean(R)
                R_std = np.std(R)
                R_dev = np.sqrt(np.mean((R-R_true[k])**2))
                R_min,R_max = np.amin(R),np.amax(R)
                hist,bin_edges = np.histogram(R,bins=bins,range=(R_min*0.95,R_max*1.05))
                x = np.zeros(bins)
                for j in range(bins):
                    x[j] = (bin_edges[j]+bin_edges[j+1])/2
                dx = x[4]-x[3]
                hist = hist/(np.sum(hist)*dx) # normalize by integral
                label = r'%1.1f, $\Delta R$: %1.1f $\mathrm{\AA}$' % (f,R_dev)
                if w < 3:
                    ax[letter[count]].plot(x,hist,color=color[jf],label=label)
                else:
                    ax[letter[count]].fill_between(x,hist,color=color[jf],alpha=0.4,label=label)
                if np.amax(hist) > hist_max:
                    hist_max = np.amax(hist)
        ax[letter[count]].plot([R_true[k],R_true[k]],[-1,10],color='grey')
        ax[letter[count]].legend(frameon=False,loc='upper right',title=r'factor multiplied on error')
        ax[letter[count]].set_ylabel('probability density [a.u.]')
        ax[letter[count]].set_yticks([])
        if k == 0:
            ax[letter[count]].set_xlabel(r'$R_\mathrm{c}$ [$\mathrm{\AA}$]')
        else:
            ax[letter[count]].set_xlabel(r'$R_{%d}$ [$\mathrm{\AA}$]' % k)
        ax[letter[count]].set_xlim(minbin[k],maxbin[k])
        ax[letter[count]].set_ylim(-hist_max*0.05,hist_max*1.05)
        ax[letter[count]].text(xtext,ytext,letter[count],transform=ax[letter[count]].transAxes)
        count += 1
    plt.tight_layout()
    plt.savefig('Figure5.pdf')
    plt.show()

print(len(R))
