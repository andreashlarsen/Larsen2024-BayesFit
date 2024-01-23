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

Figure = int(input('Figure number (6,7,8,13 for S3, or 14 for S4):')) # 6 or 7 or 8 or 13 (S3) or 14 (S4)
PLOT_R,NORMALIZE,PLOT_PRIOR,PLOT_ALL,PLOT_X,PLOT_X0 = 0,0,0,0,0,0
if Figure == 6:
    PLOT_R = 1
    NORMALIZE = 1
    Figure_name = 'Figure6.pdf'
elif Figure == 7:
    PLOT_X = 1
    Figure_name = 'Figure7.pdf'
elif Figure == 8:
    PLOT_X0 = 1
    Figure_name = 'Figure8.pdf'
elif Figure == 13:
    PLOT_R = 1
    NORMALIZE = 1
    PLOT_PRIOR = 1
    Figure_name = 'FigureS3.pdf'
elif Figure == 14:
    PLOT_R = 1
    NORMALIZE = 1    
    PLOT_ALL = 1
    Figure_name = 'FigureS4.pdf'
else:
    print('Figure must be 6, 7, 8, 13 (S3) or 14 (S4)')
    exit()

# weight
w = 0 

# factor
f = 1.0

if PLOT_X:
    plt.rcParams.update({'font.size': 11})
elif PLOT_X0:
    plt.rcParams.update({'font.size': 12})
else:
    plt.rcParams.update({'font.size': 14})

if PLOT_R:
    bins = [70,70,70,200]
elif Figure == 14:
    bins = [30,30,30,30]
else:
    bins = [40,40,40,40,40,40,40,40]

prior = ['unif$_{5\sigma}$','unif$_{4\sigma}$','unif$_{3\sigma}$','unif$_{2\sigma}$','unif$_{1\sigma}$','unif$_{0.5\sigma}$','poor','good','best'] 
color = ['firebrick','orange','green','magenta','brown','lightgreen','cyan','blue','black']
R_prior = [[],[],[],[],[],[],[[5,5],[40,10],[45,15],[90,20]],[[8,4],[35,10],[40,20],[80,20]],[[10,5],[30,10],[50,15],[70,20]]]
prior_number = [0,6,5,4,7,8,1,2,3]

xx = np.linspace(0,300,1000)
dxx = xx[4]-xx[3]
xx_max = [30,100,120,210]

if PLOT_ALL:
    prior_range = range(len(prior))
elif PLOT_PRIOR:
    prior_range = [6,7,8]
else:
    prior_range = [0,6,7,8]
    prior[0] = 'unif.'

if PLOT_PRIOR:
    RANGE = 1
else:
    RANGE = 0

if PLOT_R:
    xtext,ytext,letter = -0.14,0.96,['A','B','C','D']
    fig,ax = plt.subplots(2,2,figsize=(12,9))
    R_true = [10,30,50,70]
    minbin,maxbin = [4,23,45,69.4],[20,40,58,71]
    rows,cols = [0,0,1,1],[0,1,0,1]
    count = 0
    for k in range(len(R_true)):
        hist_max = 0
        row,col = rows[count],cols[count]
        hist_max = 0
        for i in prior_range:
            R = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (prior_number[i],w,prior_number[i],w,f),skip_header=1,usecols=[k],unpack=True)
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
            label = r'%s, $\Delta R$: %1.1f $\mathrm{\AA}$' % (prior[i],R_dev)
            if np.amax(hist) > hist_max:
                hist_max = np.amax(hist)
            if NORMALIZE:
                ax[row,col].plot(x,hist/np.amax(hist),color=color[i],label=label)
                hist_max = 1
            else: 
                ax[row,col].plot(x,hist,color=color[i],label=label)
            if PLOT_PRIOR:
                mu,sigma = R_prior[i][k]
                y = np.exp(-((xx-mu)/sigma)**2/2)/(sigma*np.sqrt(2*np.pi))
                if NORMALIZE:
                    ax[row,col].fill_between(xx,y/np.amax(y),color=color[i],alpha=0.3)
                else:
                    ax[row,col].fill_between(xx,y,color=color[i],alpha=0.3)
            
        ax[row,col].plot([R_true[k],R_true[k]],[-1,10],color='grey')
        if PLOT_PRIOR == 0:
            ax[row,col].plot([0,400],[0,0],color='black')
        if NORMALIZE:
            ax[row,col].set_ylabel('normalized probability density')
        else:
            ax[row,col].set_ylabel('probability density [a.u.]')
            ax[row,col].set_yticks([])
        if k == 0:
            ax[row,col].set_xlabel(r'$R_\mathrm{c}$ [$\mathrm{\AA}$]')
        else:
            ax[row,col].set_xlabel(r'$R_{%d}$ [$\mathrm{\AA}$]' % k)
        if RANGE:
            ax[row,col].set_xlim(0,xx_max[k])
        else:
            ax[row,col].set_xlim(minbin[k],maxbin[k])
        ax[row,col].set_ylim(0,hist_max*1.05)
        if PLOT_ALL:
            ax[row,col].legend(frameon=False,loc='upper right',fontsize=11.3)
        else:
            ax[row,col].legend(frameon=False,loc='upper right')
        ax[row,col].text(xtext,ytext,letter[count],transform=ax[row,col].transAxes)
        count += 1
    plt.tight_layout()

    #if PLOT_ALL:
    #    plt.savefig('R_folder2_all.pdf')
    #elif PLOT_PRIOR:
    #    plt.savefig('R_folder2_prior.pdf')
    #else:
    #    plt.savefig('R_folder2.pdf')
    plt.savefig(Figure_name)
    plt.show()
    print(len(R))

M = [400,50]
M_sum = np.sum(M)
K = 14
lw = 2.0

if PLOT_X:
    xtext,ytext,letter = -0.14,0.96,['A','B','C','D']
    fig,ax = plt.subplots(4,1,figsize=(5,12))
    count = 0
    for i in prior_range:
        nu = ['K','N_\mathrm{g}','0']
        color = ['blue','grey','red']
        for k in range(len(nu)):
            col=4+k
            chi2r,Ng = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (prior_number[i],w,prior_number[i],w,f),skip_header=1,usecols=[col,17],unpack=True)
            hist,bin_edges = np.histogram(chi2r,bins=bins[k],range=(0.7,1.3))
            x = np.zeros(bins[k])
            for j in range(bins[k]):
                x[j] = (bin_edges[j]+bin_edges[j+1])/2
            dx = x[4]-x[3]
            hist = hist/(np.sum(hist)*dx)
            av = np.mean(chi2r)
            label = r'$\nu = %s, \langle\chi^2_\mathrm{r}\rangle$ = %1.2f' % (nu[k],av)
            ax[count].plot(x,hist,linewidth=lw,color=color[k],label=label)
            Ng_tot = np.mean(Ng)
        
        ## theoretical distribution
        df = M_sum-Ng_tot
        xx,y = chi2r_pdf(df)
        dxx = xx[4]-xx[3]
        xx,y = chi2r_pdf(df)
        func_int = interpolate.interp1d(xx,y,kind='cubic')
        y_bin_centers = func_int(xx)
        y_norm = y_bin_centers/(np.sum(y_bin_centers)*dxx)
        ax[count].plot(xx,y_norm,lw=1,color='black',linestyle='--',zorder=10)
        ax[count].set_xlim(0.75,1.35)
        ax[count].set_ylim(0,6.5)
        ax[count].text(0.75+0.05*(1.35-0.75),5.5,'%s prior' % prior[i])
        ax[count].legend(frameon=False,loc='upper right',fontsize=10.5)
        if count == 3:
            ax[count].set_xlabel(r'$\chi_r^2$ (fit to SAXS and SANS data)')
        else:
            ax[count].set_xticks([])
        ax[count].set_ylabel('probability density [a.u.]')
        ax[count].set_yticks([])
        ax[count].text(xtext,ytext,letter[count],transform=ax[count].transAxes)
        count +=1
    plt.tight_layout()
    plt.savefig(Figure_name)
    plt.show()
    print(len(chi2r))

if PLOT_X0:
    xtext,ytext,letter = -0.14,0.96,['A','C','E','G','B','D','F','H']
    Ndatasets = 2
    skip = [7,12]
    minmax=[[0.75,1.6],[0.4,3.0]]
    yminmax = [[0,6],[0,2.5]]
    dataset_name = ['SAXS','SANS']
    color = ['blue','cyan','lime','dummy','red','grey']
    fig,ax = plt.subplots(4,2,figsize=(10,12))
    count = 0
    rows,cols = [0,1,2,3,0,1,2,3],[0,0,0,0,1,1,1,1]
    for ii in range(Ndatasets):
        nu = ['K','N_\mathrm{g}','n_\mathrm{%s}' % dataset_name[ii],'dummy(M-weighted Ng)','0','N_\mathrm{g,%s}' % dataset_name[ii]]
        for i in prior_range:
            row,col = rows[count],cols[count]
            Ng_tot,Ng_0,Ng_1 = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (prior_number[i],w,prior_number[i],w,f),skip_header=1,usecols=[17,18,19],unpack=True) 
            Ng_sum = Ng_0+Ng_1
            Ng_ii = [Ng_0,Ng_1]
            Ng_ii_rel = (Ng_sum-Ng_ii)/Ng_sum
            df_tot = M_sum-Ng_tot 
            df_ii = M[ii]-Ng_ii[ii]
            df_ii_rel = (Ng_sum-Ng_ii[ii])/(Ng_sum)
            for k in [0,1,2,4]:
                col_data=skip[ii]+k
                if k == 4:
                    chi2r = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (prior_number[i],w,prior_number[i],w,f),skip_header=1,usecols=[col_data],unpack=True)
                    X2 = chi2r*M[ii]
                    wx = df_ii_rel
                    Ng_X = Ng_ii[ii] - (Ng_sum-Ng_tot)*wx
                    chi2r = X2/(M[ii]-Ng_X)
                    hist,bin_edges = np.histogram(chi2r,bins=bins[k],range=minmax[ii])
                    x = np.zeros(bins[k])
                    for j in range(bins[k]): 
                        x[j] = (bin_edges[j]+bin_edges[j+1])/2
                    dx = x[4]-x[3]
                    hist = hist/(np.sum(hist)*dx)
                    av = np.mean(chi2r) 
                    label = r'$\nu = %s, \langle\chi^2_\mathrm{r}\rangle$ = %1.2f' % (nu[k+1],av)
                    ax[row,col].plot(x,hist,linewidth=lw,color=color[k+1],label=label)
                chi2r = np.genfromtxt('prior%d_weight%d/R_prior%d_weight%d_f%1.1f.dat' % (prior_number[i],w,prior_number[i],w,f),skip_header=1,usecols=[col_data],unpack=True)
                hist,bin_edges = np.histogram(chi2r,bins=bins[k],range=minmax[ii])
                x = np.zeros(bins[k])
                for j in range(bins[k]):
                    x[j] = (bin_edges[j]+bin_edges[j+1])/2
                dx = x[4]-x[3]
                hist = hist/(np.sum(hist)*dx)
                av = np.mean(chi2r)
                label = r'$\nu = %s, \langle\chi^2_\mathrm{r}\rangle$ = %1.2f' % (nu[k],av)
                ax[row,col].plot(x,hist,linewidth=lw,color=color[k],label=label)

            df = M[ii]-np.mean(Ng_X)
            xx,y = chi2r_pdf(df)
            dxx = xx[4]-xx[3]
            func_int = interpolate.interp1d(xx,y,kind='cubic')
            y_bin_centers = func_int(xx)
            y_norm = y_bin_centers/(np.sum(y)*dxx)
            ax[row,col].plot(xx,y_norm,lw=1,color='black',linestyle='--',zorder=10)
            ax[row,col].set_xlim(minmax[ii])
            ax[row,col].set_ylim(yminmax[ii])
            ax[row,col].legend(frameon=False,loc='upper right',fontsize=11)
            ax[row,col].text(minmax[ii][0]+0.05*(minmax[ii][1]-minmax[ii][0]),yminmax[ii][1]*0.9,'%s prior' % prior[i])
            if row < 3:
                ax[row,col].set_xticks([])
            else:
                ax[row,col].set_xlabel(r'$\chi_r^2$ (for %s, fit against SAXS and SANS data)' % dataset_name[ii])
            ax[row,col].set_yticks([])
            ax[row,col].set_ylabel('probability density [a.u.]')
            ax[row,col].text(xtext,ytext,letter[count],transform=ax[row,col].transAxes)
            count += 1
    plt.tight_layout()
    plt.savefig(Figure_name)
    plt.show()
    print(len(chi2r))
