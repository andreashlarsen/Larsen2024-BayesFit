import numpy as np
import matplotlib.pyplot as plt
from rebin import rebin2 as rb2
from rebin import rebin1 as rb1
from smooth import smooth as smo

fig,ax = plt.subplots(2,1,figsize=(4.5,8))
noise=['-2.0','-1.008','0.0','0.9935','2.08','2.992']
lognoise = np.linspace(-2,3,6)
color = ['darkred','orange','red','blue','coral','cyan']
scale = [1E15,1E12,1E9,1E6,1E3,1E0]
b = 0.0
for i in [0,2,4]:
    filename = 'data/Isim1_f1_noise%s.dat' % noise[i]
    q,I,dI,y = np.genfromtxt(filename,unpack=True)
    s=scale[i]
    ax[0].errorbar(q,s*(I+b),yerr=s*dI,linestyle='none',marker='.',color=color[i],label='%1.1f' % lognoise[i])
ax[0].set_yscale('log')
ax[0].set_xlabel(r'$q$ [$\mathrm{\AA}^{-1}$]')
ax[0].set_ylabel(r'$I(q)$ [cm$^{-1}$]')
ax[0].set_xlim(0,0.5)
ax[0].set_ylim(1E-6,1E18)
ax[0].legend(frameon=False,title='log(noise)')
ax[0].text(-0.09,1e17,'A',fontsize=14)

#sf,sl,bins = 200,800,400 (suited for ca 50,000)
sf,sl,bins = 0,0,100
N_smo = 2000

marker = ' '
for i in [11]:
    filename = 'prior%d_weight0/R_prior%d_weight0_f1.0.dat' % (i,i)
    NOISE,R1,R2,R3,R4,chi2r_K,chi2r_Ng,chi2r_M,Ng,Ng_BIFT,Nm,Ns,Dmax,M = np.genfromtxt(filename,skip_header=1+sf,skip_footer=sl,unpack=True)
    N=len(NOISE)
    print(N)
    sm = int(N/2)
    idx = np.where(Ns<23)
    n=NOISE[0+sf:-sl]
    n_rb,Ng_rb = rb2(NOISE[idx],Ng[idx],'lin',bins)
#ax[1].plot(n_rb,rb1(Ns[idx],'lin',bins),color='black',label=r'$N_\mathrm{s}$',zorder=4)
#ax[1].plot(n_rb,rb1(Ng_BIFT[idx],'lin',bins),color='grey',label=r'$N_\mathrm{g,BIFT}$',marker=marker,zorder=5)
#ax[1].plot(n_rb,rb1(Nm[idx],'lin',bins),color='blue',label=r'$N_\mathrm{m}$',marker=marker,zorder=6)
#ax[1].plot(n_rb,Ng_rb,color='cyan',label=r'$N_\mathrm{g,BayesFit}$',marker=marker,zorder=7)

ax[1].plot(NOISE,smo(Ns,N_smo,'lin'),color='black',label=r'$N_\mathrm{s}$',zorder=4)
ax[1].plot(NOISE,smo(Ng_BIFT,N_smo,'lin'),color='grey',label=r'$N_\mathrm{g,BIFT}$',marker=marker,zorder=5)
ax[1].plot(NOISE,smo(Nm,N_smo,'lin'),color='blue',label=r'$N_\mathrm{m}$',marker=marker,zorder=6)
ax[1].plot(NOISE,smo(Ng,N_smo,'lin'),color='cyan',label=r'$N_\mathrm{g,BayesFit}$',marker=marker,zorder=7)

ax[1].text(-3,24,'B',fontsize=14)

xlim,ylim = [-2.3,2.3],[0,25]
postname = ''
for x,c,zo in zip([-2,0,2],['darkred','red','coral'],[1,2,3]):
    plt.plot([x,x],ylim,linestyle='--',color=c,zorder=zo)
    
ax[1].set_xlabel('log(noise)')
ax[1].set_ylabel('Information content')
ax[1].set_xlim(xlim)
ax[1].set_ylim(ylim)
ax[1].legend()

plt.tight_layout()
plt.savefig('Figure1.pdf')
plt.show()
 
