import numpy as np
import matplotlib.pyplot as plt
from rebin import rebin2 as rb2
from rebin import rebin1 as rb1
from smooth import smooth as smo

#sf,sl,bins = 200,800,400 (suited for ca 50,000)
sf,sl,bins = 0,0,1
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

#plt.plot(n_rb,rb1(Ns[idx],'lin',bins),color='black',label=r'$N_\mathrm{s}$',zorder=4)
plt.plot(NOISE,smo(Ns,N_smo,'lin'),color='black',label=r'$N_\mathrm{s}$',zorder=4)
#plt.plot(n_rb,rb1(Ng_BIFT[idx],'lin',bins),color='grey',label=r'$N_\mathrm{g,BIFT}$',marker=marker,zorder=5)
plt.plot(NOISE,smo(Ng_BIFT,N_smo,'lin'),color='grey',label=r'$N_\mathrm{g,BIFT}$',marker=marker,zorder=5)
#plt.plot(n_rb,rb1(Nm[idx],'lin',bins),color='blue',label=r'$N_\mathrm{m}$',marker=marker,zorder=6)
plt.plot(NOISE,smo(Nm,N_smo,'lin'),color='blue',label=r'$N_\mathrm{m}$',marker=marker,zorder=6)
#plt.plot(n_rb,Ng_rb,color='cyan',label=r'$N_\mathrm{g,BayesFit}$',marker=marker,zorder=7)
plt.plot(NOISE,smo(Ng,N_smo,'lin'),color='cyan',label=r'$N_\mathrm{g,BayesFit}$',marker=marker,zorder=7)
K=9
plt.plot(n_rb,K*np.ones(len(n_rb)),linestyle='--',color='black',label=r'$K$',zorder=8)

xlim,ylim = [-4,10],[0,25]
    
plt.xlabel('log(noise)')
plt.ylabel('Information content')
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend()

plt.tight_layout()
plt.savefig('FigureS1.pdf')
plt.show()
 
