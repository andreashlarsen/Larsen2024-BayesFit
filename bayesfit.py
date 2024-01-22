"""
bayesfit
version 2.2

INSTALLATION
no installation required

RUN
python bayesfit < input

input (example for 2 datasets):
test               # output directory name
2                  # number of contrasts
Isim1.dat          # first dataset
Isim2.dat          # second dataset
coreshell4_ratio_2 # model
10 5               # mean and std dev of param 1 
30 10		       # mean and std dev of param 2
50 15              # ...
70 20              # ...
2 0.2              #
3 0.3              #
4 0.4              #
-0.1 0.01          #
0.1 0.01           #
0.05 0.005         #
0.5 0.05           #
1e-4 1e-2          #
0.8 0.08           # ...
1e-4 1e-2          # mean and std dev of param K
-10 -10 1          # alpha scan: startvalue, endvalue, number of steps
0                  # plot data (0: no, 1: yes)
4                  # weight scheme, 0: sum(chi2), 1: sum(chi2/M), 2: sum(Ng*chi2/M), 3: only first dataset, 4: only second dataset
17.74              # Ng of first dataset (from BIFT)
4.74               # Ng of second dataset (from BIFT)
"""

##############################
### IMPORT PYTHON PACKAGES ###
##############################
import numpy as np
import matplotlib.pyplot as plt
import os
from io import StringIO
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.special import rel_entr
from shutil import rmtree

##############################
### IMPORT OWN FUNCTIONS #####
##############################
from function import function
from get_header_footer import get_header_footer
import formfactors as ff

##############################
### READ USER INPUT ##########
##############################

## OUTPUT NAME
output_dir = input('-- OUTPUT DIRECTORY --\n')
os.makedirs(output_dir,exist_ok=True)
print('   output dir: %s' % output_dir)

## NUMBER OF CONTRASTS
Ncontrasts_str = input('-- NUMBER OF CONTRASTS --\n')
Ncontrasts = int(Ncontrasts_str)
print('   contrasts: %s' % Ncontrasts)

## DATA
print('-- INPUT DATA  --')
M,q,I,dI,dataname = [],[],[],[],[]
for ii in range(Ncontrasts): 
    datapath  = input()
    print('   data : %s' % datapath)
    skip_header,skip_footer = get_header_footer(datapath)
    print('        number of headerlines = %d, number of footerlines = %d' % (skip_header,skip_footer))
    # import data
    q_ii,I_ii,dI_ii = np.genfromtxt(datapath.strip(),usecols=[0,1,2],skip_header=skip_header,skip_footer=skip_footer,unpack=True)
    M_ii = len(q_ii)
    print('        number of points in dataset %d : %d' % (ii,M_ii))
    M.append(M_ii)
    q.append(q_ii)
    I.append(I_ii)
    dI.append(dI_ii)
    dataname.append(datapath.split('/')[-1])

## MODEL
modelname = input('-- MODEL --\n')
print('   model: %s' % modelname)
if modelname == 'nanodisc':
    model = ff.nanodisc
    p_name = ['Bg','c','V_l','V_t','CV_p','Nlip','T','sigmaR','Ar','eps','n_w','Rg']
elif modelname == 'sphere3':
    model = ff.sphere3
    p_name = ['radius','scale','background']
elif modelname == 'coreshell4_2':
    model = ff.coreshell4_2
    p_name = ['r1','r2','r3','r4','pX1','pX2','pX3','pX4','pN1','pN2','pN3','pN4','sX','bX','sN','bN']
elif modelname == 'coreshell3_2':
    model = ff.coreshell3_2
    p_name = ['r1','r2','r3','pX1','pX2','pX3','pN1','pN2','pN3','sX','bX','sN','bN']
elif modelname == 'coreshell3_2_constrain_p':
    model = ff.coreshell3_2_constrain_p
    p_name = ['r1','r2','r3','pX1','pX2','pX3','sX','bX','sN','bN']
elif modelname == 'coreshell3':
    model = ff.coreshell3
    p_name = ['r1','r2','r3','p1','p2','p3','sX','bX']
elif modelname == 'coreshell2_2_constrain_p':
    model = ff.coreshell2_2_constrain_p
    p_name = ['r1','r2','pX1','pX2','sX','bX','sN','bN']
elif modelname == 'coreshell2':
    model = ff.coreshell2
    p_name = ['r1','r2','p1','p2','s','b']
elif modelname == 'coreshell2_ratio':
    model = ff.coreshell2_ratio
    p_name = ['R1','R2','r2','s','b']
elif modelname == 'coreshell3_ratio':
    model = ff.coreshell3_ratio 
    p_name = ['R1','R2','R3','r2','r3','s','b']
elif modelname == 'coreshell4_ratio':
    model = ff.coreshell4_ratio
    p_name = ['R1','R2','R3','R4','r2','r3','r4','s','b']
elif modelname == 'coreshell4_ratio_2':
    model = ff.coreshell4_ratio_2
    p_name = ['R1','R2','R3','R4','rX2','rX3','rX4','rN2','rN3','rN4','sX','bX','sN','bN']
else:
    print('unknown model:  %s' % modelname)
    print('please add model to bayesfit.py')
    exit()

## PRIOR
print('-- PRIOR --')
p0_list,dp0_list = [],[]
for i in range(len(p_name)):
    tmp = input()
    p0_tmp,dp0_tmp = np.genfromtxt(StringIO(tmp),unpack=True)
    p0_list.append(p0_tmp)
    dp0_list.append(dp0_tmp)
    print('   %s = %f +/- %f' % (p_name[i],p0_tmp,dp0_tmp))
p0 = np.array(p0_list)
dp0 = np.array(dp0_list)

## ALPHA SCAN
tmp = input('-- LOGALPHA SCAN --\n')
logalpha_min,logalpha_max,logalpha_n = np.genfromtxt(StringIO(tmp),unpack=True)
print('   logalpha scan from %1.2f to %1.2f (n = %1.0f):' % (logalpha_min,logalpha_max,logalpha_n))
logalpha_scan = np.linspace(logalpha_min,logalpha_max,int(logalpha_n))

## PLOT or NOT?
PLOT = int(input('-- PLOT DATA? --\n'))
if PLOT:
    print('   yes')
else:
    print('   no')
PLOT_POST = int(input('-- PLOT PRIOR and POSTERIOR? --\n'))
if PLOT_POST:
    print('   yes')
else:
    print('   no')

## Calculate Ng for selected params:
p_sel = [0,1,2,3]

## Weight in minimization
WEIGHT = int(input('-- WEIGHT MODEL --\n'))
if WEIGHT == 0:
    print('   %d: sum(chi2)' % WEIGHT)
elif WEIGHT == 1:
    print('   %d: sum(chi2/M)' % WEIGHT)
elif WEIGHT == 2:
    print('   %d: sum(Ng*chi2/M)' % WEIGHT)
    Ng_BIFT = []
    for ii in range(Ncontrasts):
        Ng_ii = input()
        Ng_BIFT.append(float(Ng_ii))
        print('   Ng of dataset %d: %s' % (ii,Ng_ii))
elif WEIGHT == 3:
    print('   only first dataset considered')
elif WEIGHT == 4:
    print('   only second dataset considered')
else:
    print('   ERROR: WEIGHT must be either 0,1,2,3 or 4')

##############################
### STOP USER INPUT ##########
##############################

# make summary file
with open('%s/summary_%s.dat' % (output_dir,output_dir),'w') as f:
    f.write('alpha chi2 chi2r S logT P Ng\n')

K = len(p0)

fit_list,ev_list,p_list,dp_list,chi2r_list,Ng_list,Ng_ii_list,Ng_X_ii_list,Ng_all_list,Ng_sel_list = [],[],[],[],[],[],[],[],[],[]
if Ncontrasts > 1:
        chi2r_ii_matrix = [[],[]]

print('-- PROGRESS LOGALPHA SCAN --')
for logalpha in logalpha_scan:
    alpha = 10**logalpha
    a = np.sqrt(alpha)
    func = function(a,K,M,model)
     
    # merge q,I with p,dp
    q_dummy = np.ones(len(p0))*99
    q_merge,I_merge,dI_merge = [],[],[]
    for ii in range(Ncontrasts):
        w = 1
        if WEIGHT == 1:
            w = np.sqrt(1/M[ii])
        elif WEIGHT == 2:
            w = np.sqrt(Ng_BIFT[ii]/M[ii])
        elif WEIGHT == 3 and ii == 0:
            w = 1e-10 # not zero do avoid division by zero
        elif WEIGHT == 4 and ii == 1:
            w = 1e-10
        q_merge = np.concatenate((q_merge,q[ii]),axis=None)
        I_merge = np.concatenate((I_merge,I[ii]),axis=None)
        dI_merge = np.concatenate((dI_merge,dI[ii]/w),axis=None)
    x = np.concatenate((q_merge,q_dummy),axis=None)
    y = np.concatenate((I_merge,a*p0),axis=None)
    dy = np.concatenate((dI_merge,dp0),axis=None)
    
    # fit
    lb,ub = p0-5*dp0,p0+5*dp0
    for k in range(4):
        if lb[k] < 0:
            lb[k] = 0
    popt,pcov = curve_fit(func,x,y,sigma=dy,absolute_sigma=True,p0=p0,bounds=(lb,ub))
    dpopt = np.sqrt(np.diag(pcov)) 
    def get_fit(q,model,p):
        if len(q) == 1:
            fit = model(q[0],*p)
            return [fit]
        elif len(q) == 2:
            fit = model(q[0],q[1],*p)
            return fit
    fit = get_fit(q,model,popt)

    # calc chi2 and chi2r
    chi2 = 0
    for ii in range(Ncontrasts):
        R = (I[ii]-fit[ii])/dI[ii]
        chi2_ii = np.sum(R**2)
        chi2r_ii = chi2_ii/(M[ii]-K)
        if Ncontrasts > 1:
            chi2r_ii_matrix[ii].append(chi2r_ii)
        chi2 += chi2_ii 
    M_merge = len(q_merge) 
    chi2r = chi2/(M_merge-K)

    # calc S
    R_S = (p0-popt)/dp0
    S = np.sum(R_S**2)

    # calc Q 
    Q = chi2 + alpha*S
 
    def fit_merge(fit):
        fit_merge = []
        for i in range(len(fit)):
            fit_merge = np.concatenate((fit_merge,fit[i]),axis=None)
        return fit_merge

    # estimate matrices B = nabla nabla chi2
    # BB is the unitless version of B, also without factor 2
    BB = np.zeros((K,K))
    dI2 = dI_merge**2.0
    BB_ii,dI2_ii = [],[]
    for ii in range(Ncontrasts):
        BB_ii.append(np.zeros((K,K)))
        dI2_ii.append(dI[ii]**2.0)
    eps = 0.001
    for i in range(K):
        di = popt[i]*eps
        popt[i] += di
        fit_plus = get_fit(q,model,popt)
        popt[i] -= di
        dIdi = (fit_merge(fit) - fit_merge(fit_plus))/di
        dIdi_ii = []
        for ii in range(Ncontrasts):
            dIdi_ii.append((fit[ii] - fit_plus[ii])/di)
        for j in range(K):
            dj = popt[j]*eps
            popt[j] += dj
            fit_plus = get_fit(q,model,popt)
            popt[j] -= dj
            dIdj = (fit_merge(fit) - fit_merge(fit_plus))/dj
            dI2dr = 2*np.sum(dIdi*dIdj/dI2)
            dp2 = dp0[i]*dp0[j]
            BB[i,j] = dI2dr*dp2/2.0
            for ii in range(Ncontrasts):
                dIdj_ii = (fit[ii] - fit_plus[ii])/dj
                dI2dr_ii = 2*np.sum(dIdi_ii[ii]*dIdj_ii/dI2_ii[ii])
                BB_ii[ii][i,j] = dI2dr_ii*dp2/2.0

    etaBB = np.linalg.eig(BB)[0]
    Ng_all = etaBB/(alpha+etaBB) 
    Ng = np.sum(Ng_all)

    Ng_ii,Ng_X_ii = [],[]
    for ii in range(Ncontrasts):
        etaBB_ii = np.linalg.eig(BB_ii[ii])[0]
        Ng_ii.append(np.sum(etaBB_ii/(alpha+etaBB_ii)))
    Ng_sum = np.sum(Ng_ii)
    for ii in range(Ncontrasts):
        wx = (M_merge-M[ii])/M_merge
        Ng_X_ii.append(Ng_ii[ii] - (Ng_sum-Ng)*wx)

    BB_sel = BB[p_sel,:][:,p_sel]
    etaBB_sel = np.linalg.eig(BB_sel)[0]
    Ng_sel = np.sum(etaBB_sel/(alpha+etaBB_sel))

    # matrix C = nabla nabla Q
    # CC is the unitless version of C, also without factor 2
    CC = BB
    for i in range(K):
        CC[i,i] += alpha
    
    # calc detC
    detCC = np.linalg.det(CC)

    # calc logT
    # A = nabla nabla alpha*S
    # AA is the unitless version of A, also without factor 2 
    detAA = alpha**K
    logT = np.log(detCC/detAA)

    # Jeffreys prior
    jef = 2*np.log(alpha)

    # evidence
    ev = Q + logT + jef
    aS = alpha*S

    # print
    if 0:
        print('                       chi2  = %e' % chi2)
        print('                       aS = %e' % aS)
        print('                       chi2r  = %f' % chi2r)
        print('                       logT  = %f' % logT)
        print('                       Q  = %e' % Q)
        print('                       ev = %e' % ev)
        print('                       Ng = %f' % Ng)
        print('                       Ng_sel = %f' % Ng_sel)    

    ev_list.append(ev)
    p_list.append(popt)
    dp_list.append(dpopt)
    chi2r_list.append(chi2r)
    Ng_list.append(Ng)
    Ng_all_list.append(Ng_all)
    Ng_ii_list.append(Ng_ii)
    Ng_X_ii_list.append(Ng_X_ii)
    Ng_sel_list.append(Ng_sel)
    fit_list.append(fit)
    
    with open('%s/summary_%s.dat' % (output_dir,output_dir),'a') as f:
        f.write('%e %f %f %f %f %f\n' % (alpha,chi2,chi2r,S,logT,Ng)) 
    
    print('   logalpha = %1.4f, evidence = %1.4f'% (logalpha,ev))

# weighted average
ev_min = np.amin(ev_list)
delta_ev = ev_list - ev_min
Pnorm = np.exp(-delta_ev/2)
sumP = np.sum(Pnorm) 

## PLOT PROBABILITY
if PLOT and len(logalpha_scan) > 1:
    print('-- PLOT PROBABILITY --')
    print('   done')
    plt.plot(logalpha_scan,Pnorm,color='black')
    plt.plot(logalpha_scan,logalpha_scan*0,linestyle='--',color='grey')
    plt.xlabel('log(alpha)')
    plt.ylabel('Probability = exp(-ev/2)')
    plt.show()

# weighted average
def av(x,w):
    d = len(np.shape(x))
    x_av = np.average(x,weights=w,axis=-d)   
    return x_av

fit_av = []
for i in range(Ncontrasts):
    fit_tmp = [sublist[i] for sublist in fit_list]
    fit_av.append(av(fit_tmp,Pnorm))

Ng_av = av(Ng_list,Pnorm)

Ng_ii_av,Ng_X_ii_av = [],[]
for ii in range(Ncontrasts):
    tmp_ii = [item[ii] for item in Ng_ii_list]
    tmp_X_ii = [item[ii] for item in Ng_X_ii_list]
    Ng_ii_av.append(av(tmp_ii,Pnorm))
    Ng_X_ii_av.append(av(tmp_X_ii,Pnorm))
    
Ng_all_av = []
for k in range(K):
    tmp_k = [item[k] for item in Ng_all_list]
    Ng_all_av.append(av(tmp_k,Pnorm))

Ng_sel_av = av(Ng_sel_list,Pnorm)
logalpha_av = av(logalpha_scan,Pnorm)
p = av(p_list,Pnorm)
dp = av(dp_list,Pnorm)
chi2r_av = av(chi2r_list,Pnorm)

Iprior = get_fit(q,model,p0)

if Ncontrasts > 1:
    chi2r_ii_av = []
for ii in range(Ncontrasts):
    if Ncontrasts > 1:
        chi2r_ii_av.append(av(chi2r_ii_matrix[ii],Pnorm))
    with open('%s/fit_%s_dataset%d.dat' % (output_dir,output_dir,ii),'w') as f:
        f.write('# q I dI Iprior Ifit\n')
        for j in range(M[ii]):
            f.write('%f %f %f %f %f\n' % (q[ii][j],I[ii][j],dI[ii][j],Iprior[ii][j],fit_av[ii][j]))

## GOODNESS OF FIT
print('-- GOODNESS OF FIT --')
chi2r_Ng = chi2r_av*(M_merge-K)/(M_merge-Ng_av)
chi2r_M  = chi2r_av*(M_merge-K)/M_merge
print('   chi2r_K = %1.6f' % chi2r_av)
print('   chi2r_Ng = %1.6f' % chi2r_Ng)     
print('   chi2r_M = %1.6f' % chi2r_M) 
if Ncontrasts > 1:
    for ii in range(Ncontrasts):
        chi2 = chi2r_ii_av[ii]*(M[ii]-K)
        chi2r_ii_av_K = chi2/(M[ii]-K)
        chi2r_ii_av_Ng_tot = chi2/(M[ii]-Ng_av)
        chi2r_ii_av_M = chi2/M[ii]
        chi2r_ii_av_Ng_ii = chi2/(M[ii]-Ng_ii_av[ii])
        chi2r_ii_av_Ng_X_ii = chi2/(M[ii]-Ng_X_ii_av[ii])
        print('   chi2r_%d_K = %1.6f' % (ii,chi2r_ii_av_K))
        print('   chi2r_%d_Ng_tot = %1.6f' % (ii,chi2r_ii_av_Ng_tot))
        print('   chi2r_%d_Ng_%d = %1.6f' % (ii,ii,chi2r_ii_av_Ng_ii))
        print('   chi2r_%d_Ng_X_%d = %1.6f' % (ii,ii,chi2r_ii_av_Ng_X_ii))
        print('   chi2r_%d_M = %1.6f' % (ii,chi2r_ii_av_M)) 

## INFORMATION CONTENT
print('-- INFORMATION GAIN --')
print('   logalpha of minimum = %1.6f' % logalpha_av)

# Number of good parameters
print('   Ng_tot = %1.6f out of %d' % (Ng_av,K))
print('   Ng_sel = %1.6f out of %d' % (Ng_sel_av,len(p_sel)))
for ii in range(Ncontrasts):
    print('   Ng for dataset %d = %1.6f out of %d' % (ii,Ng_ii_av[ii],K))
    print('   Ng_X for dataset %d = %1.6f out of %d' % (ii,Ng_X_ii_av[ii],K))

# Kullbeck-Leibler divergence

def kl_divergence(posterior,prior):
    """
    calculate KLD by nummerical integration
    """
    idx = np.where(posterior>0)
    if len(idx[0]) == 0:
        KLD = 0
    else:
        KLD = np.sum(posterior[idx]*np.log2(posterior[idx]/prior[idx]))/np.sum(posterior[idx])
    return KLD
    

def kl_divergence0(mu_post,sig_post,mu_prior,sig_prior):
    """
    calculate KLD by analytical expression, valid for normal (prior and posterior) distributions 
    """
    var_post = sig_post**2
    var_prior = sig_prior**2
    R = mu_prior-mu_post
    var_ratio = var_post/var_prior
    KLD = (R**2/var_prior + var_ratio-np.log(var_ratio)-1)/2
    KLD /= np.log(2) # convert from nat to bit  
    return KLD

## POSTERIOR
print('-- POSTERIOR --')

for i in range(len(p_name)):

    kld = kl_divergence0(p[i],dp[i],p0[i],dp0[i])
    x = np.linspace(p0[i]-5*dp0[i],p0[i]+5*dp0[i],10000)
    posterior_pdf = norm.pdf(x,p[i],dp[i])
    prior_pdf = norm.pdf(x,p0[i],dp0[i])
    prior_pdf_uniform = np.ones(len(x))/len(x)
    kld = kl_divergence(posterior_pdf,prior_pdf)
    kld_uniform = kl_divergence(posterior_pdf,prior_pdf_uniform)
    print('%s = %f +/- %f, kld = %1.1f, kld_uniform = %1.1f' % (p_name[i],p[i],dp[i],kld,kld_uniform))

    # plot prior and posterior distributions
    if PLOT_POST:         
        plt.plot(x,prior_pdf/np.amax(prior_pdf),linestyle='--',color='grey',label='prior')
        plt.plot(x,posterior_pdf/np.amax(posterior_pdf),color='black',label='posterior, kld = %1.1f' % kld)
        plt.title(p_name[i])
        plt.legend(frameon=False)
        plt.show()  

## PLOT DATA
if PLOT:
    print('-- PLOT DATA AND FIT --')
    print('   done')    
    color = ['red','blue','green','orange']
    offset = [1,1e2,1e4,1e6]
    for i in range(Ncontrasts):
        plt.errorbar(q[i],offset[i]*I[i],yerr=offset[i]*dI[i],linestyle='none',marker='.',color=color[i],label=dataname[i],zorder=100+i)
        if i == 0:
            plt.plot(q[i],offset[i]*Iprior[i],linestyle='--',color='grey',label='prior',zorder=300+i)
            plt.plot(q[i],offset[i]*fit_av[i],color='black',label='fit',zorder=200+i)
        else: 
            plt.plot(q[i],offset[i]*Iprior[i],linestyle='--',color='grey',zorder=300+i)
            plt.plot(q[i],offset[i]*fit_av[i],color='black',zorder=200+i)

    plt.yscale('log')
    plt.xlabel(r'$q$')
    plt.ylabel(r'$I(q)$')
    plt.legend(frameon=False)
    plt.show()

## clean up
rmtree('__pycache__')
