import numpy as np
from scipy.special import jv

def psi_sphere(q,r):
    x = q*r     
    return 3*(np.sin(x)-x*np.cos(x))/x**3
    
def P_sphere(q,r):
    return psi_sphere(q,r)**2   
   
def V_sphere(r):
    return 4*r**3*np.pi/3
             
def sphere2(q,r,s):
    return s*P_sphere(q,r)  

def sphere3(q,r,s,b):
    return s*P_sphere(q,r)+b

def sin(x):
    return np.sin(x)

def cos(x):
    return np.cos(x)

def sqrt(x):
    return np.sqrt(x)

def exp(x):
    return np.exp(x)

def bessj0(x):
    """
    Bessel function of the first kind of zeroth order
    """
    return jv(0,x)

def bessj1(x):
    """"
    Bessel function of the first kind of first order
    """
    return jv(1,x)

def bessj1c(x):
    """
    bessj1(x)/x
    """
    y = np.ones(len(x))*0.5
    idx = np.where(x != 0)
    y[idx] = bessj1(x[idx])/x[idx]
    return y

def sinc(x):
    """
    function for calculating sinc = sin(x)/x
    numpy.sinc is defined as sinc(x) = sin(pi*x)/(pi*x)
    """
    return np.sinc(x/np.pi)

def psi(q,alf,r,h):
    x=q*r*sin(alf)
    y=q*h*cos(alf)/2.
    psi=2*bessj1c(x)*sinc(y)
    return psi

def V(r,h):
    return np.pi*r**2*h

def nanodisc(q,Bg,c,V_l,V_t,CV_p,Nlip,T,sigmaR,Ar,eps,n_w,Rg):
    pi = np.pi

    # Volumes (V) and scattering lengths (b)
    # specific for the DLPC/MSP1D1 nanodisc with tags
    V_p = 54293.0
    V_c = 3143.0
    V_s = 30.0

    b_p = 23473.0
    b_h = 164.0723
    b_t = 178.0440
    b_c = 1.4250e+03
    b_s = 10.0

    # constants(avogadros number,electron scattering length)
    N_A = 6.022e23
    b_e = 2.82e-13

    # derived params
    V_h = V_l - V_t
    V_p = CV_p*V_p
    V_c = CV_p*V_c

    # add 7.9 water molecules per lipid headgroup (Kucerka et al., 2005)
    b_h = b_h + n_w * b_s
    V_h = V_h + n_w * V_s

    # reparametrization (from vol to scattering contrasts)
    p_s = b_s/V_s # scattering length density of solvent
    dp_p = b_p/V_p - p_s
    dp_c = b_c/V_c - p_s
    dp_t = b_t/V_t - p_s
    dp_h = b_h/V_h - p_s

    xx=(q*Rg)**2
    P_c=(exp(-xx)+xx-1)/(xx/2.)
    P_tot=0
    F_tot=0

    b=sqrt(abs(Nlip*Ar/(2*pi*eps)))
    a=eps*b
    jmax,kmax=20,20
    dalf=pi/(2.*jmax)
    dfi=pi/(2.*kmax)

    for j in range(1,jmax+1):
        alf=j*dalf
        for k in range(1,kmax+1):
            fi=k*dfi
            r_t=sqrt((a*sin(fi))**2+(b*cos(fi))**2)
            R=sqrt( ((a+T)*sin(fi))**2 +((b+T)*cos(fi))**2)
            h_p=V_p/(pi*((a+T)*(b+T)-a*b))
            h_t=2.0*V_t/Ar
            h_h=V_h/Ar
            H=h_t+2.0*h_h

            Reff=R+abs(Rg)
            yy=q*Reff*sin(alf)
            ya=q*h_p*cos(alf)/2.

            psi_cc=(1-exp(-xx))/xx*bessj0(yy)*sin(ya)/ya

            tail=psi(q,alf,r_t,h_t)

            pro=V(R,h_p)*psi(q,alf,R,h_p)-V(r_t,h_p)*psi(q,alf,r_t,h_p)
            pro=pro/(V(R,h_p)-V(r_t,h_p))

            head=(H*psi(q,alf,r_t,H)-h_t*psi(q,alf,r_t,h_t))/(2*h_h)

            V_nd=Nlip*(V_t+V_h)+V_p
            dp_nd=(dp_t*Nlip*V_t+dp_h*Nlip*V_h+dp_p*V_p)/V_nd

            psi_nd=dp_t*(Nlip*V_t)*tail+dp_h*(Nlip*V_h)*head+dp_p*V_p*pro

            psi_nd=psi_nd/(dp_nd*V_nd)

            S_nd=psi_nd**2
            S_nd_c=psi_cc*psi_nd
            S_cc=psi_cc**2


            F=(dp_nd*V_nd)**2*S_nd+4*dp_c*V_c*dp_nd*V_nd*S_nd_c+2*(dp_c*V_c)**2*(S_cc+P_c)

            F_tot=F_tot+F*sin(alf)

    F_tot=F_tot*dalf*dfi*(2./pi)

    V_tot=V_nd+2*V_c
    dp_tot=(V_nd*dp_nd+2*V_c*dp_c)/V_tot
    P_tot=F_tot/(dp_tot*V_tot)**2
    y=c*1.e-9*N_A*F_tot*exp(-q**2*sigmaR**2)*b_e**2+Bg

    return y

def psi_coreshell2(q,r1,r2,p1,p2):
    V1,V2 = V_sphere(r1),V_sphere(r2)
    s1,s2 = psi_sphere(q,r1),psi_sphere(q,r2)
    Vs1 = V1*s1
    Vs2 = V2*s2

    A = p1*Vs1 + p2*(Vs2-Vs1)
    B = p1*V1  + p2*(V2 -V1)

    return A/B

def P_coreshell2(q,r1,r2,p1,p2):
    return psi_coreshell2(q,r1,r2,p1,p2)**2

def coreshell2(q,r1,r2,p1,p2,s,b):
    return s * P_coreshell2(q,r1,r2,p1,p2) + b

def coreshell2_2_constrain_p(q1,q2,r1,r2,p1,p2,s1,b1,s2,b2):
    y1 = coreshell2(q1,r1,r2,p1,p2,s1,b1)
    y2 = coreshell2(q2,r1,r2,p1,p2/80,s2,b2)
    return [y1,y2]

def psi_coreshell3(q,r1,r2,r3,p1,p2,p3):
    V1,V2,V3 = V_sphere(r1),V_sphere(r2),V_sphere(r3)
    s1,s2,s3 = psi_sphere(q,r1),psi_sphere(q,r2),psi_sphere(q,r3)
    Vs1 = V1*s1
    Vs2 = V2*s2
    Vs3 = V3*s3

    A = p1*Vs1 + p2*(Vs2-Vs1) + p3*(Vs3-Vs2)
    B = p1*V1  + p2*(V2 -V1)  + p3*(V3 -V2) 

    return A/B

def P_coreshell3(q,r1,r2,r3,p1,p2,p3):
    return psi_coreshell3(q,r1,r2,r3,p1,p2,p3)**2

def coreshell3(q,r1,r2,r3,p1,p2,p3,s,b):
    return s * P_coreshell3(q,r1,r2,r3,p1,p2,p3) + b

def coreshell3_2(q1,q2,r1,r2,r3,p11,p12,p13,p21,p22,p23,s1,b1,s2,b2):
    y1 = coreshell3(q1,r1,r2,r3,p11,p12,p13,s1,b1)
    y2 = coreshell3(q2,r1,r2,r3,p21,p22,p23,s2,b2)
    return [y1,y2]

def coreshell3_2_no_scale(q1,q2,r1,r2,r3,p11,p12,p13,p21,p22,p23,b1,b2):
    y1 = coreshell3(q1,r1,r2,r3,p11,p12,p13,1,b1)
    y2 = coreshell3(q2,r1,r2,r3,p21,p22,p23,1,b2)
    return [y1,y2]

def coreshell3_2_constrain_p(q1,q2,r1,r2,r3,p1,p2,p3,s1,b1,s2,b2):
    y1 = coreshell3(q1,r1,r2,r3,p1,p2,p3,s1,b1)
    y2 = coreshell3(q2,r1,r2,r3,p1*2,p2*6/10,p3*2,s2,b2)
    return [y1,y2]

def psi_coreshell4(q,r1,r2,r3,r4,p1,p2,p3,p4):
    V1,V2,V3,V4 = V_sphere(r1),V_sphere(r2),V_sphere(r3),V_sphere(r4)
    s1,s2,s3,s4 = psi_sphere(q,r1),psi_sphere(q,r2),psi_sphere(q,r3),psi_sphere(q,r4)
    Vs1 = V1*s1
    Vs2 = V2*s2
    Vs3 = V3*s3
    Vs4 = V4*s4

    A = p1*Vs1 + p2*(Vs2-Vs1) + p3*(Vs3-Vs2) + p4*(Vs4-Vs3)
    B = p1*V1  + p2*(V2-V1)   + p3*(V3-V2)   + p4*(V4-V3) 
 
    return A/B

def P_coreshell4(q,r1,r2,r3,r4,p1,p2,p3,p4):
    return psi_coreshell4(q,r1,r2,r3,r4,p1,p2,p3,p4)**2

def coreshell4(q,r1,r2,r3,r4,p1,p2,p3,p4,s,b):
    return s * P_coreshell4(q,r1,r2,r3,r4,p1,p2,p3,p4) + b

def coreshell4_2(q1,q2,r1,r2,r3,r4,p11,p12,p13,p14,p21,p22,p23,p24,s1,b1,s2,b2):
    y1 = coreshell4(q1,r1,r2,r3,r4,p11,p12,p13,p14,s1,b1)
    y2 = coreshell4(q2,r1,r2,r3,r4,p21,p22,p23,p24,s2,b2)
    return [y1,y2]

def coreshell2_ratio(q,R1,R2,r2,s,b):
    p1 = 1
    p2 = r2
    return s * P_coreshell2(q,R1,R2,p1,p2) + b

def coreshell3_ratio(q,R1,R2,R3,r2,r3,s,b):
    p1 = 1
    p2 = r2
    p3 = r3
    return s * P_coreshell3(q,R1,R2,R3,p1,p2,p3) + b

def coreshell4_ratio(q,R1,R2,R3,R4,r2,r3,r4,s,b):
    return s * P_coreshell4(q,R1,R2,R3,R4,1,r2,r3,r4) + b

def coreshell4_ratio_2(q1,q2,R1,R2,R3,R4,r12,r13,r14,r22,r23,r24,s1,b1,s2,b2):
    y1 = s1 * P_coreshell4(q1,R1,R2,R3,R4,1,r12,r13,r14) + b1
    y2 = s2 * P_coreshell4(q2,R1,R2,R3,R4,1,r22,r23,r24) + b2
    return [y1,y2]