python=python3
RESTART=0
KLD=0
PLOT_DATA=0
PLOT_FIT=0
PLOT_POST=0
WRITE=1
imax=10000
Nalpha=11
M0=400
M1=50
NOISE_MIN=0
NOISE_MAX=0
DMAX_GUESS=140

# copy bayesfit and helper scripts to folder
cp ../*.py .

f1=1.0
for f2 in 1.0
do

echo f2: $f2, where f2 is the factor multiplied on the variance of the second dataset

for j in 0
do

PRIOR=$j # 0: no prior, 1: poor prior, 2: good prior, 3: best prior, 4: unif prior (+-2sig), 5: unif prior (+-3sig), 6: unif prior (+-4sig), 7: unif prior (+-1sig), 8: unif prior (+-0.5 sig). SAXS:+10 (10-13). SANS:+20 (20-23).
echo PRIOR: $PRIOR where 0: no prior, 1: poor prior, 2: good prior, 3: best prior, 4: unif prior 2sig, 5: unif prior 3sig, 6: unif prior 4sig, 7: unif prior 1sig, 8: unif prior 0.5 sig. SAXS: plus 10 folders, 10-13. SANS: plus 20, folders 20-23.

for k in 0 1 2 3 4
do

WEIGHT=$k
echo WEIGHT: $WEIGHT, where 0: sum chi2, 1: sum chi2/M , 2: sum Ng*chi2/M , 3: chi2_0, 4: chi2_1

# Prepare
output_dir=prior${PRIOR}_weight${WEIGHT}
mkdir -p $output_dir

if [ $KLD -eq 1 ]
then
if [ $RESTART -eq 1 ]
then
cat << EOF >  $output_dir/KLD_prior$PRIOR.dat
# KLD for
# R1 R2 R3 R4 
EOF
fi
fi

if [ $WRITE -eq 1 ]
then
if [ $RESTART -eq 1 ]
then
cat << EOF > $output_dir/R_prior${PRIOR}_weight${WEIGHT}_f${f2}.dat
# R1 R2 R3 R4 chi2r_K chi2r_Ng chi2r_M chi2r_0_K chi2r_0_Ng_tot chi2r_0_Ng_0 chi2r_0_Ng_X_0 chi2r_0_M chi2r1_K chi2r_1_Ng_tot chi2r_1_Ng_1 chi2r_1_Ng_X_1 chi2r_1_M Ng_tot Ng_0 Ng_X_0 Ng_1 Ng_X_1 M_0 M_1
EOF
fi
fi

for i in $(seq 1 $imax)
do 

NOISE=$($python -c "print($NOISE_MIN+($i-0.5)*($NOISE_MAX-$NOISE_MIN)/($imax))")
if [ $NOISE_MIN -eq $NOISE_MAX ]
then
echo i = $i of $imax 
else
echo i = $i of $imax and NOISE = $NOISE
fi

cat << EOF > simulate.py
import numpy as np
import formfactors as ff

# make q
q1 = np.linspace(0.001,0.5,$M0)
q2 = np.linspace(0.001,0.5,$M1)

# make theoretical curve
R1,R2,R3,R4 = 10,30,50,70
r12,r13,r14 = 2,3,4
r22,r23,r24 = -0.1,0.1,0.05
s1,b1 = 0.5,1e-5
s2,b2 = 0.8,1e-4
y1,y2 = ff.coreshell4_ratio_2(q1,q2,R1,R2,R3,R4,r12,r13,r14,r22,r23,r24,s1,b1,s2,b2)

# simulate errors
k,c1,c2,sc1,sc2,qa = 4500,0.85,0.90,100,10,0.2
noise = 10**$NOISE
def get_sigma(q,y,c,k,qa,sc):
    I = y/y[0]*sc
    idx = np.where(q>=qa)[0][0]
    Ia  = I[idx]
    var = noise*(I + 2*c*Ia/(1-c))/(k*q)
    sigma = np.sqrt(var)
    return sigma/sc

sigma1 = get_sigma(q1,y1,c1,k,qa,sc1)
sigma2 = get_sigma(q2,y2,c2,k,qa,sc2)

# simulate data
Isim1 = np.random.normal(y1,sigma1)
Isim2 = np.random.normal(y2,sigma2)

# write data to file
with open('Isim1_f$f1.dat','w') as f:
    for i in range(len(q1)):
        f.write('%20e %20e %20e %20e\n' % (q1[i],Isim1[i],$f1*sigma1[i],y1[i]))

with open('Isim2_f$f2.dat','w') as f:
    for i in range(len(q2)):
        f.write('%20e %20e %20e %20e\n' % (q2[i],Isim2[i],$f2*sigma2[i],y2[i]))

EOF
$python simulate.py
mkdir -p data
mv Isim*_f*.dat data

if [ $PLOT_DATA -eq 1 ]
then 
cat << EOF > plot_sim.py
import numpy as np
import matplotlib.pyplot as plt

q,I,dI,y = np.genfromtxt('data/Isim1_f$f1.dat',unpack=True)
plt.errorbar(q,I,yerr=dI,linestyle='none',marker='.',color='red',label='SAXS')
plt.plot(q,y,color='red')

q,I,dI,y = np.genfromtxt('data/Isim2_f$f2.dat',unpack=True)
plt.errorbar(q,I,yerr=dI,linestyle='none',marker='.',color='green',label='SANS')
plt.plot(q,y,color='green')

plt.yscale('log')
plt.xlabel('q')
plt.ylabel('I(q)')
plt.legend(frameon=False)
plt.show()

EOF
$python plot_sim.py
fi

if [ $k -eq 2 ]
then
####################
# BIFT to find Ng
####################

cat << EOF > inputfile.dat 
data/Isim1_f${f1}.dat


5000
f$DMAX_GUESS





70

D
Y
C


EOF
cp ../results_bift/bift .
cat << EOF > read_bift.py
f = open('parameters.dat','r')
lines = f.readlines()
for line in lines:
    if 'Number of good parameters' in line:
        tmp = line.split(':')[1]
        Ng_BIFT = float(tmp.split('+-')[0])
        print('%1.2f' % Ng_BIFT)
f.close()
EOF

./bift < inputfile.dat > output_bift
NG0=$($python read_bift.py)
echo Ng SAXS from BIFT is $NG0
rm parameters.dat # to ensure that we get error if next run is not producint parameters.dat

sed -i -e "s/Isim1_f${f1}.dat/Isim2_f${f2}.dat/g" inputfile.dat
./bift < inputfile.dat > output_bift
NG1=$($python read_bift.py)
echo Ng SANS from bift is $NG1
rm parameters.dat dummy.dat gx_pr.dat gs_pr.dat st_pr.dat fit.dat pr.dat data.dat scale_factor.dat rescale.dat inputfile.dat read_bift.py output_bift bift

else
NG0=1 # not used
NG1=1 # not used 
fi

####################
# SELECT PRIOR
####################
if [ $PRIOR -eq 0 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 5
30 10
50 15
70 20
2 0.2
3 0.3
4 0.4
-0.1 0.01
0.1 0.01
0.05 0.005
0.5 0.05
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 1 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
5 5
40 10
45 15
90 20
2 0.2
3 0.3
4 0.4
-0.1 0.01
0.1 0.01
0.05 0.005
0.5 0.05
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-1.5 2 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 2 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
8 4
35 10
40 20
80 20
2 0.2
3 0.3
4 0.4
-0.1 0.01
0.1 0.01
0.05 0.005
0.5 0.05
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-1 2 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 3 ]
then

cat << EOF > input 
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 5
30 10
50 15
70 20
2 0.2
3 0.3
4 0.4
-0.1 0.01
0.1 0.01
0.05 0.005
0.5 0.05
1e-5 1e-2
0.8 0.08
1e-4 1e-2
1 6 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 4 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 2
30 4
50 6
70 8
2 0.08
3 0.12
4 0.16
-0.1 0.004
0.1 0.004
0.05 0.002
0.5 0.01
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 5 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 3
30 6
50 9
70 12
2 0.12
3 0.18
4 0.24
-0.1 0.006
0.1 0.006
0.05 0.003
0.5 0.01
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 6 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 4
30 8
50 12
70 16
2 0.16
3 0.24
4 0.32
-0.1 0.008
0.1 0.008
0.05 0.004
0.5 0.01
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 7 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 1
30 2
50 3
70 4
2 0.04
3 0.06
4 0.08
-0.1 0.002
0.1 0.002
0.05 0.001
0.5 0.01
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 8 ]
then

cat << EOF > input
$output_dir
2
data/Isim1_f$f1.dat
data/Isim2_f$f2.dat
coreshell4_ratio_2
10 0.5
30 1
50 1.5
70 2
2 0.02
3 0.03
4 0.04
-0.1 0.001
0.1 0.001
0.05 0.0005
0.5 0.01
1e-5 1e-2
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
$NG1
EOF

elif [ $PRIOR -eq 10 ]
then

cat << EOF > input
$output_dir
1
data/Isim1_f$f1.dat
coreshell4_ratio
10 5
30 10
50 15
70 20
2 0.2
3 0.3
4 0.4
0.5 0.05
1e-5 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
EOF

elif [ $PRIOR -eq 20 ]
then

cat << EOF > input
$output_dir
1
data/Isim2_f$f2.dat
coreshell4_ratio
10 5
30 10
50 15
70 20
-0.1 0.01
0.1 0.01
0.05 0.005
0.8 0.08
1e-4 1e-2
-10 -10 1
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG1
EOF

elif [ $PRIOR -eq 11 ]
then

cat << EOF > input
$output_dir
1
data/Isim1_f$f1.dat
coreshell4_ratio
5 5
40 10
45 15
90 20
2 0.2
3 0.3
4 0.4
0.5 0.05
1e-5 1e-2
-1.5 2 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
EOF

elif [ $PRIOR -eq 21 ]
then

cat << EOF > input
$output_dir
1
data/Isim2_f$f2.dat
coreshell4_ratio
5 5
40 10
45 15
90 20
-0.1 0.01
0.1 0.01
0.05 0.005
0.8 0.08
1e-4 1e-2
-1.5 2 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG1
EOF

elif [ $PRIOR -eq 12 ]
then

cat << EOF > input
$output_dir
1
data/Isim1_f$f1.dat
coreshell4_ratio
8 4
35 10
40 20
80 20
2 0.2
3 0.3
4 0.4
0.5 0.05
1e-5 1e-2
-1.5 2.5 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
EOF

elif [ $PRIOR -eq 22 ]
then

cat << EOF > input
$output_dir
1
data/Isim2_f$f2.dat
coreshell4_ratio
8 4
35 10
40 20
80 20
-0.1 0.01
0.1 0.01
0.05 0.005
0.8 0.08
1e-4 1e-2
-1.5 2.5 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG1
EOF


elif [ $PRIOR -eq 13 ]
then

cat << EOF > input 
$output_dir
1
data/Isim1_f$f1.dat
coreshell4_ratio
10 5
30 10
50 15
70 20
2 0.2
3 0.3
4 0.4
0.5 0.05
1e-5 1e-2
-1 6 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG0
EOF

elif [ $PRIOR -eq 23 ]
then

cat << EOF > input 
$output_dir
1
data/Isim2_f$f2.dat
coreshell4_ratio
10 5
30 10
50 15
70 20
-0.1 0.01
0.1 0.01
0.05 0.005
0.8 0.08
1e-4 1e-2
-1 6 $Nalpha
$PLOT_FIT
$PLOT_POST
$WEIGHT
$NG1
EOF

else

echo you must choose PRIOR 0-8, 10-13 or 20-23
rm -f input

fi
$python bayesfit.py < input > output

################
# READ OUTPUT
################

cat << EOF > read_kld.py
f = open('output','r')
lines = f.readlines()
for line in lines:
    if 'R1 = ' in line:
        try:
            tmp = line.split('kld =')[1]
            if $PRIOR in [0,10,20]:
                KLD_R1 = float(tmp.split(', kld_uniform = ')[1])
            else:
                KLD_R1 = float(tmp.split(',')[0])
        except:
            pass #print('skip')
    elif 'R2 = ' in line:
        try:
            tmp = line.split('kld =')[1]
            if $PRIOR in [0,10,20]:
                KLD_R2 = float(tmp.split(', kld_uniform = ')[1])
            else:
                KLD_R2 = float(tmp.split(',')[0])
        except:
            pass #print('skip')
    elif 'R3 = ' in line:
        try:
            tmp = line.split('kld =')[1]
            if $PRIOR in [0,10,20]:
                KLD_R3 = float(tmp.split(', kld_uniform = ')[1])
            else:
                KLD_R3 = float(tmp.split(',')[0])
        except:
            pass #print('skip')
    elif 'R4 = ' in line:
        try:
            tmp = line.split('kld =')[1]
            if $PRIOR in [0,10,20]:
                KLD_R4 = float(tmp.split(', kld_uniform = ')[1])
            else:
                KLD_R4 = float(tmp.split(',')[0])
        except:
            pass #print('skip')

f.close
with open('$output_dir/KLD_prior$PRIOR.dat','a') as f:
    f.write('%f %f %f %f \n' % (KLD_R1,KLD_R2,KLD_R3,KLD_R4))
EOF
if [ $KLD -eq 1 ]
then
$python read_kld.py
fi

cat << EOF > read_radii.py
f = open('output','r')
lines = f.readlines()
for line in lines:
    if 'R1 = ' in line:
        tmp = line.split('=')[1]
        R1 = float(tmp.split('+/-')[0])
    elif 'R2 = ' in line:
        tmp = line.split('=')[1]
        R2 = float(tmp.split('+/-')[0])
    elif 'R3 = ' in line:
        tmp = line.split('=')[1]
        R3 = float(tmp.split('+/-')[0])
    elif 'R4 = ' in line:
        tmp = line.split('=')[1]
        R4 = float(tmp.split('+/-')[0])

    elif 'chi2r_K =' in line:
        tmp = line.split('=')[1]
        chi2r_K = float(tmp)
    elif 'chi2r_Ng =' in line:
        tmp = line.split('=')[1]
        chi2r_Ng = float(tmp)
    elif 'chi2r_M =' in line:
        tmp = line.split('=')[1]
        chi2r_M = float(tmp)
    
    elif 'chi2r_0_K =' in line:
        tmp = line.split('=')[1]
        chi2r_0_K = float(tmp)
    elif 'chi2r_0_Ng_tot =' in line:
        tmp = line.split('=')[1]
        chi2r_0_Ng_tot = float(tmp)
    elif 'chi2r_0_Ng_0 =' in line:
        tmp = line.split('=')[1]
        chi2r_0_Ng_0 = float(tmp)
    elif 'chi2r_0_Ng_X_0 =' in line:
        tmp = line.split('=')[1]
        chi2r_0_Ng_X_0 = float(tmp)
    elif 'chi2r_0_M =' in line:
        tmp = line.split('=')[1]
        chi2r_0_M = float(tmp)

    elif 'chi2r_1_K =' in line:
        tmp = line.split('=')[1]
        chi2r_1_K = float(tmp)
    elif 'chi2r_1_Ng_tot =' in line:
        tmp = line.split('=')[1]
        chi2r_1_Ng_tot = float(tmp)
    elif 'chi2r_1_Ng_1 =' in line:
        tmp = line.split('=')[1]
        chi2r_1_Ng_1 = float(tmp)
    elif 'chi2r_1_Ng_X_1 =' in line:
        tmp = line.split('=')[1]
        chi2r_1_Ng_X_1 = float(tmp)
    elif 'chi2r_1_M =' in line:
        tmp = line.split('=')[1]
        chi2r_1_M = float(tmp)

    elif 'Ng_tot =' in line:
        tmp = line.split('=')[1] 
        tmp2 = tmp.split('out of')[0]
        Ng_tot = float(tmp2)
    elif 'Ng for dataset 0 =' in line:
        tmp = line.split('=')[1]
        tmp2 = tmp.split('out of')[0]
        Ng_0 = float(tmp2)
    elif 'Ng for dataset 1 =' in line:
        tmp = line.split('=')[1]
        tmp2 = tmp.split('out of')[0]
        Ng_1 = float(tmp2)    
    elif 'Ng_X for dataset 0 =' in line:
        tmp = line.split('=')[1]
        tmp2 = tmp.split('out of')[0]
        Ng_X_0 = float(tmp2)
    elif 'Ng_X for dataset 1 =' in line:
        tmp = line.split('=')[1]
        tmp2 = tmp.split('out of')[0]
        Ng_X_1 = float(tmp2)

    elif 'number of points in dataset 0 :' in line:
        tmp = line.split(':')[1]
        M_0 = int(tmp)
    elif 'number of points in dataset 1 :' in line:
        tmp = line.split(':')[1]
        M_1 = int(tmp)

f.close()

with open('$output_dir/R_prior${PRIOR}_weight${WEIGHT}_f${f2}.dat','a') as f:
    try:
        # if two datasets
        f.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d\n' % (R1,R2,R3,R4,chi2r_K,chi2r_Ng,chi2r_M,chi2r_0_K,chi2r_0_Ng_tot,chi2r_0_Ng_0,chi2r_0_Ng_X_0,chi2r_0_M,chi2r_1_K,chi2r_1_Ng_tot,chi2r_1_Ng_1,chi2r_1_Ng_X_1,chi2r_1_M,Ng_tot,Ng_0,Ng_X_0,Ng_1,Ng_X_1,M_0,M_1))
    except:
        # if one dataset
        f.write('%f %f %f %f %f %f %f %f \n' % (R1,R2,R3,R4,chi2r_K,chi2r_Ng,chi2r_M,Ng_tot))
EOF
if [ $WRITE -eq 1 ]
then
$python read_radii.py
fi

done
done
done
done

## clean up
rm input output 
rm bayesfit.py formfactors.py function.py get_header_footer.py read_kld.py read_radii.py simulate.py
