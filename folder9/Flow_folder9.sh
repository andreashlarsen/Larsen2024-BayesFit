python=python3
RESTART=1
PLOT_DATA=0
PLOT_FIT=0
PLOT_POST=0
WRITE=1
imax=50000
Nalpha=11
M0=400
NOISE_MIN=-5
NOISE_MAX=11
DMAX_GUESS=140

# copy bayesfit and helper scripts to folder
cp ../*.py .

f1=1
for f2 in 1.0
do

echo f2: $f2, where f2 is the factor multiplied on the variance of the second dataset

for j in 11
do

PRIOR=$j # 0: no prior, 1: poor prior, 2: good prior, 3: best prior, 4: uniform prior (+-2sigma), 5: uniform prior (+-3sigma), 6: uniform prior (+-4sigma), 7: uniform prior (+-1sigma), 8: uniform prior (+-0.5 sigma)
echo PRIOR: $PRIOR where 0: no prior, 1: poor prior, 2: good prior, 3: best prior, 4: uniform prior 2sigma, 5: uniform prior 3sigma, 6: uniform prior 4sigma, 7: uniform 1sigma, 8: uniform 0.5sigma

for k in 0
do

WEIGHT=$k
echo WEIGHT: $WEIGHT, where 0: sum chi2, 1: sum chi2/M , 2: sum Ng*chi2/M , 3: chi2_0, 4: chi2_1

# Prepare
output_dir=prior${PRIOR}_weight${WEIGHT}
mkdir -p $output_dir

if [ $WRITE -eq 1 ]
then
if [ $RESTART -eq 1 ]
then
cat << EOF > $output_dir/R_prior${PRIOR}_weight${WEIGHT}_f${f2}.dat
# NOISE R1 R2 R3 R4 chi2r_K chi2r_Ng chi2r_M Ng Ng_BIFT Nm Ns Dmax M
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

# make theoretical curve
R1,R2,R3,R4 = 10,30,50,70
r12,r13,r14 = 2,3,4
s1,b1 = 0.5,1e-5
y1 = ff.coreshell4_ratio(q1,R1,R2,R3,R4,r12,r13,r14,s1,b1)

# simulate errors
k,c1,sc1,qa = 4500,0.85,100,0.2
noise = 10**$NOISE
def get_sigma(q,y,c,k,qa,sc):
    I = y/y[0]*sc
    idx = np.where(q>=qa)[0][0]
    Ia  = I[idx]
    var = noise*(I + 2*c*Ia/(1-c))/(k*q)
    sigma = np.sqrt(var)
    return sigma/sc

sigma1 = get_sigma(q1,y1,c1,k,qa,sc1)

# simulate data
Isim1 = np.random.normal(y1,sigma1)

# write data to file
with open('Isim1_f${f1}_noise$NOISE.dat','w') as f:
    for i in range(len(q1)):
        f.write('%20e %20e %20e %20e\n' % (q1[i],Isim1[i],sigma1[i],y1[i]))
        
with open('Isim1_f${f1}_noise${NOISE}_shanum.dat','w') as f:
    for i in range(len(q1)):
        f.write('%20e %20e %20e\n' % (q1[i],Isim1[i],sigma1[i]))

EOF
$python simulate.py
mkdir -p data
mv Isim*_f*.dat data

if [ $PLOT_DATA -eq 1 ]
then 
cat << EOF > plot_sim.py
import numpy as np
import matplotlib.pyplot as plt

q,I,dI,y = np.genfromtxt('data/Isim1_f${f1}_noise$NOISE.dat',unpack=True)
plt.errorbar(q,I,yerr=dI,linestyle='none',marker='.',color='red',label='SAXS')
plt.plot(q,y,color='red')

plt.yscale('log')
plt.xlabel('q')
plt.ylabel('I(q)')
plt.legend(frameon=False)
plt.show()

EOF
$python plot_sim.py
fi

####################
# BIFT
####################

cat << EOF > inputfile.dat 
data/Isim1_f${f1}_noise$NOISE.dat


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

####################
# Shanum
####################

#/home/andreas/ATSAS-3.2.1-1/bin/shanum data/Isim1_f${f1}_noise${NOISE}_shanum.dat --dmax=$DMAX_GUESS
#/home/andreas/ATSAS-3.2.1-1/bin/shanum data/Isim1_f${f1}_noise${NOISE}_shanum.dat --dmax=$DMAX_GUESS > output_shanum
#/home/andreas/ATSAS-3.0.3-1/bin/shanum data/Isim1_f${f1}_noise${NOISE}_shanum.dat $DMAX_GUESS
/home/andreas/ATSAS-3.0.3-1/bin/shanum data/Isim1_f${f1}_noise${NOISE}_shanum.dat $DMAX_GUESS > output_shanum

####################
# SELECT PRIOR
####################
if [ $PRIOR -eq 10 ]
then

cat << EOF > input
$output_dir
1
data/Isim1_f${f1}_noise$NOISE.dat
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

elif [ $PRIOR -eq 11 ]
then

cat << EOF > input
$output_dir
1
data/Isim1_f${f1}_noise$NOISE.dat
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

elif [ $PRIOR -eq 12 ]
then

cat << EOF > input
$output_dir
1
data/Isim1_f${f1}_noise$NOISE.dat
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

elif [ $PRIOR -eq 13 ]
then

cat << EOF > input 
$output_dir
1
data/Isim1_f${f1}_noise$NOISE.dat
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

else

echo you must choose PRIOR 10-13
rm -f input

fi
$python bayesfit.py < input > output

################
# READ OUTPUT
################

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

    elif 'Ng_tot =' in line:
        tmp = line.split('=')[1] 
        tmp2 = tmp.split('out of')[0]
        Ng_tot = float(tmp2)

    elif 'number of points in dataset 0 :' in line:
        tmp = line.split(':')[1]
        M_0 = int(tmp)

f.close()

f = open('parameters.dat','r')
lines = f.readlines()
for line in lines:
    if 'Number of good parameters' in line:
        tmp = line.split(':')[1]
        Ng_BIFT = float(tmp.split('+-')[0])
    elif 'Number of Shannon channels' in line:
        Ns = float(line.split(':')[1])
    elif 'Maximum diameter' in line:
        tmp = line.split(':')[1]
        Dmax = float(tmp.split('+-')[0])

f.close()

f = open('output_shanum','r')
lines = f.readlines()
for line in lines:
    if 'Nopt=' in line:
        Nm = float(line.split('=')[1])

f.close()

with open('$output_dir/R_prior${PRIOR}_weight${WEIGHT}_f${f2}.dat','a') as f:
    f.write('$NOISE %f %f %f %f %f %f %f %f %f %f %f %f %d\n' % (R1,R2,R3,R4,chi2r_K,chi2r_Ng,chi2r_M,Ng_tot,Ng_BIFT,Nm,Ns,Dmax,M_0))

EOF
if [ $WRITE -eq 1 ]
then
$python read_radii.py
fi
## clean up bift and shanumfiles
rm parameters.dat dummy.dat gx_pr.dat gs_pr.dat st_pr.dat fit.dat pr.dat data.dat scale_factor.dat rescale.dat inputfile.dat read_bift.py output_bift bift
rm output_shanum

done
done
done
done

## clean up
rm input output 
rm bayesfit.py formfactors.py function.py get_header_footer.py read_radii.py simulate.py
#rm -r data
