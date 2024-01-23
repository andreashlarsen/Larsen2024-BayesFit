# Larsen2024-BayesFit
scripts for paper Larsen2024-BayesFit

## manuscript
manuscript in preparation    
preprint available: https://arxiv.org/abs/2311.06408    

## overview of files
#### bayesfit.py  v2.3
bayesfit. takes in multiple datafiles and fits with an analytical model    

#### formfactors.py  function.py   get_header_footer.py    
helpfunctions, called by bayesapp.py

 #### Flow_template.sh Flow.sh
 bash script that calls bayesapp with various settings, as explained in the paper    

 ## installation 
 bayesfit.py is a python script (python3), requirements: python3 with numpy, matplotlib, scipy

 ## running
 python bayesfit < input

input (example for 2 datasets):
test               # output directory name
2                  # number of contrasts
Isim1.dat          # first dataset
Isim2.dat          # second dataset
coreshell4_ratio_2 # model
10 5               # mean and std dev of param 1 
30 10		            # mean and std dev of param 2
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
17.74              # Ng of first dataset (from BIFT, used to calculated the reduced chi-square, default: 2)
4.74               # Ng of second dataset (from BIFT, used to calculate the reduced chi-square, default: 2)
