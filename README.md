# Larsen2024-BayesFit
scripts for paper Larsen2024-BayesFit

## manuscript
manuscript in preparation    
preprint available: https://arxiv.org/abs/2311.06408    

## overview of files
#### bayesapp script: bayesfit.py
bayesfit version 2.3. The program fits multiple SAXS or SANS data simultaneously, with analytical models. Priors can be applied to each parameter in the model. 

#### supporting scripts: formfactors.py  function.py   get_header_footer.py    
helpfunctions, imported and called by bayesapp.py    

 #### bash scripts: folderX/Flow_folderX.sh
 bash script that calls bayesapp with various settings (number of points in each dataset, different priors etc)     
 overview of folders can be found in the LOOGBOOK txt file    
 these scripts can be generated by running Flow_master.sh (which uses Flow_template.sh)    
 
 #### plotting scripts: folderX/FigureX.py
 python scripts for generating the figures of the paper    

 ## installation 
 bayesfit.py is a python script (python3), requirements: python3 with numpy, matplotlib, scipy

 ## running bayesfit
 python bayesfit < inputfile

#### inputfile (example for 2 datasets):
```
output_dir         # output directory name    
2                  # number of contrasts    
Isim1.dat          # first dataset    
Isim2.dat          # second dataset    
coreshell4_ratio_2 # model    
10 5               # mean and std dev of param 1     
30 10              # mean and std dev of param 2    
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
17.74              # Ng of first dataset (from BIFT, used to calculated the reduced chi-square)    
4.74               # Ng of second dataset (from BIFT, used to calculate the reduced chi-square)
```

The number of lines in the inputfile depends on the number of datasets (number of contrasts) and the number of parameters of the model in use.     
In this example, there are two datasets/contrasts, and the model is coreshell4_ratio_2 with 14 paramters (K=14).     
Available models can be found in the bayesfit.py script, which calls the formfactor.py script.
The last two lines of the inputfile are only used for weight scheme 2. Can be left blank if this scheme is not used    

