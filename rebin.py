import sys
import numpy as np

def rebinSAS(q,I,dI,option,b):
    """                   
    rebin small-angle scattering data
    q rebinned with simple average to q_RB
    I rebinned with simple average to I_RB
    dI rebinned with square average (error propagation) to dI_RB
    option: 'log' or 'lin' for logarithmic or linear rebinning
    b: binsize (if option = 'lin'), factor to increase binsize logarithmically (if option = 'log')
    """

    q_RB = rebin(q,option,b,'simple')
    I_RB = rebin(I,option,b,'simple')
    dI_RB = rebin(dI,option,b,'square')
    
    print('REBIN:%d data points rebinned to %d data points' % (len(q),len(q_RB)))

    return q_RB,I_RB,dI_RB

def rebin3(x1,x2,x3,option,b):
    """
    rebin 3D element
    """

    x1_RB = rebin(x1,option,b,'simple')
    x2_RB = rebin(x2,option,b,'simple')
    x3_RB = rebin(x3,option,b,'simple')

    return x1_RB,x2_RB,x3_RB

def rebin2(x1,x2,option,b):
    """ 
    rebin 2D element
    """

    x1_RB = rebin(x1,option,b,'simple')
    x2_RB = rebin(x2,option,b,'simple')
   
    return x1_RB,x2_RB

def rebin1(x1,option,b):
    """ 
    rebin 1D element
    """

    x1_RB = rebin(x1,option,b,'simple')

    return x1_RB

def rebin(x_vec,option,b,average_type):
    """ 
    rebin 1D element

    input:
    x_vec : vector to rebin
    option: 
        'log' : logarithmic rebinning
        'lin' : linear rebinning
    b (lin): binsize, integer. E.g. 10
    b (log): factor to increase binsize logarithmically, float. E.g. 1.04
    average_type: 
        'simple' : simple average, sum(x_i)/N
        'square' : squared average, sqrt(sum(x_i**2))/N
 
    output:
    x_RB: rebinned vector
    """
    
    # set bin size
    if option == 'lin':
        binsize = b
    elif option == 'log':
        binsize = 1
    else:
        print('\nERROR!! rebin(): option should be \'log\' or \'lin\'\n')
        sys.exit()
    
    # initialize rebinned vector
    x_RB = []
    # reset sum and count for first bin
    sum,count = 0,0

    # rebin 
    for x in x_vec:
        if count >= binsize:
            # append bin average value to x_RB
            average = get_average(sum,count,average_type)
            x_RB.append(average)
            # reset sum and count for next bin
            sum,count = 0,0
            # increase binsize if logarithmic rebinning
            if option == 'log':
                binsize *= b
        sum,count = add_x(sum,count,average_type,x)        

    # last bin
    if average_type == 'square':
        average = np.sqrt(sum)/count
    elif average_type == 'simple':
        average = sum/count    
    
    x_RB.append(average)
    
    return np.array(x_RB)


def add_x(sum,count,average_type,x):
    """
    add x or x-square to sum and add 1 to count
    """
    if average_type == 'square':
        sum += x**2
    elif average_type == 'simple':
        sum += x
    else:
        print('\nERROR!! rebin(): average_type should be \'simple\' or \'square\'\n')
        sys.exit()
    count += 1
    return sum,count

def get_average(sum,count,average_type):
    """
    calculate simple or square average from sum and count
    """
    if average_type == 'square':
        average = np.sqrt(sum)/count
    elif average_type == 'simple':
        average = sum/count
    return average
