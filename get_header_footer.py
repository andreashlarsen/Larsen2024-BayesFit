import numpy as np

def check_line(skip_line,line,CONTINUE):
    """
    check if a line is a header or footer line
    also check for zero and nan
    """
    
    tmp = line.split()
    imax = 3
    if len(tmp) < imax:
            imax = len(tmp)
    try:
        NAN = 0
        for i in range(imax):
            1/float(tmp[i]) # divide to ensure non-zero values
            if np.isnan(float(tmp[i])):
                NAN = 1
        if NAN:
            skip_line += 1
        else:
            CONTINUE = False
    except:
        skip_line+=1

    return skip_line,CONTINUE
    
def get_header_footer(file):
    """
    get number of headerlines and footerlines
    """

    f = open(file)
    try:
        lines = f.readlines()
    except:
        print('Error: cannot read lines of file. Do you have some special characters in the file? Try removing them and rerun')
        print('file: %s' % file)
    
    CONTINUE_H,CONTINUE_F = True,True
    header,footer,j = 0,0,0

    while CONTINUE_H or CONTINUE_F:
        
        # check if next line from top/bottom of file is a header/footer (or contains zero or nan)
        header,CONTINUE_H = check_line(header,lines[j],CONTINUE_H)
        footer,CONTINUE_F = check_line(footer,lines[-1-j],CONTINUE_F)
        
        # stop if there are no more lines, else continue to next line
        j += 1
        if j == len(lines):
            CONTINUE_H = False
            CONTINUE_F = False

    return header,footer