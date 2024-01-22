import numpy as np

def calc_y(x,a,M,model,p):
    K = len(p)
    if len(M) == 1:
        q = x[:-K]
        y1 = model(q,*p)
        yS = a*p  
        y = np.concatenate((y1,yS),axis=None)
    elif len(M) == 2: 
        q1 = x[:-(K+M[1])]
        q2 = x[M[0]:-K]
        y1,y2 = model(q1,q2,*p)
        yS = a*p  
        y = np.concatenate((y1,y2,yS),axis=None)
    return y 

def function(a,K,M,model):
    if K == 1:
        def func(x,p1):
            p = np.array([p1])
            return calc_y(x,a,M,model,p)
    elif K == 2:
        def func(x,p1,p2):
            p = np.array([p1,p2])
            return calc_y(x,a,M,model,p)
    elif K == 3:
        def func(x,p1,p2,p3):
            p = np.array([p1,p2,p3])
            return calc_y(x,a,M,model,p)
    elif K == 4:
        def func(x,p1,p2,p3,p4):
            p = np.array([p1,p2,p3,p4])
            return calc_y(x,a,M,model,p)
    elif K == 5:
        def func(x,p1,p2,p3,p4,p5):
            p = np.array([p1,p2,p3,p4,p5])
            return calc_y(x,a,M,model,p)
    elif K == 6:
        def func(x,p1,p2,p3,p4,p5,p6):
            p = np.array([p1,p2,p3,p4,p5,p6])
            return calc_y(x,a,M,model,p)
    elif K == 7:
        def func(x,p1,p2,p3,p4,p5,p6,p7):
            p = np.array([p1,p2,p3,p4,p5,p6,p7])
            return calc_y(x,a,M,model,p)
    elif K == 8:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8])
            return calc_y(x,a,M,model,p)
    elif K == 9:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9])
            return calc_y(x,a,M,model,p)  
    elif K == 10:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10])
            return calc_y(x,a,M,model,p)
    elif K == 11:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11])
            return calc_y(x,a,M,model,p)
    elif K == 12:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12])
            return calc_y(x,a,M,model,p)
    elif K == 13:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13])
            return calc_y(x,a,M,model,p)    
    elif K == 14:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14])
            return calc_y(x,a,M,model,p)
    elif K == 15:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15])
            return calc_y(x,a,M,model,p)
    elif K == 16:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16])
            return calc_y(x,a,M,model,p)
    elif K == 17:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17])
            return calc_y(x,a,M,model,p)
    elif K == 18:
        def func(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18):
            p = np.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18])
            return calc_y(x,a,M,model,p)
    
    return func
