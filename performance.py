import sys
import math

def get_cm(filename, th=1.0e-3, pe=2, pr=1 ):
    cm = [[0,0],[0,0]]
    f = open(filename)
    for line in f:
        v = line.rstrip().split()
        e_value = float(v[pe])
        r = int(v[pr])
        if e_value <= th:
            p = 1
        else:
            p = 0
        cm[p][r] = cm[p][r]+1
    return cm
    

def get_q2(cm):
    n = float(cm[0][0]+cm[1][1]+cm[1][0]+cm[0][1])
    return (cm[0][0]+cm[1][1])/n

def get_mcc(cm):
    d = math.sqrt((cm[0][0]+cm[1][0])*(cm[0][0]+cm[0][1])*(cm[1][1]+cm[1][0])*(cm[1][1]+cm[0][1]))
    return (cm[0][0]*cm[1][1]-cm[0][1]*cm[1][0])/d

def get_tpr(cm):
    return float(cm[1][1])/(cm[1][0]+cm[1][1])

def get_ppv(cm):
    return float(cm[1][1])/(cm[1][1]+cm[0][1])


if __name__ == "__main__":
    filename = sys.argv[1]
    th = float(sys.argv[2])
    cm = get_cm(filename, th)
    q2 = get_q2(cm)
    mcc = get_mcc(cm)
    tpr = get_tpr(cm)
    ppv = get_ppv(cm)
    # print('TN = ', cm[0][1], 'FN = ', cm[0][1])
    # print('FP = ', cm[1][0], 'TP = ', cm[1][1])
    print('TH = ', th, 'Q2 = ', q2, 'MCC = ', mcc, 'TPR = ', tpr, 'PPV = ', ppv, 'CM=', cm)
    
    