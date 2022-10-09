import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations

epsilon_m=np.loadtxt("./parameters/SpCas9_epsilon_modifed.txt")
forward_rates=np.loadtxt("./parameters/SpCas9_forward_rates.txt")


def calcu_cleavage_rate(mismatch_pos,epsilon_m,forward_rates):
    if mismatch_pos is None:
        return 0
    if type(mismatch_pos)!=type([]):
        mismatch_pos=list(mismatch_pos)
    if len(mismatch_pos)>4:
        return .0
    epsilon_m_cp=epsilon_m.copy()
    for pos in mismatch_pos:
        epsilon_m_cp[pos+1:23] = epsilon_m_cp[pos+1:23] - epsilon_m[pos+1+20]
    bkward_rate= []
    for i in range(21):
        bkward_rate.append(forward_rates[i]*np.exp(epsilon_m_cp[i])/np.exp(epsilon_m_cp[i+1]))
    gama=[]
    for i in range(21):
        gama.append(bkward_rate[i]/forward_rates[i+1])
    result=1/(1+np.sum(np.cumprod(gama)))
    return result


def get_mismatch_positions(seq1,seq2,PAMlength=3):
    seq1=seq1[:-PAMlength]
    seq2=seq2[:-PAMlength]
    mis=[]
    cnt=0
    for i in range(20):
        if seq1[i]!=seq2[i]:
            mis.append(i)
            cnt+=1
            if cnt>4:
                return None
    return 20-np.array(mis)