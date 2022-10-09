from asyncore import write
from csv import writer
from msilib import sequence
import multiprocessing as mp
from multiprocessing import pool
import numpy as np
import matplotlib.pyplot as plt
import os
import regex
from functools import partial
import pandas as pd

import datetime

def findAllpSite(data):
    pattern1=regex.compile(r"[ATCG]{21}GG")
    pattern2=regex.compile(r"CC[ATCG]{21}")
    res1=regex.findall(pattern1,data,overlapped=True)
    res2=regex.findall(pattern2,data,overlapped=True)
    return res1,res2
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

def calc_fenergey(epsilon,guide_length=20):
    new_epsilon = epsilon.copy()
    energies = -1*new_epsilon[0:(guide_length+1)] # convention: epsC>0 means downward slope
    energies[0] = new_epsilon[0]                 # convention: epsPAM>0 means upward slope
    energy_on_target=np.sum(energies)
    punishments = new_epsilon[(guide_length+1):]
    return energy_on_target,punishments
def get_energies_state20(energy_on_target,punishments,mismatch_positions, guide_length=20,rel_concentration=1.):
    if type(mismatch_positions)==type([]):
        mismatch_positions = np.array(mismatch_positions)
    punishment = 0
    if len(mismatch_positions)>0:
         punishment = np.sum(punishments[(mismatch_positions.astype(int)-1)])
    energy = energy_on_target-np.log(rel_concentration)+punishment
    return energy
def seq35_to_seq53(seq):
    remap={ord("A"):"T",
           ord("T"):"A",
           ord("C"):"G",
           ord("G"):"C"
           }
    if type(seq)==type("str"):
        return seq.translate(remap)[::-1]
    elif type(seq)==type([]):
        res=[]
        for sg in seq:
            res.append(sg.translate(remap)[::-1])
        return res

def find_sgrna_pSite(target):
    if len(target)<23:
        return []
    sgrna_53,sgrna_35=findAllpSite(target)
    sgrna_35=seq35_to_seq53(sgrna_35)
    sgrna_psite=list(set(sgrna_53)|set(sgrna_35))
    return sgrna_psite

def find_pSite(target,genome_dir,genome_pSite_dir):
    pSite_53,pSite_35=[],[]
    file_list=os.listdir(genome_dir)
    flag=0
    for file in file_list:
        if file.startswith("chr"):
            with open(os.path.join(data_dir,file)) as f:
                f.readline()
                g=f.read().replace("\n","").replace(" ","")
                f.close()
            if regex.search(target,g) is not None:
                flag=1
                g_list=g.split(target)
                for i in g_list:
                    psite1,psite2=findAllpSite(i)
                    pSite_53.extend(psite1)
                    pSite_35.extend(psite2)
            elif regex.search(seq35_to_seq53(target),g) is not None:
                flag=1
                g_list=g.split(seq35_to_seq53(target))
                for i in g_list:
                    psite1,psite2=findAllpSite(i)
                    pSite_53.extend(psite1)
                    pSite_35.extend(psite2)
            else:
                chr_name=os.path.join(genome_pSite_dir,file.split(".")[0]+"_AllpSite.txt")
                with open(chr_name,"r") as f:
                    lines = f.readlines()
                    pSite_53.extend([line.strip() for line in lines])
                    f.close()
                    del lines
            del g
    if flag==0:
        return None
        
    else:
        pSite_35=seq35_to_seq53(pSite_35)
        pSite_53.extend(pSite_35)
        return pSite_53
    
rel_concentration = 1.   
info=""

epsilon = np.loadtxt('parameters/SpCas9_epsilon.txt')
bases=np.array(["A","T","C","G"])
energy_on_target,punishments=calc_fenergey(epsilon)
score_on_target=np.exp(-energy_on_target)
index_max=np.argsort(punishments)[-4:]+1
energy_max=get_energies_state20(energy_on_target,punishments,index_max, guide_length=20,rel_concentration=1)
score_max=-np.log(np.exp(-energy_max)/(np.exp(-energy_max)+score_on_target))

def calcu_scores(pSite,sgrna):
    
    res=0
    pos=[0,0,0,0,0]
    for site in pSite:
        mismatch_positions=get_mismatch_positions(sgrna,site)
        if mismatch_positions is not None:
                pos[len(mismatch_positions)]+=1
                energy=get_energies_state20(energy_on_target,punishments,mismatch_positions, guide_length=20,rel_concentration=1.)
                res+=np.exp(-energy)
    if res==0:
        score=score_max
    else:
        score=-np.log(res/(res+score_on_target))
    score=score/score_max*100
    return score,pos


if __name__=="__main__":
    st=datetime.datetime.now()
    
    #writer=pd.ExcelWriter("-result-"+st.strftime("%Y-%m-%d-%H-%M")+".xlsx")
    
    data_dir="./GCA_000002595.3\\ncbi_dataset\\data\\GCA_000002595.3"
    genome_pSite_dir="./genome_AllpSite"
    
    gene_sequence_dir="./sequence"
    gene_dir_list=os.listdir(gene_sequence_dir)
    for gene in gene_dir_list:
        print(gene)
        writer=pd.ExcelWriter(gene+"-result-"+st.strftime("%Y-%m-%d-%H-%M")+".xlsx")
        gene_path=os.path.join(gene_sequence_dir,gene)
        for exon in os.listdir(gene_path):
            if not exon.startswith("Exon"):
                continue
            exon_path=os.path.join(gene_path,exon)
            with open(exon_path,"r") as f:
                target=f.read().replace("./n","")
                f.close()
            
            sgrna_pSite=find_sgrna_pSite(target)
            if len(sgrna_pSite)==0:
                df=pd.DataFrame({"sgrna":[],"score":[],"mismatch":[]})
                df.to_excel(writer,sheet_name=exon.split(".")[0],index=False)
                
                continue
            
            pSite=find_pSite(target,data_dir,genome_pSite_dir)
            if pSite is None:
                print("the input sequence is not found in the genome")
                df=pd.DataFrame({"sgrna":[],"score":[],"mismatch":[]})
                df.to_excel(writer,sheet_name=exon.split(".")[0],index=False)
                
                continue
            else:
                print("the input sequence is found in the genome")
                
                # 多进程
                pool=mp.Pool(processes=4)
                para=partial(calcu_scores,pSite)
                res=pool.map(para,sgrna_pSite)
                res=list(zip(*res))
                scores=res[0]
                mismatchPosCnt=res[1]
                    
                    
                ed=datetime.datetime.now()
                print("===================")
                print("耗时:{}".format(ed-st))
                print("===================") 

                df=pd.DataFrame({"sgrna":sgrna_pSite,"score":scores,"mismatch":mismatchPosCnt})
                df.sort_values(by="score",ascending=False,inplace=True)
                df.to_excel(writer,sheet_name=exon.split(".")[0],index=False)
                
                pool.close()
        writer.close()