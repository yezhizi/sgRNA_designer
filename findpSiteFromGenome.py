from gettext import find
import regex
import os
import tools


def findpSiteFromGenome(file_dir,save_dir="./genome_AllpSite"):
    file_list=os.listdir(file_dir)
    for file_path in file_list:
        if file_path.startswith("chr"):
            path=os.path.join(file_dir,file_path)
            with open(path,"r") as f:
                f.readline()
                chr=f.read().replace("\n","").replace(" ","")
                f.close()
            seq1,seq2=tools.findAllpSite(chr)
            seq1.extend(tools.seq35_to_seq53(seq2))
            if not os.path.isdir(save_dir):
                os.mkdir(save_dir)
            save_name=os.path.join(save_dir,file_path.split(".")[0]+"_AllpSite.txt")
            with open(save_name,"w") as f:
                f.write("\n".join(seq1))
                f.close()
                
        
if __name__=="__main__":
    findpSiteFromGenome("./GCA_000002595.3")