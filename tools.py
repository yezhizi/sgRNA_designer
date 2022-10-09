import regex
def findAllpSite(data):
    pattern1=regex.compile(r"[ATCG]{21}GG")
    pattern2=regex.compile(r"CC[ATCG]{21}")
    res1=regex.findall(pattern1,data,overlapped=True)
    res2=regex.findall(pattern2,data,overlapped=True)
    return res1,res2

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