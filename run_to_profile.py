import cProfile, pstats
from io import StringIO
from utilities import *
a = Strand(['PConst','attB-1-GG','GFP','attB-1-AA','attB-2-GG','attP-1-AA',\
            'PConst','Term','attP-2-GG','attP-1-GG',\
            'GENOME'],[1,1,-1,1,1,-1,-1,1,1,-1,0])
time = np.linspace(0, 18*3, 1001)
trackGenes = ['GFP', 'RFP', 'BFP']
def mkseq(seq, allele, incu):
    #helper function to make the {int:lambda: inducer time} hash
    #inputs are the sequence of inducers e.g. ['2','1'], 
    #list of all inducers ['1','2','3'], and time used per inducer 
    out = {}; c = 0;
    lambdas = [lambda t:int(t>=0 and t<incu), lambda t:int(t>=incu and t<2*incu),\
               lambda t:int(t>=incu*2 and t<incu*3)]
    for i in seq:
        out.update({i:lambdas[c]})
        c += 1
    for i in allele:
        if i not in out:
            out.update({i:lambda t:0})
    return out
pr = cProfile.Profile()
pr.enable()
os, og = serial_wrapper(a, mkseq(list("12"), ["1","2"],18), 30, {}, trackGenes, time)
pr.disable()
s = StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())