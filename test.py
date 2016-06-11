import multiprocessing as mp
import time
import os

def f(x):
    return x.work()

def run():
	p = mp.Pool(processes=4)
	foo1 = Foo(); foo1.p = lambda x: x + 5
	foo2 = Foo(); foo2.p = lambda x: x + 10
	#r3 = p.map(f, [lambda x: x+5])
	#print(list(r3))
	#r3 = [p.apply(f, args = (i,)) for i in [lambda x: x+5]]
	r3 = [p.apply_async(f, args = (i,)) for i in [foo1, foo2]]
	ro3 = [j.get() for j in r3]
	for a in ro3:
		print(a)
import numpy as np
from utilities import *
if __name__ == '__main__':
	#seeded RNG
	np.random.seed(seed=123456)
	N = 5000
	#set up the initial condition of the circuit
	victoria = Strand(['RFP','attB-1','TERM','attB-2','attP-1','PConst','attP-2', 'GFP','GENOME'], [-1,1,-1,-1,-1,-1,-1,1,0])
	#set up time
	t = np.linspace(0, 40, 20)
	#figure out which strands should I be tracking giving a sequence of inputs
	trackStrands = sequences(victoria, {'sa':['1'], 'sab':['1','2'], 'sb':['2']})

	#set up inducer times
	inducer = {'1':[(0, 100)], '2':[(0,100)]}
	#run the simulation
	outstrand, outgene = parallel_wrapper(victoria, inducer, N, trackStrands, [], t, proc = 10)
	print(outstrand)
	#run()
def main():
	andreGenome = Strand(['attB-1'],[1])
	andreFuel = Strand(['attP-1','attB-1'],[1,1])
	victoria = Strand(['RFP','attB-1','TERM','attB-2','attP-1','PConst','attP-2','GFP'], [-1,1,-1,-1,-1,-1,-1,1])
	lu2a = Strand(['attB-1-GG','attB-1-AA','attB-2','attP-1-AA','attP-2','attP-1-GG'],[1,1,1,-1,1,-1])
	init = DNAState([lu2a])
	#init = DNAState([victoria])
	#init = DNAState([andreGenome, andreFuel], genomic = [0])

	print('---- init ----')
	print(init)
	print('---- possible states ----')
	getStates([init]);
	#[print(x) for x in getStates([init])];
	return;

if __name__ == '__main__':
	main()
