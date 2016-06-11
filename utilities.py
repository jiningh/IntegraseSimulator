import numpy as np
from mia import *
import multiprocessing as mp
#some utilities for single stranded systems

def simulate_wrapper(strand, intseq, N, time, simgene):
	"""
	Helper function for parallel wrapper, same logic as serial wrapper
	"""
	os = []; og = [];
	a = Model([[strand, 1]], simgene = simgene)
	#unpack the integrase inducer turn on times
	for key, val in intseq.items(): a.param.pI[key]['active'] = val
	#for N simulations
	for i in range(N):
		#run simulation from model for tottime time
		a.run_gillespie(time[-1])
		#get the output trajectory in terms of {strand object, or integrase, or gene (strand object, string, string):timeseries}
		b, gene, inte = a.dsc.gettraj(time)
		os.append(b); og.append(gene)
		a.refresh()
	return [os, og]

def parallel_wrapper(strand, intseq, N, trackStrands, trackGene, time, proc = 3):
	"""
	Same wrapper as serial_wrapper (see below), except this is parallel
	"""
	if N % proc != 0: raise RuntimeError('N/proc is not an exact division')
	#setup outputs, with blank timeseries
	outgene = {}; outstrand = {}
	for i in trackGene:
		outgene.update({i: np.zeros(time.shape)})
	for i, s in trackStrands.items():
		outstrand.update({i: np.zeros(time.shape)})
	#easier to keep track of {strand object:strand name} as strand object is the output of the simulation
	invtrack = {v:k for k,v in trackStrands.items()}
	p = mp.Pool(processes = proc)
	simgene = (not len(trackGene) == 0)
	threads = [p.apply_async(simulate_wrapper, (strand, intseq, int(N/proc), time, simgene,)) for i in range(proc)]
	results = [j.get() for j in threads]
	c = 0
	for result in results:
		oss = result[0]; osg = result[1];
		c += len(oss)
		for i in oss:
			for strand, ts in i.items():
				if strand in invtrack:
					outstrand[invtrack[strand]] += np.array(ts)
		for i in osg:
			for g, ts in i.items():
				if g in outgene:
					outgene[g] += np.array(ts)

	return outstrand, outgene

def serial_wrapper(strand, intseq, N, trackStrands, trackGene, time):
	"""
	Generic Wrapper for genomic circuits
	Inputs:
		strand = strand object, initial circuit
		intseq = integrase sequence, {integrase id (string):lambda of inducer}
		N = number of simulations (int)
		trackStrands = which strands should be outputted, {strand name (arbitrary string): strand object}
		trackGene = genes to keep track of, [gene (string)]
		time = simulation time array, list of timepoints, the model itself does not inherently keep linearily spaced time, but rather the time of the Gillespie simulation
	Outputs:
		outstrand = strand timeseries output {strand name (string, as given in trackstrand):[time series as a list (int or double)]}
		outgene = gene timeseries output {gene name (string, as given in trackgene):[time series as a list (int or double)]}
	"""
	#create model object of a single genetic strand
	a = Model([[strand, 1]], simgene = (not len(trackGene) == 0))
	#unpack the integrase inducer turn on times
	for key, val in intseq.items():
		a.param.pI[key]['active'] = val
	#setup outputs, with blank timeseries
	outgene = {}; outstrand = {}
	for i in trackGene:
		outgene.update({i: np.zeros(time.shape)})
	for i, s in trackStrands.items():
		outstrand.update({i: np.zeros(time.shape)})
	#easier to keep track of {strand object:strand name} as strand object is the output of the simulation
	invtrack = {v:k for k,v in trackStrands.items()}
	#for N simulations
	for i in range(N):
		#run simulation from model for tottime time
		a.run_gillespie(time[-1])
		#get the output trajectory in terms of {strand object, or integrase, or gene (strand object, string, string):timeseries}
		b, gene, inte = a.dsc.gettraj(time)
		#for every strand if the strand is in tracked strands, add this timeseries to the total timeseries
		for strand, ts in b.items():
			if strand in invtrack:
				outstrand[invtrack[strand]] += np.array(ts)
		#for every gene if the gene is in tracked genes, add this timeseries to the total timeseries for that gene
		for g, ts in gene.items():
			if g in outgene:
				outgene[g] += np.array(ts)
		#erase the counters so we can run the simulation again
		#this entire loop can be parallelized if you want
		a.refresh()
	#for g, b in outgene.items(): 
	#	outgene[g] = b.tolist()
	#for g, b in outstrand.items(): 
	#	outstrand[g] = b.tolist()
	return outstrand, outgene

def getStates(inits, debugint):
	"""
	Get the states that is the consequence of adding a specific integrase to a specific strand
	This is helpful to get the states that are the expected outcomes in a specific sequence of inducer
	Input:
		inits = DNAStrand object of the initial strand
		debugint = string, integrase used to flip inits
	Outputs:
		DNAState object of the final state after there are no more possible integrase moves left
	"""
	#set up initial states
	pastState = inits
	while True:
		#possible output states
		out = []
		for state in pastState:
			#explore all possible states given the presense of integrase debugint
			states, actions = state.explore(debug = debugint, fliponly=False)
			#for all the states deemed available, remove strands that gets excised:
			for state in states:
				removelist = []
				#remove strands that do not have an origin
				for strand in state.strands:
					if strand.origin() is None: removelist.append(strand)
				#add to current output the reduced state (no strands with no origin)
				out.append(state.reduce(removelist))
		#if there are no more availabel states, we're done
		if len(out) == 0: break
		#if there are more available states, run until no more available states
		pastState = list(set(out));
	#should only be 1 state left standing
	if len(pastState) > 1: raise RuntimeError('Ambigous DNA States')
	#return that state
	return pastState[0]

def sequence(init, seq):
	#given a sequence of inducers, get the state at the end of the sequence of inducrs
	past = DNAState([init])
	for i in seq:
		#keep getting states given new inputs
		past = getStates([past], i)
	return past.strands

def sequences(init, seqs):
	#given a series of sequences of inducers, get the output for the series of sequence of inducers
	out = {'init':init}
	for name, seq in seqs.items():
		#essentially just runs the previous function and labels the strands with names of the sequences
		strands = sequence(init, seq)
		if len(strands) > 1: raise RuntimeError('Ambigous Strands')
		out.update({name:strands[0]})
	return out


def getAllStates(init):
	#given an initial DNA state, explore all possible states given all possible integrases
	pastState = init
	totState = []; totState += init;
	while True:
		out = []
		for state in pastState:
			states, actions = state.explore()
			out += states
			print('----' + str(state))
			[print(states[b], actions[b], states[b].getTranscription()) for b in range(len(states))]
		print('=====')
		#[print(x) for x in out]
		if len(out) == 0: break
		totState += out
		pastState = out;
		#print(out);
		if len(totState) > 10: break
	return totState

