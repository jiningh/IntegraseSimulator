import random
#Plasmid or Genomic Strand Object
class Strand:
	#Strand Elements, Tuple, loaded on init, immutable object
	#Possible elements:
	#	Genes in the allowableGenes list
	#	Promoters Terminators, Origins, in their respective lists
	#	Attachment sites in form of attB-bxb1-GG or attB-bxb1, where bxb1 is the 
	#   integrase and GG is the dinucleotide sequence (optional) 
	se = None
	#Direction vector, same length as strand elements, -1 for back, 1 for forward, 
	#0 for things that don't matter (attR, attL sites, origins)
	di = None
	#Is this sequence a genome, specifically does it contain a GENOME origin
	genome = None
	#After determining the origin of this strand once, the origin is stored here
	storedorigin = None
	#Things that belong to specific categories, feel free to append this list
	allowableGenes = ['GFP','YFP','RFP','MCHERRY','XIS','BFP']
	allowablePromoters = ['PCONST']
	allowableTerminators = ['TERM']
	allowableOrigins = ['COLE1', 'GENOME']
	#List of genes that can be transcribed from this strand, after it is generated once, it is kep here
	openForTranscription = None
	#Allows for "advanced" origins, for example if this is set to {'COLE1':['COL','attB-1','E1']}, the logic
	#will look for 'COL','attB-1','E1' in strand elements and if it is found, it will impose a COLE1 origin
	originSeq = None
	#Load all the things, specifically strand elements and direction
	def __init__(self, s, d, originSeq = {}):
		if len(s) !=  len(d): raise RuntimeError('Strand element length not equal to direction length')
		self.se = tuple(s); self.di = tuple(d); self.originSeq = originSeq
	#Equivalency is definately useful
	def __eq__(self, other):
		return (self.se == other.se) and (self.di == other.di)
	#Comparison method for Cython
	#def __richcmp__(self, other, int op):
	#	if op == 0: return (self.se, self.di) < (other.se, other.di) 
	#	elif op == 1: return (self.se, self.di) <= (other.se, other.di)
	#	elif op == 2: return (self.se, self.di) == (other.se, other.di)
	#	elif op == 3: return (self.se, self.di) != (other.se, other.di)
	#	elif op == 4: return (self.se, self.di) > (other.se, other.di)
	#	elif op == 5: return (self.se, self.di) >= (other.se, other.di)
	#Making it all look pretty
	def __str__(self):
		if len(self.originSeq) == 0:
			return 'Strand(' + str(list(self.se)) + ', ' + str(list(self.di)) + ')'
		return 'Strand(' + str(list(self.se)) + ', ' + str(list(self.di)) + ', originSeq = ' + str(self.originSeq) + ')'
	def __repr__(self):
		if len(self.originSeq) == 0:
			return 'Strand(' + str(list(self.se)) + ', ' + str(list(self.di)) + ')'
		return 'Strand(' + str(list(self.se)) + ', ' + str(list(self.di)) + ', originSeq = ' + str(self.originSeq) + ')'
	#Hashing is required to be a key in a dictionary
	def __hash__(self):
		return hash((self.se, self.di))
	#Find the origin of this sequence, will return the first origin it finds. Multiple origin plasmids not supported
	def origin(self):
		if self.storedorigin is not None: return self.storedorigin
		if 'GENOME' in [i.upper() for i in self.se]: return 'GENOME'
		#go through and try to find the origin
		for i in self.se:
			if i.upper() in self.allowableOrigins: 
				self.storedorigin = i.upper()
				return i.upper()
		#if nothing found impose the originSeq rules
		for key, seq in self.originSeq.items():
			i = -1
			while True:
				if self.se[i + 1:len(self.se) - len(seq) + 1].count(seq[0]) == 0: break
				i = self.se.index(seq[0], i + 1, len(self.se) - len(seq) + 1)
				if seq == self.se[i:i+len(seq)]: 
					self.storedorigin = key
					return key
		return None
	#Is this strand a genomic strand (does it contain a genomic origin)
	def isGenome(self):
		if self.genome is not None: return self.genome
		return self.origin() == 'GENOME'
	#Given position i in strand elements, check if site is P site or B site
	def checkPB(self, i):
		if self.se[i].find('attB') == 0: return 'b'
		if self.se[i].find('attP') == 0: return 'p'
		return None
	#Given positions of attachment sites at selr and serr, 
	#excise the strand, returns a list of 2 strands that this strand turns into
	def excise(self, selr, serr):
		sel = min(selr, serr); ser = max(selr, serr);
		left = 'attL-' if self.checkPB(sel) == 'b' else 'attR-'
		right = 'attR-' if self.checkPB(ser) == 'p' else 'attL-'
		id1 = self.getID(sel)
		di1 = self.getDI(sel)
		se = list(self.se); di = list(self.di);
		s1 = Strand(se[0:sel] + [left + id1 + di1] + se[ser+1:], di[0:sel] + [0] + di[ser+1:], originSeq = self.originSeq)
		s2 = Strand(se[sel+1:ser] + [right + id1 + di1], di[sel+1:ser] + [0], originSeq = self.originSeq)
		return [s1, s2]
	#Given the positions of attachment sites at selr and serr,
	#flip the strand, returns one strand
	def flip(self, selr, serr):
		se = list(self.se); di = list(self.di);
		sel = min(selr, serr); ser = max(selr, serr);
		id1 = self.getID(sel)
		di1 = self.getDI(sel)
		left = 'attL-' if self.checkPB(sel) == 'b' else 'attR-'
		right = 'attR-' if self.checkPB(ser) == 'p' else 'attL-'
		return Strand(se[0:sel] + [left + id1 + di1] + se[sel+1:ser][::-1] + [right + id1 + di1] + se[ser+1:], \
				di[0:sel] + [0] + [x*-1 for x in di[sel+1:ser][::-1]] + [0] + di[ser+1:], originSeq = self.originSeq)
	#Given another strand and position sel on this strand and ser on the other strand
	#integrate the 2 strands together
	def integrate(self, other, sel, ser):
		sse = list(self.se); sdi = list(self.di); ose = list(other.se); odi = list(other.di);
		id1 = self.getID(sel)
		di1 = self.getDI(sel)
		left = 'attL-' if self.checkPB(sel) == 'b' else 'attR-'
		right = 'attR-' if other.checkPB(ser) == 'p' else 'attL-'
		return [Strand(sse[0:sel] + [left + id1 + di1] + ose[ser+1:] + ose[0:ser] + [right + id1 + di1] + sse[sel+1:],\
			sdi[0:sel] + [0] + odi[ser+1:] + odi[0:ser] + [0] + sdi[sel+1:], originSeq = self.originSeq)]
	#Get all genes in strand elements (not just transcribable genes)
	def getAllGenes(self):
		return list(set(self.se).intersection(self.allowableGenes))
	#Get all transcribable genes (promoter followed by gene n elements down stream without terminator in between)
	def getTranscription(self):
		if self.openForTranscription is not None: return self.openForTranscription
		state = 0
		out = []
		c = -1
		for i in self.se:
			c += 1
			if self.di[c] != 1: 
				if i.upper() in self.allowableGenes and state == 1: state = 0
				continue
			if i.upper() in self.allowablePromoters and state == 0: state = 1
			if i.upper() in self.allowableTerminators and state == 1: state = 0
			if i.upper() in self.allowableGenes and state == 1:
				state = 0; out.append(i)
		c = len(self.se)
		state = 0
		for i in self.se[::-1]:
			c -= 1
			if self.di[c] != -1: 
				if i.upper() in self.allowableGenes and state == 1: state = 0
				continue
			if i.upper() in self.allowablePromoters and state == 0: state = 1
			if i.upper() in self.allowableTerminators and state == 1: state = 0
			if i.upper() in self.allowableGenes and state == 1:
				state = 0; out.append(i)
		self.openForTranscription = out
		return self.getTranscription()
	#Get all integrases that can interact with attachment sites on this strand
	def integrase(self):
		out = []
		for i in range(len(self.se)):
			if self.checkAtt(i): out.append(self.getID(i))
		return out
	#Given attachment site at position i, get intergrase which interacts with it
	def getID(self, i):
		return self.se[i].split('-')[1]
	#Check if attachment site at i and j interact with the same integrase and have
	# the same dinucleotide site
	def checkSameInt(self, i, j):
		return self.getID(i) == self.getID(j) and self.getDI(i) == self.getDI(j)
	#check if attachment site at i on this strand and j on other strand interact with the same
	# integrase and dinucleotide site
	def checkSameInt(self, other, i, j):
		return self.getID(i) == other.getID(j) and self.getDI(i) == other.getDI(j)
	#check if site i is an attachment site, or a b site or a p site
	def checkAtt(self, i, bp = 'BOTH'):
		if bp.upper() == 'BOTH': return self.se[i].find('att') == 0
		if bp.upper() == 'B': return self.se[i].find('attB') == 0
		if bp.upper() == 'P': return self.se[i].find('attP') == 0
	#get the dinueleotide configuration at site i
	def getDI(self, i):
		if len(self.se[i].split('-')) == 2: return ''
		return '-' + self.se[i].split('-')[2]
	#List all the attachment site B and site P positions on this strand
	def listBP(self):
		attP = []; attB = [];
		for i in range(len(self.se)):
			if self.checkAtt(i, bp = 'p'): attP.append(i)
			if self.checkAtt(i, bp = 'b'): attB.append(i)
		return attP, attB
	#return the direction of strand element i (string)
	def dir(self, i):
		return self.di[i]
	#return strand element i (string)
	def str(self, i):
		return self.se[i]
	#count how many of ele strand elements are there. If ele is none, return size of strand elements vector
	def count(self, ele):
		if ele is None: return len(self.se)
		return self.se.count(ele)

#Current state of the DNA system
class DNAState:
	#vector of strand objects
	strands = None
	#Number of plasmid species present given an origin, form of {'COLE1':5} (5 strands with COLE1 origin)
	plasmidSpecies = None
	#initalize
	def __init__(self, s):
		self.strands = list(set(s))
		self.plasmidSpecies = {}
	#output own list of strands
	def gets(self): return self.strands
	#must be hashable to be used as a dictionary key
	def __hash__(self):
		return hash((tuple(self.strands)))
	#equivalency is dependent on strands vector
	def __eq__(self, other):
		return (self.strands == other.strands)
	#If using cython, can use this function instead of __eq__
	#def __richcmp__(self, other, int op):
	#	if op == 0: return self.strands < other.strands
	#	elif op == 1: return self.strands <= other.strands
	#	elif op == 2: return self.strands == other.strands
	#	elif op == 3: return self.strands != other.strands
	#	elif op == 4: return self.strands > other.strands
	#	elif op == 5: return self.strands >= other.strands
	#String representation of object
	def __str__(self):
		return 'DNAState(' + str(self.strands) + ')'
	def __repr__(self):
		return 'DNAState(' + str(self.strands) + ')'

	#Get all possible genes (regardless of its on) in this DNA system
	def getAllGenes(self):
		out = []
		for i in self.strands:
			out += i.getAllGenes()
		return list(set(out))

	#Output a list of deltas for plasmid regeneration and degradation
	#Action items or Deltas are of the form [['action',item,[(item,-1),]]]
	#Action can be pdeg, pregen (plasmid degrade and regeneration), and the plasmid
	#will either increase in count or decrease in count
	def regendegrade(self):
		states = {}; actions = []
		for i in self.strands:
			if i.isGenome(): continue
			actions.append(['pdeg',i,[(i,-1)]])
			actions.append(['pregen',i,[(i,1)]])
			nStrands = self.strands[:]
			nStrands.remove(i)
			states.update({i:DNAState(nStrands)})
		return states, actions

	#given that strands in the strandlist are at a zero count, reduce the DNA system to
	#not include those strands anymore
	def reduce(self, strandlist):
		nStrands = self.strands[:]
		for i in strandlist:
			nStrands.remove(i)
			if i.isGenome(): raise RuntimeError('Reducing a genomic strand')
		return DNAState(nStrands)

	#Given strand at i[0] at position i[1] and strand at j[0] at position j[1]
	#try and execute any integrase actions
	def flip(self, i, j, bimole = False):
		sl, sel = i; sr, ser = j;
		strands = self.strands;
		if not strands[sl].checkAtt(sel) or not strands[sr].checkAtt(ser): 
			raise RuntimeError('flipping strands that are not attachment sites');
		if not strands[sl].checkSameInt(strands[sr], sel, ser): 
			raise RuntimeError('not same integrase for both att sites')
		id1 = strands[sl].getID(sel); 
		actionStrand = []; action = ''
		nStrands = strands[:];
		#action id uni is unimolecular action, bi is bimolecular action, action list is as follows
		#['actionID (uni or bi)','integrase ID',[(strand, delta),]]
		if sl == sr and not bimole:
			action = 'uni'
			if self.strands[sl].dir(sel) == strands[sr].dir(ser):
				if strands[sl].isGenome():
					nStrands = strands[0:sl] + strands[sl+1:]
				nStrands += strands[sl].excise(sel, ser)
				actionStrand.append((nStrands[-1], 1))
				actionStrand.append((nStrands[-2], 1))
				actionStrand.append((strands[sl], -1))
			else:
				nStrands[sl] = strands[sl].flip(sel, ser)
				actionStrand.append((nStrands[sl], 1))
				actionStrand.append((strands[sl], -1))
				if not strands[sl].isGenome():
					nStrands.append(strands[sl]); 
		else:
			action = 'bi'
			if strands[sl].dir(sel) == strands[sr].dir(ser):
				if (strands[sl].isGenome()) and (strands[sr].isGenome()):
					if sl == sr: return None
					s = min(sl, sr)
					l = max(sl, sr)
					nStrands = strands[0:s] + strands[s+1:l] + strands[l+1:]
				elif (strands[sl].isGenome()) or (strands[sr].isGenome()):
					u = sl if strands[sl].isGenome() else sr
					nStrands = strands[0:u] + strands[u+1:]
				nStrands += strands[sl].integrate(strands[sr], sel, ser)
				actionStrand.append((nStrands[-1], 1))
				actionStrand.append((strands[sr], -1))
				actionStrand.append((strands[sl], -1))
			else: return None
		#Return new DNAState, as well as its associated action list
		return DNAState(nStrands), [action, id1, actionStrand]

	#Given this DNA system, explore all neighboring DNA systems
	#Debug allows one to only execute reactions by a specific integrase
	#flip only will only execute unimolecular strands
	#Note that plasmids with a non-genomic origin will produce a DNA state
	#where a new strand could exist but the previous strand does not disappear
	#due to possible multiple copy numbers
	def explore(self, debug = None, fliponly = False):
		strands = self.strands
		states = []; attB = []; attP = []; actions = []
		for i in range(len(strands)):
			tP, tB = strands[i].listBP()
			attB += [(i, j) for j in tB]; attP += [(i, j) for j in tP];
		for i in attB:
			for j in attP:
				if strands[i[0]].checkSameInt(strands[j[0]],i[1],j[1]):
					if debug is not None:
						if strands[i[0]].getID(i[1]) != debug:
							continue
					if fliponly:
						if i[0] != j[0]: continue
					if i[0] == j[0]:
						out = self.flip(i,j,bimole = True)
						if out is not None:
							state, a = out
							#print(state, self)
							#if state != self and state not in states:
							states.append(state)
							actions.append(a)
					out = self.flip(i,j)
					if out is not None:
						state, a = out
						#if state != self and state not in states:
						states.append(state)
						actions.append(a)
		return states, actions;
	#Get all genes that is open to transcription in all plasmids in this state
	def getTranscription(self):
		out = []
		for i in self.strands:
			out += i.getTranscription()
		return out
	#count the number of plasmid species given an input origin
	def countPlasmidSpecies(self, origin):
		if origin not in self.plasmidSpecies:
			out = 0
			for i in self.strands:
				if i.origin() == origin:
					out += 1
			#return out
			self.plasmidSpecies.update({origin:out})
		#print(self.plasmidSpecies)
		return self.plasmidSpecies[origin]

#Parameter set, subclass of model, change parameters only in this class
class ParameterSet:
	#integrase parameters, {'integrase':params}, for params see template in init function
	pI = {}
	#Gene parameters, {'gene':params}, for params see genetemplate in init function
	genes = {}
	#bimolecular reaction rate given kflip, kD of integrase and n number of integrases
	diffStrand = lambda s, kflip, n, Kd: kflip*n**2*(n-1)**2/(Kd**4+2*Kd**3*n+2*Kd**2*n*(n-1)+Kd**2*n**2+2*Kd*n**2*(n-1)+n**2*(n-1)**2)
	#unimolecular reaction rate given the same parameters as above
	sameStrand = lambda s, kflip, n, Kd: kflip*n*(n-1)*(n-2)*(n-3)/(Kd**4+Kd**3*n+Kd**2*n*(n-1)+Kd*n*(n-1)*(n-2)+n*(n-1)*(n-2)*(n-3))
	#origin regeneration parameters for plasmids
	origins = {'COLE1':{'kprod':50.,'kd':100}}
	#plasmid degradation parameter
	plasmiddeg = 0.3
	#plasmid production rate given plasmid parameters and n number of plasmids
	pprode = lambda s, kprod, kd, n: kprod*(n/kd)/(1+(n/kd))
	#plasmid degradation rate given plasmid parameter and n counts
	pdege = lambda s, kdeg, n: kdeg*n
	#loading default parameters
	def __init__(self, integrases, genes):
		#active is when the inducer is on, it can either be a lambda, or a list of tuples bopunding on times
		#ex. [(0, 10),(20,50)]
		template = {'kflip1':0.4, 'kflip2':0.4, 'kprod':50, 'kleak':0.00*50, 'kdeg':0.3, 'Kd':10, 'active':lambda t: 1}
		genetemplate = {'kprod':50, 'kdeg':0.3}
		for i in integrases:
			self.pI.update({i:template.copy()})
		for i in genes:
			self.genes.update({i:genetemplate.copy()})
	#plasmid production given origin, number of plasmid species
	def pprod(self, n, origin, nspecies):
		if origin is None: return 0.0
		return self.pprode(self.origins[origin]['kprod']/float(nspecies), self.origins[origin]['kd']/float(nspecies), n)
		#return self.pprode(self.origins[origin]['kprod'], self.origins[origin]['kd'], n)
	#plasmid degradation
	def pdeg(self, n):
		return self.pdege(self.plasmiddeg, n)
	#gene production
	def gprod(self, gene): return self.genes[gene]['kprod']
	#gene degradation
	def gdeg(self, gene, n): return self.genes[gene]['kdeg'] * n
	#integrase production
	def prod(self, integ, time):
		return self.pI[integ]['kprod']*self.induce(integ, time) + self.pI[integ]['kleak']
	#integrase degradation
	def deg(self, integ, n):
		return self.pI[integ]['kdeg'] * n
	#unimolecular integrase interactions
	def uni(self, integ, n, s1):
		out = self.sameStrand(self.pI[integ]['kflip1'], n, self.pI[integ]['Kd']) * s1
		return out
	#bimolecular integrase interactions
	def bi(self, integ, n, s1, s2):
		out = self.diffStrand(self.pI[integ]['kflip2'], n, self.pI[integ]['Kd']) * s1 * s2
		return out
	#helper function to check if inducer is a lambda
	def checkLambda(self, f):
		c = lambda:0
		return isinstance(f, type(c)) and c.__name__ == f.__name__
	#helper function to work on if the inducer is on for a specific integrase at a specific time
	def induce(self, integ, time):
		if self.checkLambda(self.pI[integ]['active']): return self.pI[integ]['active'](time)
		for i in self.pI[integ]['active']:
			if i[0] <= time < i[1]: return True
		return False

#Counter object for Gillespie simulation
class DNAStateCount:
	#time vector
	time = [0]
	#Current State, should be a DNAState object
	curState = None
	#Integrase counts {'integrase'[timeseries],}
	integrase = {}
	#Gene counts, same format as integrase counts
	genes = {}
	#Strand counts, same format as integrase counts
	strands = {}
	#If there are no more available states. In the case of system without genes
	#This is when there is no more neighboring integrase states
	#In system with genes, this is where the sum of propensities is zero
	noMoreStates = False
	#Initialize with strand counts, and if a gene vector should be loaded
	def __init__(self, strandCounts, simgene):
		self.time = [0]; self.curState = None; self.integrase = {}; 
		self.noMoreStates = False;
		self.strands = {}; self.genes = {}; self.integrase = {}
		pI = []; genes = []
		#strandCounts in form of [[strand, copies],]
		for i in range(len(strandCounts)):
			self.strands.update({strandCounts[i][0]:[strandCounts[i][1]]})
			pI += strandCounts[i][0].integrase()
			genes += strandCounts[i][0].getAllGenes()
		self.curState = DNAState(list(self.strands.keys()))
		for j in pI: self.integrase.update({j:[0]})
		if simgene: 
			for j in genes: self.genes.update({j:[0]})
			
	#(act, int, [[strand, +-x]])
	#Given an action and a dt, update the state of the counter
	def update(self, action, dt):
		strands = self.strands
		self.time.append(self.time[-1] + dt)
		for k in self.integrase: self.integrase[k].append(self.integrase[k][-1])
		for k in self.genes: self.genes[k].append(self.genes[k][-1])
		for k in self.strands: self.strands[k].append(self.strands[k][-1])
		if action[0] == 'int':
			self.integrase[action[1]][-1] = self.integrase[action[1]][-1] + action[2]
			return
		if action[0] == 'gene':
			self.genes[action[1]][-1] = self.genes[action[1]][-1] + action[2]
			return
		#print(action)
		for strand in self.curState.strands:
			if strand not in strands: 
				strands.update({strand:[0]*(len(self.time))})
			if len(strands[strand]) != len(self.time): raise RuntimeError('Updater broken')
		for a in action[2]:
			strands[a[0]][-1] += a[1]

	#gives the count of strand at last time point, outputs int
	def count(self, strand):
		return self.strands[strand][-1]

	#Get the state of everything at the last timepoint, outputs {strand:count,}
	def getend(self):
		out = {}
		for key, val in self.strands.items():
			out.update({key:val[-1]})
		return out

	#gets the state of everything at last timepoint, outputs {gene:count,}
	def getendgene(self):
		out = {}
		for key, val in self.genes.items():
			out.update({key, val[-1]})
		for key, val in self.integrase.items():
			out.update({key, val[-1]})
		return out

	#Get the full trajectory given a linearily or logspace timeseries
	#automatically convert times generated by Gillespie to the nice formated timeseries inputted
	#outputs {strand or gene or integrase:[timeseries],}
	def gettraj(self, t):
		outstrand = {}; outgene = {}; outint = {};
		for i in self.strands: outstrand.update({i:[]})
		for i in self.genes: outgene.update({i:[]})
		for i in self.integrase: outint.update({i:[]})
		timeindex = 0;
		#print(self.strands)
		for i in range(len(t)):
			while True:
				if timeindex == len(self.time) - 1: break
				if self.time[timeindex] > t[i]: raise RuntimeError('Timer circuit broken')
				if self.time[timeindex] <= t[i] and self.time[timeindex + 1] > t[i]: break
				else: timeindex += 1
			for s, ts in self.strands.items(): 
				outstrand[s].append(self.strands[s][timeindex])
			for k, ts in self.genes.items(): outgene[k].append(self.genes[k][timeindex])
			for k, ts in self.integrase.items(): outint[k].append(self.integrase[k][timeindex])
		return outstrand, outgene, outint

#Main model class, contains Gillespie algorithm
class Model:
	#dsc is the counter (previous class), param is the parameter set class (2nd previous class)
	#the cache is the explore variable, {DNAState:(states, actions (corresponding to states), rdstates (states given plasmid degradation if degradation reaches 0), plasmid degradation and regen actions)}
	dsc = None; param = None;
	explore = {}; initCondition = None
	#ran will be True if the dsc needs to be refreshed, simgene is if genes are simulated
	ran = False; simgene = False
	#initialize strandCounts (see previous class, same format), and if genes are replicated
	def __init__(self, strandCounts, simgene = True):
		self.explore = {}
		self.initCondition = strandCounts
		self.precompute = False; self.dsc = None; self.param = None;
		self.ran = False
		pI = []; genes = []
		for i in range(len(strandCounts)):
			pI += strandCounts[i][0].integrase()
			genes += strandCounts[i][0].getAllGenes()
		self.simgene = simgene
		self.param = ParameterSet(list(set(pI)), list(set(genes)))
		self.refresh()

	#Refreshes the model, and recycles the counter object
	def refresh(self):
		self.dsc = None
		del self.dsc
		a = DNAStateCount(self.initCondition[:], self.simgene)
		self.dsc = a
		self.ran = False

	#It's easy to precompute things, this was implemented in utilities rather than here
	def precompute(self): 
		raise NotImplemented('Precompute not yet supported')

	#Take a Gillespie step
	def gillespie_step(self):
		#fetch all possible Integrase actions, if current state is not in the cache, explore states
		if self.dsc.curState not in self.explore:
			states, action = self.dsc.curState.explore()
			rdstates, rdaction = self.dsc.curState.regendegrade()
			self.explore.update({self.dsc.curState:(states, action, rdstates, rdaction)})
		states, stateaction, rdstates, rdactions = self.explore[self.dsc.curState]
		
		#if no neighboring state exist and genes are not simulated, stop the simulation
		if len(states) == 0 and self.simgene == False: 
			self.dsc.noMoreStates = True
			return
		
		#get propensities
		actions, propensity = self.gillespie_propensity(stateaction[:], rdactions[:])
		
		#if no more reactions, stop the simulation
		if sum(propensity) == 0: 
			self.dsc.noMoreStates = True
			return

		#draw a time and an action to perform
		dt, aid = self.draw(propensity)
		#figure out if state needs to change (if action is bi or uni it probably will need to change)
		#It wont need to change if only plasmids are integrated
		if aid < len(states):
			self.dsc.curState = states[aid] 
		#perform the update given dt and the action
		self.dsc.update(actions[aid], dt)
		#figure out which plasmids reached 0 count, and remove them from the state
		removelist = []
		for strand in self.dsc.curState.strands:
			if self.dsc.count(strand) == 0: removelist.append(strand)
		self.dsc.curState = self.dsc.curState.reduce(removelist)

	#draw 2 numbers for dt and which action to take
	def draw(self, propensity):
		sp = sum(propensity)
		dt = random.expovariate(sp)#np.random.exponential(1.0/sp)
		t = 0
		r = random.uniform(0, sp)#np.random.rand() * sp
		for i in range(len(propensity)):
			t += propensity[i] 
			if t >= r: break
		return dt, i

	#get propensity given a list of DNA actions and plasmid degradation/regeneration actions
	def gillespie_propensity(self, DNAStateAction, rdactions):
		actions = DNAStateAction; propensity = []
		for i in actions:
			loss = [j[0] for j in i[2] if (j[1] == 0 or j[1] == -1)]
			if len(loss) > 2: raise RuntimeError('More than 2 decreasing states')
			if i[0] == 'uni': 
				propensity.append(self.param.uni(i[1], self.dsc.integrase[i[1]][-1], self.dsc.count(loss[0])))
			if i[0] == 'bi': 
				propensity.append(self.param.bi(i[1], self.dsc.integrase[i[1]][-1], self.dsc.count(loss[0])-int(loss[0] == loss[1]), self.dsc.count(loss[1])))
		#print(propensity)
		actions += rdactions
		for i in rdactions:
			if i[0] == 'pdeg':
				#if self.param.pdeg(self.dsc.count(i[1])) != 0: print(i[1])
				propensity.append(self.param.pdeg(self.dsc.count(i[1])))
			if i[0] == 'pregen':
				#if self.param.pprod(self.dsc.count(i[1]), i[1].origin()) != 0: print(self.param.pprod(self.dsc.count(i[1]), i[1].origin()))
				propensity.append(self.param.pprod(self.dsc.count(i[1]), i[1].origin(), self.dsc.curState.countPlasmidSpecies(i[1].origin())))
		for i in self.dsc.integrase:
			actions += [('int', i, 1), ('int', i, -1)]
			propensity += [self.param.prod(i, self.dsc.time[-1]), self.param.deg(i, self.dsc.integrase[i][-1])]
		for i in self.dsc.genes:
			actions += [('gene', i, 1), ('gene', i, -1)]
			availGene = self.dsc.curState.getTranscription()
			#print(i)
			propensity += [self.param.gprod(i)*int(i in availGene), self.param.gdeg(i, self.dsc.genes[i][-1])]
		#[print(c) for c in actions]
		#print(propensity)
		return actions, propensity

	#run the gillespioe code for a certain amount of time
	def run_gillespie(self, t):
		if self.ran: raise RuntimeError('Needs to refresh Model')
		while self.dsc.time[-1] < t and not self.dsc.noMoreStates:
			self.gillespie_step()
			#print(self.dsc.time[-1])
		self.ran = True