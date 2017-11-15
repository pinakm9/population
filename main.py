import numpy
import time

# Timing wrapper
def timer(func):
	def new_func(*args,**kwargs):
		start = time.time()
		val = func(*args,**kwargs)
		end = time.time()
		print('Time taken by function {} is {} seconds'.format(func.__name__, end-start))
		return val
	return new_func

# Locates each population group in the geno file
def locate_groups():
	d, d_ = {}, {}
	with open('./../data/assignment/india_2009_2013.unrelated.ind') as file:
		for i, line in enumerate(file):
			words = line.split()
			d_[i] = words[2]
			if 'Ignore' not in words:
				if words[2] not in d:
					d[words[2]] = []
				d[words[2]].append(i)
	return d, d_

locations, loc = locate_groups()
groups = locations.keys()

# Finds stop-point to avoid chromosome 23 
def filter23():
	col = []
	with open('./../data/assignment/india_2009_2013.snp') as file:
		for i, line in enumerate(file):
			if line.split()[1] == '23':
				break
	return i

stop23 = filter23()

class Group:

	def __init__(self, name, lines):
		self.name = name
		self.cols = lines
		self.rows = ['' for i in range(stop23)]
		self.n = len(self.cols)
		self.p = {}

	@timer # Computes $\sum_{i}X_i$ and $n_1$ for a group for each locus
	def estimate(self):
		for i in range(stop23):
			self.p[i] = 1 - self.rows[i].count('0')/float(self.n)
 
G = {}
for g in groups:
	G[g] = Group(g, locations[g])

@timer
def set_data():
	with open('./../data/assignment/india_2009_2013.geno') as file:
		for j, line in enumerate(file):
			if j == stop23:
				break
			for i in range(378):
				if loc[i] != 'Ignore':
					G[loc[i]].rows[j] += line[i] 


@timer # Extracts c-th column from geno file
def get_column(c):
	col = ''
	with open('./../data/assignment/india_2009_2013.geno') as file:
		for line in file:
			col += line[c]
	return col

# Variance across
def Nl(g1, g2, l):
	x1, x2 = g1.p[l], g2.p[l]
	n1, n2 = g1.n, g2.n
	#print(x1, x2, n1, n2) 
	return (x1-x2)**2 - x1*(1-x1)/(n1-1) - x2*(1-x2)/(n2-1)

# Total variance
def Dl(g1, g2, l):	
	x1, x2 = g1.p[l], g2.p[l] 
	return (x1-x2)**2 + x1*(1-x1) + x2*(1-x2)

# Fixation index
def fixation(g1, g2):
	N, D = 0.0, 0.0
	for l in range(stop23):
		N += Nl(g1, g2, l)
		D += Dl(g1, g2, l)
	return N/D

set_data()

for g in groups:
	G[g].estimate()

with open('fixation.txt', 'w') as file:
	for i, g1 in enumerate(groups):
		for j, g2 in enumerate(groups):
			if i > j:
				s = '{}\t{}\t{}\n'.format(G[g1].name, G[g2].name, fixation(G[g1], G[g2]))
				print(s)
				file.write(s)

with open('fixation.txt') as file:
	s = 0
	for line in file:
		s += float(line.split()[2])
	print('Average fixation = {}'.format(s/float(26*51)))