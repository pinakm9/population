

# Creates a set of population groups
def get_groups():
	names = []
	with open('./../data/assignment/india_2009_2013.unrelated.ind') as file:
		for line in file:
			words = line.split()
			if 'Ignore' not in words:
				names.append(words[2])
	return list(set(names)) 

groups = get_groups()

# Locates each population group in the geno file
def locate_groups():
	d = {}
	for g in groups:
		d[g] = []
	with open('./../data/assignment/india_2009_2013.unrelated.ind') as file:
		for i, line in enumerate(file):
			words = line.split()
			if 'Ignore' not in words:
				d[words[2]].append(i)
	return d

locations = locate_groups()

# Extracts c-th column from geno file
def get_column(c):
	col = []
	with open('./../data/assignment/india_2009_2013.geno') as file:
		for line in file:
			col.append(int(line[c]))
	return col

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

	def __init__(self, name):
		self.name = name
		self.p, self.n = self.estimate()   

	# Computes $\sum_{i}X_i$ and $n_1$ for a group for each locus
	def estimate(self):
		d = {}
		with open('./../data/assignment/india_2009_2013.geno') as file:
			cols = list(map(get_column, locations[self.name]))
		n = len(cols)
		for i in range(stop23):
			d[i] = sum([cols[k][i] for k in range(n)])/float(n)
		return d, n

# Variance across
def Nl(g1, g2, l):
	x1, x2 = g1.p[l], g2.p[l]
	n1, n2 = g1.n, g2.n
	print(x1, x2, n1, n2) 
	return (x1-x2)**2 #- x1*(1-x1)/(n1-1) - x2*(1-x2)/(n2-1)

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

Groups = []
for g in groups:
	Groups.append(Group(g))

for i, g1 in enumerate(Groups):
	for j, g2 in enumerate(Groups):
		if i > j:
			print('{}\t{}\t{}'.format(g1.name, g2.name, fixation(g1, g2)))


