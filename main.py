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

def set_chr_range():
	d, column = {}, []
	with open('./../data/assignment/india_2009_2013.snp') as file:
		for line in file:
			column.append(int(line.split()[1]))
	for i in range(1,24):
		d[i] = []
		d[i].append(column.index(i))
	for i in range(1,23):
		d[i].append(d[i+1][0]-1)
	return d

# Gets chromosome number for an individual
def get_chrom():
	pass

# Computes $\sum_{i}X_i$ and $n_1$ for a group for each locus
def estimate(group):
	d = {'X': [], }
	with open('./../data/assignment/india_2009_2013.geno') as file:
		cols = list(map(get_column, locations[group]))

print set_chr_range()