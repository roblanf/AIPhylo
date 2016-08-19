# An alignment generator using PyCogent
# It's not really a generator yet, but it has the basics

# the approach follows Phillipe et al, in that we generate two 
# alignments on trees with matching topologies but different branch 
# lenths. 

from cogent.core import alignment, tree
from cogent.evolve import substitution_model
from cogent import LoadTree
import random
from cogent.evolve.models import JC69



def simulate_alignment(p, q, r, s, t, tree = 'random'):
	'''
	A function to simulate alignments based on Philippe et al 2005
	http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-5-50

	Input variables

	p: branch length p in figure 2b (long branch)
	q: branch length q in figure 2b (short branch)
	r: branch length r in figure 2b (internal branch)
	s: length of alignment (s is for sites)
	t: tau in figure 2b. Amount of heterotachy [0, 1] = [none, lots]
	w: omega in Philippe et al, relative weight of partitions
	tree: ['random', 'fixed']
		  'random': choose a random 4 taxon tree
		  'fixed': use the tree ((a,b),(c,d))

	Output

	A PyCogent alignment
	A PyCogent tree which is the true tree with true branch lengths

	'''

	# relative weight of two trees
	# must be 0.5 to get the branch lengths averaging nicely
	w = 0.5 

	# Define our alignment lengths according to s and w
	t1_length = int(round(s * w))
	t2_length = s - t1_length

	# Define the branch lengths of the two trees
	t1_bl1 = (1 + t) * p
	t1_bl2 = (1 - t) * q
	t2_bl1 = (1 - t) * p
	t2_bl2 = (1 + t) * q

	# choose a tree from all possible unrooted 4 taxon trees
	all_trees = ['((a:%f, b:%f):%f,(c:%f,d:%f):%f);',
				 '((a:%f, c:%f):%f,(b:%f,d:%f):%f);',
				 '((a:%f, d:%f):%f,(b:%f,c:%f):%f);'
				]
	if tree == 'random':
		tree_string = random.choice(all_trees)
	elif tree == 'fixed':
		tree_string = all_trees[0]
	else:
		raise ValueError('Unrecognised option for "tree". Check')

	# the true tree has branch lengths p and q as long as w = 0.5
	true_tree_bl = tree_string %(p, q, r/2.0, p, q, r/2.0)
	true_tree = LoadTree(treestring = true_tree_bl)

	# build our two trees
	t1 = build_tree(tree_string, t1_bl1, t1_bl2, r)
	t2 = build_tree(tree_string, t2_bl1, t2_bl2, r)

	# simulate the alignments
	a1 = get_alignment(t1, t1_length)
	a2 = get_alignment(t2, t2_length)

	aln = a1 + a2

	return(aln, true_tree)


def build_tree(tree_string, bl1, bl2, r):
	'build a PyCogent tree object from a string and branch lengths'
	# we use r/2.0 because PyCogent defaults to adding a branch of 
	# length 1 if you don't explicitly specify it
	# having 2 branches of r/2.0 keeps our internal branch at r
	tree_string_bl = tree_string %(bl1, bl2, r/2.0, bl1, bl2, r/2.0)
	t = LoadTree(treestring = tree_string_bl)
	return t

def get_alignment(tree, N_sites):
	'build a PyCogent alignment object from a tree and length'
	sm = JC69()
	lf = sm.makeLikelihoodFunction(tree)
	lf.setConstantLengths()
	aln = lf.simulateAlignment(sequence_length = N_sites)
	return(aln)


