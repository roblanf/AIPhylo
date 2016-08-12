# assess all possible 4 taxon trees with ML

from cogent.core import alignment, tree
from cogent.evolve import substitution_model
from cogent import LoadTree
from cogent.evolve.models import JC69

def ml4(aln, true_tree):
	'''
	Input a true tree and an alignment
	Calculate the likelihood of all possible unrooted 4-taxon trees
	Return True if the ML tree is the true tree
	Return False otherwise
	'''

	# all trees with unit branch lengths
	all_trees = [LoadTree(treestring = '((a,b),(c,d))'),
				 LoadTree(treestring = '((a,c),(b,d))'),
				 LoadTree(treestring = '((a,d),(b,c))')
				]

	# optimise lf for all trees
	sm = JC69()

	results = []
	for t in all_trees:
		lf = sm.makeLikelihoodFunction(t)
		lf.setAlignment(aln)
		lf.optimise(local=True)
		results.append(lf.getLogLikelihood())

	# get the ml tree and compare to true tree
	ml_tree = all_trees[results.index(max(results))]
	
	return ml_tree.sameTopology(true_tree)