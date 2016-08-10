# An alignment generator using PyCogent
# It's not really a generator yet, but it has the basics

# the approach follows Phillipe et al, in that we generate two 
# alignments on trees with matching topologies but different branch 
# lenths. 

from cogent.core import alignment, tree
from cogent.evolve import substitution_model
from cogent.parse.tree import DndParser
from cogent.seqsim.tree import RangeNode
from cogent import LoadTree

# parameters correspond to figure 2 of Phillipe et al 2005: BMC Evo Bio 5:50
p = 1
q = 0.1 
t = 0.5 # the 'amount of heterotachy', 0 = none, 1 = a lot
r = 0.1 # the internal branch
w = 0.5 # amount of data from each tree. Phillipe et al always have w = 0.5

# define our alignment lengths according to w
# NB phillippe et al use an alignment length of 10000
total_alignment_length = 10000
t1_length = int(round(total_alignment_length * w))
t2_length = total_alignment_length - t1_length

# these are the branch lengths of our two trees, see fig 2 of Phillipe et al
t1_bl1 = (1 + t) * p
t1_bl2 = (1 - t) * q
t2_bl1 = (1 - t) * p
t2_bl2 = (1 + t) * q

# build our trees
# we use r/2 so we can root the tree and define where evolution begins,
# the inference will work on unrooted trees though.
# Note that here we assume the tree structure is ((a,b),(c,d))
# we will presumably have to randomise this to train an AI, otherwise
# it's a very simple problem!
t1_string = '((a:%f, b:%f):%f,(c:%f,d:%f):%f)root;' %(t1_bl1, t1_bl2, r/2, t1_bl1, t1_bl2, r/2)
t2_string = '((a:%f, b:%f):%f,(c:%f,d:%f):%f)root;' %(t2_bl1, t2_bl2, r/2, t2_bl1, t2_bl2, r/2)

t1 = LoadTree(treestring = t1_string)
t2 = LoadTree(treestring = t2_string)

# simulte our two alignments using a Jukes Cantor model
# this assumes that all state frequences are equal at 0.25, and 
# that transition rates are also equal (e.g. A<->T == C<->G etc)

sm = substitution_model.Nucleotide()

lf1 = sm.makeLikelihoodFunction(t1)
lf1.setConstantLengths()
aln1 = lf1.simulateAlignment(sequence_length = t1_length)

lf2 = sm.makeLikelihoodFunction(t2)
lf2.setConstantLengths()
aln2 = lf2.simulateAlignment(sequence_length = t2_length)


# now we just join together our two alignments
aln = aln1 + aln2