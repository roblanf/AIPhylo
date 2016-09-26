from simulator import different_tree_simulate_alignment
from ml import *

'''
       all_trees are all of the possible unrooted tree for species
       tree_information_list is a list contain different alignment information
       tree_information is a list with parameters p,q,r,s,t
	
        p: branch length p in figure 2b (long branch)
	q: branch length q in figure 2b (short branch)
	r: branch length r in figure 2b (internal branch)
	s: length of alignment (s is for sites)
	t: tau in figure 2b. Amount of heterotachy [0, 1] = [none, lots]

	the example shows alignment with frist 20 sequences follow (a,c),(b,d)
	then 100 sequences follow (a,b),(c,d)
	and then 10 sequences follow (a,d),(b,c)

	note: I think it's better to just add all_trees in the test.
	if we just want to do experiment with 4 type of species
	then it's no difference. Otherwise we can change it in test
'''
        
all_trees = ['((a:%f, b:%f):%f,(c:%f,d:%f):%f);',
	    '((a:%f, c:%f):%f,(b:%f,d:%f):%f);',
	    '((a:%f, d:%f):%f,(b:%f,c:%f):%f);'
	    ]
tree_information_list=[]
tree_information1=[0.1, 0.5, 0.1,10, 0.0, 1]
tree_information_list.append(tree_information1)
tree_information2=[0.1, 0.5, 0.1,100,0.0, 0]
tree_information_list.append(tree_information2)
tree_information3=[0.1, 0.5, 0.1,10 ,0.0, 2]
tree_information_list.append(tree_information3)

aln,tree=different_tree_simulate_alignment(tree_information_list,all_trees)

print aln
print tree.asciiArt()
print ml4(aln,tree)
