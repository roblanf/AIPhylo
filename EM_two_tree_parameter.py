from cogent import LoadTree, LoadSeqs,DNA
from cogent.evolve.models import JC69
from simulator import simulate_alignment



##input tree structure and alignment,two different p,q,r parameter
######### aln is the sequence for four species 
######### tree is the already known tree structure like '((a:%f, b:%f):%f,(c:%f,d:%f):%f);'
######### tree_parameter is a list contain two different tree parameters
######### tree_parameter=[[p1,q1,r1],[p2,q2,r2]]

##output which alignment belong to which tree
########result=[1,1,2,2,1] indicates the first site belong to tree with parameter1


def expectation_singlesite(aln,tree1,tree2):
    modle=JC69()
    result=[]
    lf1=modle.makeLikelihoodFunction(tree1)
    lf2=modle.makeLikelihoodFunction(tree2)
    for i in range(len(aln)):
        lf1.setAlignment(aln[i])
        prob1=lf1.getLogLikelihood()
        lf2.setAlignment(aln[i])
        prob2=lf2.getLogLikelihood()
        if(prob1>prob2):
            result.append(1)
        if(prob1<prob2):
            result.append(2)
        if(prob1==prob2):
            result.append(0)
    return result
        
def optimization(result,aln,tree1,tree2,tree):
    aln1=LoadSeqs(data=[('a', ''), ('c', ''), ('b', ''), ('d', '')],moltype=DNA)
    aln2=LoadSeqs(data=[('a', ''), ('c', ''), ('b', ''), ('d', '')],moltype=DNA)
    for i in range(len(aln)):
        if (result[i]==1):
            aln1=aln1+aln[i]
        if (result[i]==2):
            aln2=aln2+aln[i]
        if (result[i]==0):
            aln1=aln1+aln[i]
            aln2=aln2+aln[i]
    print "aln1:"    
    print aln1
    print "aln2:"
    print aln2
    tree_parameter=[[],[]]
    modle=JC69()
    lf1=modle.makeLikelihoodFunction(tree)
    lf1.setAlignment(aln1)
    lf1.optimise(local=True)
    p1=(lf1.getParamValue('length','a')+lf1.getParamValue('length','c'))/2.0
    q1=(lf1.getParamValue('length','b')+lf1.getParamValue('length','d'))/2.0
    r1=lf1.getParamValue('length','edge.1')+lf1.getParamValue('length','edge.0')
    lf2=modle.makeLikelihoodFunction(tree)
    lf2.setAlignment(aln2)
    lf2.optimise(local=True)
    p2=(lf2.getParamValue('length','a')+lf2.getParamValue('length','c'))/2.0
    q2=(lf2.getParamValue('length','b')+lf2.getParamValue('length','d'))/2.0
    r2=lf2.getParamValue('length','edge.1')+lf2.getParamValue('length','edge.0')
    tree_parameter[0]=[p1,q1,r1]
    tree_parameter[1]=[p2,q2,r2]
    print "tree parameter in this time"
    print tree_parameter
    return tree_parameter
    


# simulate a 1000bp alignment p1=0.19 q1=0.05,r1=r2=0.1,p2=0.01 q2=0.95
aln, tree = simulate_alignment(0.1, 0.5, 0.1, 500, 0.9, 'fixed')
print aln    
print tree.asciiArt()

known_tree='((a:%f, b:%f):%f,(c:%f,d:%f):%f);'
tree_parameter=[[1,0.01,0.5],[0.01,1,0.5]]
tree=LoadTree(treestring = '((a,b),(c,d))')

for i in range(5):
     tree1=LoadTree(treestring=known_tree%(tree_parameter[0][0],tree_parameter[0][1],tree_parameter[0][2]/2.0,tree_parameter[0][0],tree_parameter[0][1],tree_parameter[0][2]/2.0))
     tree2=LoadTree(treestring=known_tree%(tree_parameter[1][0],tree_parameter[1][1],tree_parameter[1][2]/2.0,tree_parameter[1][0],tree_parameter[1][1],tree_parameter[1][2]/2.0))
     result= expectation_singlesite(aln,tree1,tree2)
     print result
     tree_parameter=optimization(result,aln,tree1,tree2,tree)


print "final tree parameter for tree1 and tree2"
print tree_parameter

 

