from cogent import LoadTree, LoadSeqs, DNA
from cogent.evolve.models import JC69
from simulator import simulate_alignment
import matplotlib.pyplot as plt
import numpy as np
import copy

# input tree structure and alignment,two different p,q,r parameter
# aln is the sequence for four species
# tree is the already known tree structure like '((a:%f, b:%f):%f,(c:%f,d:%f):%f);'
# tree_parameter is a list contain two different tree parameters
# tree_parameter=[[p1,q1,r1],[p2,q2,r2]]

# output which alignment belong to which tree
# result=[1,1,2,2,1] indicates the first site belong to tree with parameter1


def expectation_singlesite(aln, tree1, tree2):
    modle = JC69()
    result = []
    lf1 = modle.makeLikelihoodFunction(tree1)
    # make lilelihood function for tree1,tree2
    lf2 = modle.makeLikelihoodFunction(tree2)
    for i in range(len(aln)):  # for each site,compare the likelihood for it belong to tree1,tree2
        # and assign it to the one with larger likelihood
        lf1.setAlignment(aln[i])
        prob1 = lf1.getLogLikelihood()
        lf2.setAlignment(aln[i])
        prob2 = lf2.getLogLikelihood()

        # if(prob1 > prob2):
        #     result.append(1)
        # if(prob1 < prob2):
        #     result.append(2)
        # # if it is the same, assign to both of the two trees
        # if(prob1 == prob2):
        #     result.append(0)

        _max = max(prob1, prob2)
        prob1 -= _max
        prob2 -= _max
        exp_prob1 = np.exp(prob1)
        exp_prob2 = np.exp(prob2)
        _sum = exp_prob1 + exp_prob2
        random_prob = np.random.random() * _sum
        if exp_prob1 < random_prob:
            result.append(1)
        else:
            result.append(2)

    return result


# input two trees, sites and alignment
# aln is the sequence for four species
# tree1 is tree with branch lengths p1 q1 r1
# result is the output of expectation

# output likelihood,new parameters
# likelihood is the likelihood follow the input tree branch lengths and input sites alignment
# tree_parameters are the parameters after optimization

def optimization(result, aln, tree1, tree2):

    # get the sites for each tree according to the assignments
    aln1 = LoadSeqs(data=[('a', ''), ('c', ''),
                          ('b', ''), ('d', '')], moltype=DNA)
    aln2 = LoadSeqs(data=[('a', ''), ('c', ''),
                          ('b', ''), ('d', '')], moltype=DNA)
    for i in range(len(aln)):
        if (result[i] == 1):
            aln1 = aln1 + aln[i]
        if (result[i] == 2):
            aln2 = aln2 + aln[i]
        if (result[i] == 0):
            aln1 = aln1 + aln[i]
            aln2 = aln2 + aln[i]
    tree_parameter = [[], []]
    modle = JC69()

    # calculate the likelihood and do optimization. optimise will generates
    # new tree parameters
    lf1 = modle.makeLikelihoodFunction(tree1)
    lf1.setAlignment(aln1)
    lf1.optimise(local=True)
    likelihood1 = lf1.getLogLikelihood()


    # new tree parameters generates by optimise. As tree1/2 is symmetric, get
    # p,q,r from 6 branch lengths
    p1 = (lf1.getParamValue('length', 'a') +
          lf1.getParamValue('length', 'c')) / 2.0
    q1 = (lf1.getParamValue('length', 'b') +
          lf1.getParamValue('length', 'd')) / 2.0
    r1 = lf1.getParamValue('length', 'edge.1') + \
        lf1.getParamValue('length', 'edge.0')

    lf2 = modle.makeLikelihoodFunction(tree2)
    lf2.setAlignment(aln2)
    lf2.optimise(local=True)
    likelihood2 = lf2.getLogLikelihood()
    p2 = (lf2.getParamValue('length', 'a') +
          lf2.getParamValue('length', 'c')) / 2.0
    q2 = (lf2.getParamValue('length', 'b') +
          lf2.getParamValue('length', 'd')) / 2.0
    r2 = lf2.getParamValue('length', 'edge.1') + \
        lf2.getParamValue('length', 'edge.0')

    # return the new tree_parameter. As likelihood is in log, so plus together
    # get the total likelihood for the whole sites
    tree_parameter[0] = [p1, q1, r1]
    tree_parameter[1] = [p2, q2, r2]
    likelihood = likelihood1 + likelihood2


    return tree_parameter, likelihood


p1 = []
q1 = []
r1 = []
p2 = []
q2 = []
r2 = []  # to keep a record for the final parameters for several results--> plots may draw later
likelihood = []
start=[]

# number in range of j indicates times of experiments,for each experiment, generate new sites
##########start of several experiments#################
aln, tree = simulate_alignment(0.1, 0.5, 0.1, 500, 0.9, 'fixed')

for j in range(1):

    best_ll = -np.inf
    ##########start of one experiment#################
    tree_parameter = []
    tree1_parameter = []
    tree2_parameter = []  # tree parameter generates by random
    

    # simulate a 1000bp alignment p1=0.19 q1=0.05,r1=r2=0.1,p2=0.01 q2=0.95
    

    # we assume that tree topology is known in advance

    known_tree = '((a:%f, b:%f):%f,(c:%f,d:%f):%f);'
    tree = LoadTree(treestring='((a,b),(c,d))')

   # generates random parameters at first
    for m in range(3):
        tree1_parameter.append(np.random.exponential(-1.0 / np.log(0.05)))
        tree2_parameter.append(np.random.exponential(-1.0 / np.log(0.05)))
    tree_parameter.append(tree1_parameter)
    tree_parameter.append(tree2_parameter)
    start.append(tree_parameter)

   # iteration between expectation and maximization
    for i in range(100):
        # print tree_parameter

        # expectation with parameters generations from the last iteration
        tree1 = LoadTree(treestring=known_tree % (tree_parameter[0][0], tree_parameter[0][1], tree_parameter[
                         0][2] / 2.0, tree_parameter[0][0], tree_parameter[0][1], tree_parameter[0][2] / 2.0))
        tree2 = LoadTree(treestring=known_tree % (tree_parameter[1][0], tree_parameter[1][1], tree_parameter[
                         1][2] / 2.0, tree_parameter[1][0], tree_parameter[1][1], tree_parameter[1][2] / 2.0))
        result = expectation_singlesite(aln, tree1, tree2)
        # print result

        # optimization with assignments from expectation
        optimization_result = optimization(result, aln, tree1, tree2)
        tree_parameter = optimization_result[0]

        # keep a record of the likelihood, may save and draw plots later for
        # the learning curve
        
        #likelihood.append(optimization_result[1])
        #print optimization_result[1]
        if optimization_result[1] > best_ll:
            best_ll = optimization_result[1]
            best_param = copy.copy(tree_parameter)
    likelihood.append(best_ll)

    print best_param

    #print likelihood  # to see the change of likelihood over 7 iterations

    # to keep q1<q2 for drawing assignment diagram where blue and red should
    # be on one side for random assignments/start parameters
    np.savetxt('random_param\likelihood10_11.txt', likelihood)
    
##    times=range(1,101)
##    plt.plot(times,likelihood)
##    plt.ylabel('log likelihood')
##    plt.xlabel('iteration times')
##    plt.title('random guessing 12')
##    plt.show()
    
    if(best_param[0][1] < best_param[1][1]):
        p1.append(best_param[0][0])
        q1.append(best_param[0][1])
        r1.append(best_param[0][2])
        p2.append(best_param[1][0])
        q2.append(best_param[1][1])
        r2.append(best_param[1][2])
    else:
        p2.append(best_param[0][0])
        q2.append(best_param[0][1])
        r2.append(best_param[0][2])
        p1.append(best_param[1][0])
        q1.append(best_param[1][1])
        r1.append(best_param[1][2])

## #########end of one experiment##########
##
print p1, q1, r1, p2, q2, r2
##np.savetxt('random_param\p1.txt', p1)
##np.savetxt('random_param\q1.txt', q1)
##np.savetxt('random_param\\r1.txt', r1)
##np.savetxt('random_param\p2.txt', p2)
##np.savetxt('random_param\q2.txt', q2)
##np.savetxt('random_param\\r2.txt', r2)
##########end of several experiments#############

