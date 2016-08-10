# Pseudocode for testing phylogenetic methods
# based on the 4 taxon case to start with...

########################### Simulation ##############################

# 1. Define a parameter space where we know traditional methods don't
#    always get the right answer. This we can get from Phillipe's ms
#    This should include: p, q, r, w, t, and alignment_length
#    

# 2. Write a generator to produce a simulated alignment given 
#    p, q, r, w, t, alignment_length, an unrooted tree of 4 taxa.
#	 We could also consider number_of_taxa as a parameter, but it is
#	 simpler to just start with the 4 taxon case. The existing 
# 	 methods are well characterised here, and it's good to start 
#	 simple. 

############################# Testing ###############################

# I'm assuming we have at least 4 methods:
#		1. Maximum Likelihood (ML)
#		2. Neighbour Joining (NJ)
#		3. Maximum Parsimony (MP)
#		4. A trained Artificial Intelligence (AI) 

# To test the methods, we use BL50 following Phillipe et al:

# 1. Pick some number of points in parameter space using:
#    p, q, w, t, alignment_length

# 2. For each point in parameter space: 

#	 	Generate N (at least 100) alignments with different r, from
#    	1/N to 1 in increments of 1/N

#		For each method:
#			Calculate BL50 (the minimum r for which the method gets
#			the correct tree in at least 50% of cases). Use linear
#			regression.

# Note: there are good existing implementations of ML, NJ, and MP.
# In Python, we can use BioPython for NJ and MP, and PyCogent for
# ML. 

# Now, for any given point in parmater space, we simply prefer the
# method with the best (i.e. lowest) BL50. 

# The question here is whether it's possible to train an AI to do 
# better than the current approaches.