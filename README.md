# AIPhylo

This is the beginning of an attempt to train an AI to infer phylogenies. 

# Dependencies
PyCogent

# Simulating alignments

We follow the scheme set out in Philippe et al 2005: http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-5-50

```{python}
from simulator import simulate_alignment

# simulate a 100bp alignment
aln, tree = simulate_alignment(0.1, 0.5, 0.1, 100, 1.0, 'random')
print aln
print tree.asciiArt()
```

Should give you something like this:

```
>a
TATCATTCTATACGTTGACCGCTTTTAGCGCAAGTTAAGCGGTTCTGTAATATCTGGTCTTGTAGTACTGTTACGAGGAGCTAGCTACCCTTTAGGGCTA
>c
GGCCATGGTAAAAAAGATGTTGCCAACCCCGGCGCTACCATGATAAATGGGTTGAGACTCGTCATTTGCCCGGATCAGGTTGTCGTCGCTGCCGACGACG
>b
GGCCATGGTAACAAAGATGTGGCCACGCCCGGCGCTACCATGATAAATGCTATGAGATGTCAGTCTCGAGCGGAAAGCGAGTCGCTACTTGCCCAGGCTT
>d
TATCATTCTATACGTTGACTGCTTTTCGCGCACGTTAAAGGGTTATGTAATCTCCGGTCTGGTGGGATTCTGAGGTCGGCAAGCCAAGCTTTTAGCGTCG


                    /-a
          /edge.0--|
         |          \-d
-root----|
         |          /-b
          \edge.1--|
                    \-c



```

The tree seems to be arbitrarily midpoint rooted. I think that's just a PyCogent thing for viewing.

# Is the ML tree the true tree?

We can answer this with the ml4 function, which finds the ML tree by getting the likelihood (under the JC69 model) of all 3 unrooted 4-taxon trees, and returns ```True``` if the true tree is the ML tree, an d ```False``` otherwise:

```
from simulator import *
from ml import *

aln, tree = simulate_alignment(0.675, 0.15, 0.16, 10000, 0.9, 'random')
ml4(aln, tree) 
```

# A quick test
In Phillipe et al figure 3, the suggestion is that when p = 0.675, q = 0.15, and t = 0.9, the BL50 is ~0.20. We can roughly test this by simulating a set of datasets with r = 0.175 and seeing if we get roughly 50% of datatsets where the ML tree is the true tree, and 50% where it isn't.

```
input = [simulate_alignment(0.675, 0.15, 0.2, 10000, 0.9, 'random') for _ in range(20)]
results = [ml4(i[0], i[1]) for i in input]

# percent True
float(sum(results))/float(len(results))
```

When I do this, I get 0.65, which is close enough to suggest our code works as expected.