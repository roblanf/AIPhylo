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

The tree seems to be arbitrarily midpoint rooted. I think that's just a PyCogent thing for viewing