# GoPrime
GoPrime was designed by [Emma Howson](https://www.pirbright.ac.uk/users/dr-emma-howson) (Pirbright Institute) and [Richard Orton](https://www.gla.ac.uk/researchinstitutes/iii/staff/richardorton/) (University of Glasgow), in collaboration with [Don King](https://www.pirbright.ac.uk/users/dr-don-king) (Pirbright Institute) and [Veronica Fowler](https://www.researchgate.net/profile/Veronica_Fowler).  

**Citation:** a manuscript describing GoPrime has recently been submitted:\
GoPrime: development of an in silico framework to predict the performance of real-time PCR primers and probes using foot-and-mouth disease virus as a model. Emma L A Howson, Richard J Orton, Valerie Mioul, Sarah Cleaveland, Tiziana Lembo, Donald P King and Veronica L Fowler

**Abstract**
Real-time PCR (rPCR) is a widely accepted diagnostic tool for the detection and quantification of nucleic acid targets. In order for these assays to achieve high sensitivity and specificity, primer and probe-template complementarity is essential; however mismatches are often unavoidable and can result in false-negative results and errors in quantifying target sequences. Primer and probe sequences therefore require continual evaluation to ensure they remain fit for purpose. This paper describes the development of a linear model and associated computational tool (GoPrime) designed to predict the performance of rPCR primers and probes across multiple sequence data. Empirical data were generated using DNA oligonucleotides (n = 90) that systematically introduced variation in the primer and probe target regions of an assay used to detect foot-and-mouth disease virus (FMDV). These assays revealed consistent impacts of patterns of substitutions in primer and probe-sites on rPCR cycle threshold (CT) and limit of detection (LOD). These data were used to populate GoPrime, which was subsequently used to predict rPCR results for DNA templates (n=7) representing the natural sequence variability within FMDV. GoPrime was also applicable to other areas of the FMDV genome, with predictions for the likely targets of a FMDV-typing assay consistent with published experimental data. Although further work is required to improve these tools and assess the broader impact of mismatches for other assays, these data support the use of mathematical models for rapidly predicting the performance of rPCR primers and probes in silico.

**Usage:**\
GoPrime is written in the Java programming language, the Java jar file can be run with the following command:

```
java -jar GoPrime.jar primer_sequences.fasta target_sequences.fasta
```

**primer_sequences.fasta**\
This file should contain the 2 primer sequences (forward and reverse) and single probe sequence in FASTA format, 5’-3’ direction, and in the below order:
```
>FwdPrimer
ACTGGGTTTTACAAACCTGTGA
>Probe
TCCTTTGCACGCCGTGGGAC
>ReversePrimer
GCGAGTCCTGCCACGGA
```
The order is essential as GoPrime uses the order of the sequences rather than their names to determine which is the probe, forward and reverse primers. GoPrime can handle ambiguity codes in both primers and probes.

**target_sequences.fasta**\
This file contains the target sequences, such as full or partial viral genomes, that you went to evaluate your primer and probe set against. This should be in FASTA format, and in 5’-3’ direction. GoPrime can handle ambiguity codes in the target sequences.

GoPrime operates by taking each target sequence in turn. It then searches the entire sequence for potential primer/probe binding sites that meet the minimum priming criteria, see Howson *et al.* (submitted) for full details but briefly:
```
Primers
Maximum of 2 mismatches between the primers and template at the 3’ end (last 4 nucleotides)
Minimum of 82.05% match between primers and template

Probes
Minimum of 85% match between the probe and the template
```
After finding candidate priming sites for the probe and forward and reverse primers individually, these are then evaluated to see if any combinations of these can form a successful set that are the correct orientation (with respective to each other) and that there are no overlaps. If a successful set is found, GoPrime will then evaluate if it will prime successfully (based on the above rules), and if success GoPrime will calculate an expected deltaCT and deltaLOD score based on the penalties calculated from the linear model analysis, see Howson *et al.* (submitted) for full details of the mismatches penalties.

GoPrime will output two text files:

**target_sequences.fasta_cts.txt**
This file contains the deltaCT and deltaLOD scores for each target sequence. Fields are:
```
1)	Target sequence name
2)	deltaCT
3)	deltaLOD
4)	FwdPrimer – deltaCT and start/end co-ordinates of primer binding site on the target sequence
5)	Probe– deltaCT and start/end co-ordinate of probe binding site on the target sequence
6)	RevPrimer – deltaCT and start/end co-ordinates of primer binding site on the target sequence
7)	Length of PCR product
```

**target_sequences.fasta_out.txt**
This file contains the output messages from GoPrime as it processes each target sequence in turn. It will output the details of all candidate binding sites for primers and probe individually, as well as evaluation details of each potential primer/probe set. An example output is below, which is evaluating a single target sequence (TestSeq):
```
Evaluating seq 1 >TestSeq
23 position is a candidate for 5'-Fwd-3'-fwd 95.45% 97.44% [%Match %MatchPair], 1 1 [TotMis Tot1-4]
72 position is a candidate for 5'-Probe-3'-probe 100% 100% [%Match %MatchPair], 0 0 [TotMis Tot1-4]
92 position is a candidate for 3'-Rev-5'-rev 100% 100% [%Match %MatchPair], 0 0 [TotMis Tot1-4]
1-[1]-1 Fwd-[Probe]-Rev individual candidate primer/probe positions found in expected orientation
0-[1]-0 Fwd-[Probe]-Rev individual candidate primer/probe positions found in opposite orientation
Set-1 RT-PCR Success Fwd=2/23 Probe=53/72 Rev=108/92
deltaCT=3.4 [fCT=3.4 pCT=0 rCT=0]
```
In the above output GoPrime identifies a single potential binding site for the forward (postion 23), probe (position 72), and reverse (position 92) primers. In each case, the %Match between the primer/probe and template is reported (and for primers a %MatchPair when considering the combined fwd+rev primer length), and the totoal number of mismatches across the whole primer (TotMis) vs the total in the 3’ end (last 4 nucleotides) of the primer (Tot1-4). A single potential fwd/probe/rev set is found (1-[1]-1) in the expected orientation (fwd going fwd in direction, rev going in rev direction), and this has a predicted deltaCT of 3.4. This is a result of mismatches between the forward primer and the template, as the probe and reverse primer are perfect matches.

GoPrime is idealy suited to evaluate the sensitivity and specificty of genotype specific real time PCR primers and probes for viral data sets.

# GoPrimeTree
GoPrimeTree is a simple program to visualise the deltaCT results generated by GoPrime (in the target_sequences.fasta_cts.txt output file) on a phylogenetic tree. GoPrimeTree expects a phylogenetic tree in nexus (.nex) format, that contains all the target sequences with the same sequence names being used; we recommend generating this in [https://github.com/rambaut/figtree](FigTree), which saves trees as nexus format by default. GoPrimeTree will then add [https://github.com/rambaut/figtree](FigTree) displayable colour coding to the tree nexus file, so the result deltaCT results can be visualised in [https://github.com/rambaut/figtree](FigTree). GoPrimeTree using different shades of the colours green, orange and red are  to represent target sequences with ΔCT‘s between 0-10, 10-20, and 20-30, respectively, whilst black is used for sequences that fail to amplify..

**Usage:**\
GoPrimeTree is written in the Java programming language, the Java jar file can be run with the following command:

```
java -jar GoPrimeTree.jar target_sequences.fasta_cts.txt  tree.nex
```
GoPrimeTree will then output a new nexus tree file called: tree_col.nex, which can be opened in [https://github.com/rambaut/figtree](FigTree)

# Example data
We provide two example data sets, both relating to Foot-and-Mouth Disease Virus (FMDV)




