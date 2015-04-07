ParsePhyloBayes
========================================================

Overview
---------

This program is developed for the manuscript ``Relaxing the Molecular Clock to Different Degrees for Different Substitution Types'' in Molecular and Evolutionary Biology. 
It reads a sample of substitution histories (mappings) from PhyloBayes. 
It then computes the substitution lengths for context-independent or context-dependent substitutions for branches and the variance-covariance matrices of substitution lengths between branches. 
It produces the outputs in the format of the inputs for Multidivtime so that each 
substitution type can be treated as a single gene (locus) in Multidivtime.


Getting Started
---------------

### PhyloBayes 

1. Prepare the nucleotide sequence alignment in the PHYLIP format. Notice that the alignment should include the sequences from the outgroup species.
   ```
   <number_of_taxa> <number_of_sites>
   taxon1 sequence1 ...
   taxon2 sequence2 ...
   ...
   ```

2. Prepare the unrooted phylogenetic tree topology in NEWICK format. Save it as a .tre file.

3. Download and install the modified version of PhyloBayes MPI 1.5a from
   [http://megasun.bch.umontreal.ca/People/lartillot/www/index.htm](http://megasun.bch.umontreal.ca/People/lartillot/www/index.htm).
   
4. Run PhyloBayes analyses at least twice using the following command.
   ```
   mpirun -np <n> pb_mpi  -gtr -cat -dgam 4 -d <dataset> -T <treetopology> <chainname>
   ```
   More commands for PhyloBayes are available at [http://megasun.bch.umontreal.ca/People/lartillot/www/pb_mpiManual1.5.pdf](http://megasun.bch.umontreal.ca/People/lartillot/www/pb_mpiManual1.5.pdf).
   
5. Assess the MCMC convergence by comparing two independent chains using **tracecomp** from PhyloBayes.
   
6. Generate substitution histories (mappings) with the following command.
   ```
   mpirun -np <n> readpb_mpi -map <chainname>
   ```
   This will produce a series of .map files, each is the substitution histories of a single site.
   
   
### ParsePhyloBayes

1. Download **parsePhyloBayes.jar** and copy it to the directory that contains the .map files generated from PhyloBayes.

2. Prepare a text document that records the outgroup information. The outgroup information should be in the following format.
   ```
   <number_of_outgroup_species>
   <outgroup_species1>
   <outgroup_species2>
   ...
   ```
3. Run parsePhyloBayes.jar using the following command.
```
java -jar parsePhyloBayes.jar <number_of_sites> <number_of_samples_per_site> <chainname> <context_dependent_or_independent> <outgroup_file>
```
  * The `<number_of_sites>` argument is the length of the sequence alignment.
  * The `<number_of_samples_per_site>` argument is the number of substitution history samples generated per site. This number can be checked from the .map files.
  * The `<chainname>` argument is the chain name specified in PhyloBayes analyses.
  * The `<context_independent_or_dependent>` argument can be either 1 or 2, where 1 is the context-independent case, and 2 is the context-dependent case, respectively. The definition of substitution types for context-independent or context-dependent substitutions can be found in the manuscript Table 1 and Table 2.

4. The outputs of parsePhyloBayes.jar include
   * o.estb.type 
   * o.estb.group 
   * substitutionLength.txt
   
   o.estb.type and o.estb.group can be used as the inputs to Multidivtime. The analyses in the manuscript used o.estb.group as inputs to Multidivtime.
   
   o.estb.type for single nucleotide substitutions are for 12 different substitution types.
   o.estb.group for single nucleotide substitutions are for 6 different substitution types (combine strand symmetric substitution types).
   
   o.estb.type for triple site nucleotide substitutions are for 4 different substitution types (non-CpG/CpG transversion/transition).
   o.estb.group for triplet site nucleotide substitutions are for 9 different substitution types (defined in manuscript).
   
   substitutionLength.txt file stores the estimated substitution lengths for each strand symmetric substitution type. This makes it easier to report results.

### Multidivtime

1. Prepare a .tree file that contains the rooted tree topology including the outgroup species.

2. Download and compile the Multidivtime software from [http://statgen.ncsu.edu/thorne/multidivtime.html](http://statgen.ncsu.edu/thorne/multidivtime.html). More details are 
in the website given above.

3. Edit the `multicntrl.dat` file. One example multicntrl.dat file is as follows.
```
/* the following lines are all needed in multicntrl.dat ...
   do not add or delete lines but change entry on left of each
   line as you see fit ...  */
9species.tree
9 ... number of genes ... FOLLOWING LINES CONTAIN ONLY NAMES OF DATA FILES
o.estb.group0
o.estb.group1
o.estb.group2
o.estb.group3
o.estb.group4
o.estb.group5
o.estb.group6
o.estb.group7
o.estb.group8
10000 ... numsamps: How many times should the Markov chain be sampled?
1000 ... sampfreq: How many cycles between samples of the Markov chain?
10000000 ... burnin: How many cycles before the first sample of Markov chain?
44.2 ... rttm: a priori expected number of time units between tip and root
0.1 ... rttmsd: standard deviation of prior for time between tip and root
0.0012 ... rtrate: mean of prior distribution for rate at root node
0.0012 ... rtratesd: standard deviation of prior for rate at root node
0.0 ... brownmean: mean of prior for brownian motion constant "nu"
0.0 ... brownsd: std. deviation of prior for brownian motion constant "nu"
/* the following lines are all needed (i.e., do not delete them) but you may 
   not want to alter entries unless you are familiar with the computer code */
1.0  ... minab: parameter for beta prior on proportional node depth
0.1 ... newk: parameter in Markov chain proposal step
0.5 ... othk: parameter in Markov chain proposal step
0.5 ... thek: parameter in Markov chain proposal step
100.0 ... bigtime: number higher than time units between tip and root could be in your wildest imagination
/* the program will expect the entry below to be the number of constraints
   and then the specified number of constraints should follow on
   subsequent lines */
0 ... number of constraints on node times
0 ... number of tips which are not collected at time 0
0  ... nodata: 1 means approximate prior, 0 means approximate posterior
0 ...commonbrown: 1 if all genes have same tendency to change rate, 0 
```

    Note: the `commonbrown` parameter should be set to 0 so that different genes (substitution types) are allowed to change rates independently.
    
4. Copy ``o.estb.group*`` files to the folder with compiled multidivtime. Run the Multidivtime analysis with the following commend:
```
multidivtime <name> > & out.<name> &
```
