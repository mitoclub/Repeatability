# Repeatability

We introduce, estimate, and analyze changes in repeatability — a property of each SNP that quantifies its effect on increasing or decreasing the abundance of direct repeats overlapping the DNA sequence at that position. This property can be biologically relevant for processes such as deletion formation and other DNA alterations.

# STEP 1:

We define a greedy R function which retrieves a redundant list of repeats overlapping a given position with a given nucleotide. 
Correctness of this function should be tested (`RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence.R`).

This function produces repeats, some of which are not biologically adequate:

- Self-repeats: for example, when the motif (motif here is a run of DNA overlapping the position with the nucleotide) equals the repeat (run of DNA identical to the motif but not overlapping the position);
- Biologically not meaningful cases: when mismatches are located at the start or end (since we allow some degree of degradation);
- Biologically redundant: when short repeats are nested completely within longer ones.

All these problems are acceptable and serve as quality control, confirming that the function works correctly from mathematical and algorithmic viewpoints. That is why we keep all such biologically inadequate cases and filter them out later when deriving various metrics of repeatability.

To do:

- Check the logic of the function (compare outputs, verify manually);
- Ask a python/C++ expert to code it for faster performance.

# STEP 2:

Run this function on all mtDNA positions and all nucleotides, and save outputs.

# STEP 3:

Derive repeatability metrics and perform descriptive statistics of repeatability for different positions, sites, and model mutations (mutations known to be associated with specific phenotypes).

# STEP 4:

Consider extensions:

- Visualisation of potential deletions (e.g., using R Shiny, circus plots);
- Analysis of ssDNA and ssRNA viral genomes;
- When fast code becomes available, analysis of interactions (epistasis) between different variants (probably focusing on closely located synonymous mutations separately, but never together).
- Analyse other mammalian/vertebrate mtDNAs

