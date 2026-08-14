# REPEATABILITY README

We introduce, estimate, and analyze changes in repeatability — a property of each SNP that quantifies its effect on increasing or decreasing the abundance of direct repeats overlapping the DNA sequence at that position. This property can be biologically relevant for processes such as deletion formation and other DNA alterations.

This pipeline identifies sequence repeats overlapping a specified mtDNA position, runs the search for selected or all mtDNA positions and nucleotides, archives the resulting files, and then extracts summary tables containing the longest perfect repeats. The three files are intended to be run in sequence: function definition → repeat calculation → repeat summarization.

## STEP 1: Define and test the repeat-retrieval function

**File:** `RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence.R`

We define a greedy R function which retrieves a redundant list of repeats overlapping a given position with a given nucleotide. Correctness of this function should be tested.

### Function logic

The function loads `Biostrings` and `stringdist`, then defines `RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence()` with default settings: `max_flank_left = 20`, `max_flank_right = 20`, `min_length = 5`, `max_length = 41`, and `mismatch_fraction = 0.2`. If the nucleotide at the queried position differs from `fixed_nuc`, it is replaced in a copy of the sequence before searching.

The helper function `get_motif_around_position()` extracts a motif containing the selected position with asymmetric left and right flanks, automatically truncated at sequence boundaries. The helper `find_approximate_repeats()` generates every k-mer in the sequence with the same length as the motif, calculates Hamming distances, and retains matches within the allowed mismatch threshold.

The function tests motif lengths from 5 to 41 bases. For each length, it tests every possible division of flanking bases around the position, provided that the left and right flanks do not exceed their configured maxima. For each motif, the permitted number of mismatches is calculated as `floor(mismatch_fraction × motif length)`, normally corresponding to 20% of the motif length.

Matching repeats are annotated with motif sequence, repeat sequence, coordinates, Hamming distance, nucleotide, and `effective.length`, calculated as motif length minus the Hamming distance. All matches are combined, prefixed with the queried position, and sorted by decreasing `effective.length`. If matches exist, they are written as a tab-delimited file named `01KP.<position>.<nucleotide>.txt`; if no matches exist, no file is created and a message is printed.

### Inputs

- A `DNAString` sequence.
- A 1-based genomic position.
- One nucleotide: `A`, `T`, `G`, or `C`.
- An output directory.
- Optional motif, flank, and mismatch settings.

### Outputs

A data frame containing all detected approximate repeats and, when non-empty, a tab-delimited output file with columns including `pos`, `nuc`, `motif.seq`, `motif.length`, motif coordinates, repeat sequence and coordinates, Hamming distance, and `effective.length`.

### Biological adequacy and quality control

This function produces repeats, some of which are not biologically adequate:

- **Self-repeats:** when the motif (a run of DNA overlapping the position with the nucleotide) equals the repeat (run of DNA identical to the motif but not overlapping the position);
- **Biologically not meaningful cases:** when mismatches are located at the start or end (since we allow some degree of degradation);
- **Biologically redundant:** when short repeats are nested completely within longer ones.

All these problems are acceptable and serve as quality control, confirming that the function works correctly from mathematical and algorithmic viewpoints. That is why we keep all such biologically inadequate cases and filter them out later when deriving various metrics of repeatability.

### To do

- Check the logic of the function (compare outputs, verify manually);
- Ask a Python/C++ expert to code it for faster performance.

## STEP 2: Run the function on all mtDNA positions and nucleotides

**File:** `01B.KP.RunTheFunctionOnWholeMtDna.Rmd`

Run this function on all mtDNA positions and all nucleotides, and save outputs.

### Workflow

The R Markdown document loads `knitr`, `dplyr`, `tidyr`, `ggplot2`, `Biostrings`, `stringdist`, and `IRanges`, and sources the repeat-search function. It reads `../data/1_raw/Homo_sapients.mtDNA.fasta` using `readDNAStringSet()` and selects the first sequence as `mtDNA_seq`.

The script first runs the function for positions `8251`, `8472`, `8473`, `12705`, `14798`, and `16223`. At each position, all four nucleotides (`T`, `A`, `G`, and `C`) are tested using the default function settings, and results are written to `/home/popa/Repeatability/data/2_derived/01KP`.

Then it loops through every position in the mtDNA sequence and each nucleotide in `A`, `T`, `G`, and `C`. For each combination, it reports progress and elapsed time, calls the repeat-search function, and creates one output file when repeats are found. For the whole-genome scan, up to four result files can be produced per mtDNA position, with filenames of the form `01KP.<position>.<nucleotide>.txt`. The full scan is computationally expensive and is documented as taking approximately two days.

### Archive and cleanup

The script defines `create_tar_parts()`, which collects the generated text files, estimates each file's archive size, groups files below a target size, and writes numbered gzip-compressed tar archives. The example archives `/home/popa/Repeatability/data/2_derived/01KP` into `/home/popa/Repeatability/data/2_derived/01KP.ARCHIVE` using `target_mib = 300`.

Files are named `repeatability_data_001.tar.gz`, `repeatability_data_002.tar.gz`, and so on. The original text files are not deleted or modified by this block. The README instructs the user to manually delete the individual `.txt` files from `/home/popa/Repeatability/data/2_derived/01KP` before pushing the repository. This keeps the large individual result files out of GitHub while retaining the compressed archives.

## STEP 3: Derive repeatability metrics and perform descriptive statistics

**File:** `01C.KP.TryRepeatabilityMetrics.Rmd`

Derive repeatability metrics and perform descriptive statistics of repeatability for different positions, sites, and model mutations (mutations known to be associated with specific phenotypes).

### Extract archived results

The script defines `extract_tar_archives()`, which finds archives matching `repeatability_data_<number>.tar.gz`, sorts them, and extracts them into an output directory. The destination must be empty; archives are retained by default because `delete_archives = FALSE`. It extracts archives from `/home/popa/Repeatability/data/2_derived/01KP.ARCHIVE` into `/home/popa/Repeatability/data/2_derived/01KP.RESTORED`.

### Longest perfect repeats across the whole mtDNA

The script reads every `.txt` file in the restored directory using the expected column types. It retains records where the repeat is different from the motif, the repeat is a perfect match (`repeat.hamming.distance == 0`), and the queried position lies strictly inside the motif rather than at its boundary.

For each input file, the script keeps all rows having the maximum `effective.length`. It combines results from all files, removes duplicate rows, sorts by position, and writes them to `/home/popa/Repeatability/data/2_derived/01C.TheLongestPerfecInWholeMtDna.txt`. This produces a summary table containing the longest direct perfect repeats across the complete mtDNA sequence.

### Longest perfect repeats within the major arc

The analysis is repeated but restricts both the motif and repeat to the major-arc interval from position `5798` through `16568`. The same perfect-match, non-identical-repeat, and strictly-internal-position filters are applied. The output is written to `/home/popa/Repeatability/data/2_derived/01C.TheLongestPerfecWithinMajorArc.txt`, containing the longest qualifying perfect repeats located entirely within the major arc.

## STEP 4: Consider extensions

- Visualisation of potential deletions (e.g., using R Shiny, circus plots);
- Analysis of ssDNA and ssRNA viral genomes;
- When fast code becomes available, analysis of interactions (epistasis) between different variants (probably focusing on closely located synonymous mutations separately, but never together);
- Analyse other mammalian/vertebrate mtDNAs.
