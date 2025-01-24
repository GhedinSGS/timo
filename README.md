<!-- Center the main title -->
<h1 align="center">timo</h1>

<h3 align="center">A single nucleotide variant (SNV) calling tool used to identify viral variants in sequencing data</h3>

---
<h2 align="center"> Table of Contents </h2>

- [Quick Start](#quick-start)
- [Setup and Requirements](#setup-and-requirements)
- [Example input/output data](#example-inputoutput-data)
- [What you can do with the outputs](#what-you-can-do-with-the-outputs)
- [Manuscript information](#manuscript-information)

## Quick Start

1. Clone repo
2. Required inputs
    ```
    --ref (-r): path and filename to reference
    --infile (-i): path and filename to bam file 
    --qual (-q): phred quality cutoff (default: 25) 
    --cutoff (-c): minor frequency variant cutoff (default: 0.01)
    --strain (-T): virus identifier (ex: H1N1, H3N2, SARS2, ZIKA)
    --outputdir (-o): path to where output file should be saved (default: ./)
    --covercutoff (-C): coverage cutoff to filter variants (default: 200X)
    ```
3. Run code with defaults:
   ```
   python3 timo.v4.py --infile bamfiles/filename.sorted.bam --ref reference/ref.name.fa
   ```

4. Outputs: 'snplist' csv file
    -  name: name of input sample
    - segment: name of segment/CHROM in reference
    - ntpos: nucleotide position in segment/CHROM
    - major: major nucleotide at ntpos (most abundant nt at a given position)
    - majorfreq: relative frequency of major variant
    - minor: minor nucleotide (<50%) at ntpos 
    - minorfreq: frequency of minor variant
    - binocheck: binomial check to confirm similar abundance of minor variant in fwd and rev reads
    - A, G, C, T, - : number of reads at the given ntpos with nucleotide or deletion (-)
    - totalcount: total read depth at given position
    - aapos: amino acid position (if aligning to coding sequence only)
    - majoraa: amino acid of major nucleotide sequence (see 'major' column)
    - majorcodon: codon sequence of major. Uppercase nucleotides indicate position of major nt in codon.


## Setup and Requirements

1. **Dependencies and Packages**:
    - python3
    - python libraries: numpy, scipy, pysam, time, oeprator, argparse, sys
  
## Example input/output data
- see example data

## What you can do with the outputs
- SNV variant population dynamics (dN/dS, Shannon entropy, nucleotide diversity, consensus sequence)

---

## Manuscript information

**For more information see: _Optimized Quantification of Intrahost Viral Diversity in SARS-CoV-2 and Influenza Virus Sequence Data_**

**Main contact**: [elodie.ghedin@nih.gov](mailto:main.author@example.com)

**Link to paper**: [mBio link](https://journals.asm.org/doi/10.1128/mbio.01046-23)
