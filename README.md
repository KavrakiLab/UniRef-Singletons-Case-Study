# UniRef Singletons Arising from Protein Sequence Clustering Heuristics: A Case Study

# Abstract

With protein databases now containing billions of protein sequences, protein sequence clustering helps to tame the protein sequence database deluge by attempting to collate proteins that are functionally related. Computationally, this facilitates sequenced based comparisons to only the cluster representative, and not all of the cluster members, such as in the case of UniRef. 
Clustered databases underpin a range of downstream applications, including protein structure and function prediction. To cluster sequences quickly, recent algorithms utilize several heuristics which can introduce functional inconsistencies in clusters. This can occur when sequences are either incorrectly excluded from protein clusters containing functionally related proteins or, conversely, are erroneously grouped with proteins of unrelated function. 
By examining modern clustering algorithms, we show that the high number ($\sim$60\%) of single-member clusters has a considerable impact on the quality of functional annotation metrics associated with clusters, the so-called clustering statistics. We identify three cases that illustrate the broader impacts of modern fast clustering algorithms to downstream tasks such as protein structure prediction.

# Requirements
External protein clustering tools are needed to replicate benchmarking results. In particular, [MMSeqs2 Software Suite](https://github.com/soedinglab/MMseqs2) and [CD-HIT](github.com/weizhongli/cdhit). Please follow the installation steps provided by their respective software vendors to install their software. Commands do not need to be added to the path for this code to run.

# Installation

Created from a conda environment.
```bash
conda env create -f environment.yml
conda activate UniRef-Singletons-Case-Study
```

Running Code:

```bash
python main.py {dataset_year} {path_to_mmseqs_executable} {path_to_cd_hit_executable} {gene ontology obo file}
```
