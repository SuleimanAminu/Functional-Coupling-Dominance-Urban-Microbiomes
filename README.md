# Functional Coupling, Dominance Hierarchies, and Resilience Landscapes in Urban Microbiomes

This repository contains the computational analyses used to examine ecological
principles of functional coupling, dominance hierarchies, 
and resilience in complex microbial systems.

Urban microbiomes are used here as a model system to examine how
antimicrobial resistance (ARG) and virulence factor (VF) mechanisms are
organized, integrated, and redistributed across environments. Rather
than focusing solely on taxonomic turnover, the analysis evaluates how
functional structure reveals generalizable ecological rules governing
stability, override dynamics, and hierarchical dominance.

All analyses operate at the level of harmonized functional mechanisms
derived from genome-resolved metagenomic data.

Environments analyzed include:

-   Ambulance\
-   Hospital interior\
-   Hospital sewage\
-   Public transport

------------------------------------------------------------------------

## Conceptual Framework

This workflow evaluates three central ecological dimensions:

### 1. Functional Coupling Theory

Resistance and virulence functions 
are conceptualized as interacting ecological modules.\
To quantify their integration, the framework computes:

-   Shannon entropy of ARG mechanisms\
-   Shannon entropy of VF mechanisms\
-   Total functional entropy\
-   Functional Coupling Ratio (FCR)\
-   Integration Index

These metrics describe how tightly resistance and virulence functions
are balanced and co-structured within environments.

### 2. Dominance Hierarchies

Environment-level functional centroids are compared using predictive
separability and perturbation modeling.\
Dominance is operationalized as the capacity of one functional 
configuration to override another under controlled compositional perturbation

### 3. Resilience Landscapes

Synthetic donor--recipient blending is used to construct
functional stability gradients.\
For each synthetic state:

-   Functional entropy is recalculated\
-   Prediction confidence is evaluated\
-   Override thresholds are identified

These dynamics define resilience boundaries and quantify functional 
fragility within and across environments.
------------------------------------------------------------------------

## Analytical Workflow

| Analytical Stage | Data Inputs | Ecological Outputs |
|------------------|------------|--------------------|
| Genome Processing | ARG annotations (CARD), VF annotations (VFDB), MAG taxonomy | Harmonized mechanism-level functional table |
| Batch Correction & CLR | Mechanism TPM matrix | CLR matrix, batch-corrected matrix |
| Ecological Modeling | CLR matrix + metadata | Diversity metrics, PERMANOVA results, differential mechanisms |
| Functional Integration | Mechanism abundances | ARG entropy, VF entropy, FCR, Integration Index |
| Machine Learning | CLR features + integration metrics | Classification performance, feature importance rankings |
| Synthetic Perturbation | Environment-level functional centroids | Override thresholds, dominance matrix, resilience gradients |

------------------------------------------------------------------------

## Repository Structure

    01_genome_processing/
    02_batch_correction/
    03_ecology_analysis/
    04_machine_learning/
    05_synthetic_ecology/

Scripts follow the analytical progression described above.

------------------------------------------------------------------------

## Requirements

### Python ≥ 3.9

-   pandas\
-   numpy\
-   scipy\
-   scikit-learn\
-   matplotlib\
-   seaborn

### R ≥ 4.3

-   tidyverse\
-   vegan\
-   limma\
-   compositions\
-   ComplexHeatmap\
-   FSA\
-   ggpubr

------------------------------------------------------------------------

## Reproducibility

All analyses use fixed random seeds.\
Batch correction preserves biological signal.\
Model validation prevents information leakage between training and
holdout data.

------------------------------------------------------------------------

## License

MIT License — free to use, adapt, and cite with attribution. 

This repository forms part of a manuscript currently under consideration for publication.
We respectfully request that the code not be redistributed, repackaged, or used in derivative publications until the peer-review process is complete.

If you intend to build upon this work prior to publication, please contact suleiman.aminu@um6p.ma.

---

 
