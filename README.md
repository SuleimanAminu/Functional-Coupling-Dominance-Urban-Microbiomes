# Coordinated Resistance–Virulence Organization Across Urban Microbiomes
This repository contains the computational workflows, ecological analyses, and machine-learning frameworks used in the study:
## "Coordinated resistance–virulence organization across urban microbiomes reveals context-dependent ecological structuring"
### currently under review at npj Biofilms and Microbiomes.

The project investigates how antimicrobial resistance (ARG) and virulence factor (VF) mechanisms are coordinated within genome-resolved urban microbiomes across built environments, and how these coordinated functional architectures relate to ecological differentiation, perturbation dynamics, and predictive environmental organization.

Rather than focusing exclusively on taxonomic turnover or independent gene inventories, the framework evaluates how higher-order functional organization reflects broader ecological processes including:
-   environmental filtering
-   environmental filtering
-   persistence-associated ecological tendencies
-   ecological interaction and dispersal
-   context-dependent multifunctionality
-   perturbation-transition dynamics
-   instability landscapes across urban systems
All analyses operate on harmonized mechanism-level functional profiles derived from genome-resolved metagenomic data.

## Built Environments Analyzed
The study evaluates microbial communities from four major built-environment systems:
- Ambulance interiors
- Hospital interior environments
- Hospital sewage systems
- Public transport infrastructure

These environments represent contrasting ecological regimes differing in:
- antimicrobial exposure
- environmental filtering
- ecological connectivity
- disturbance frequency
- microbial exchange potential
- infrastructural heterogeneity
------------------------------------------------------------------------

## Conceptual Framework

TThe analytical framework integrates genome-resolved metagenomics, ecological theory, compositional data analysis, machine learning, and synthetic perturbation modeling to evaluate coordinated resistance–virulence organization across urban microbiomes.

### 1. Coordinated Resistance–Virulence Organization
Resistance and virulence functions are interpreted as coordinated ecological components embedded within shared genomic and microbial community contexts.
Resistance and virulence functions.
To quantify functional organization, the framework computes:

-   Shannon entropy of ARG mechanisms
-   Shannon entropy of VF mechanisms
-   Total functional entropy
-   Functional Coupling Ratio (FCR)
-   Integration Index

These metrics characterize the relative balance and organization of resistance-associated and virulence-associated ecological processes within environments.
Importantly, the framework interprets "integration" in an ecological and statistical sense describing coordinated organization and recurrent co-occurrence rather than direct mechanistic coupling between individual genes.

### 2. Environment-Associated Ecological Organization

Functional mechanisms are abstracted into higher-order ecological tendencies including:
- persistence-associated tendencies
- interaction-associated tendencies
- context-dependent / multifunctional tendencies
These ecological tendencies are interpreted as overlapping and context-dependent rather than mutually exclusive biological states.
The framework evaluates how different built environments exhibit distinct but partially overlapping ecological configurations shaped by:
- environmental filtering
- microbial dispersal
- ecological connectivity
- anthropogenic pressures
- stress adaptation

### 3. Predictive Ecological Organization
Machine-learning analyses evaluate whether coordinated functional architectures encode reproducible environment-associated ecological structure.

Classification analyses incorporate:
- mechanism-level CLR-transformed features
- resistance–virulence integration metrics
- entropy-derived ecological variables

The framework compares:
- species-level ecological organization
- mechanism-level ecological organization
- combined species–function ecological structure
to evaluate the contribution of shared taxonomic and functional organization.

### 4. Synthetic Perturbation and Transition-State Ecology
Synthetic donor–recipient mixing experiments are used to evaluate perturbation behavior of coordinated resistance–virulence architectures.
Environment-level functional centroids are interpolated across donor fractions (α) to generate synthetic ecological states.
For each perturbation trajectory:
- prediction confidence is evaluated
- entropy landscapes are calculated
- first-transition thresholds are identified
- intermediate ecological states are quantified
The perturbation framework emphasizes:
- non-linear ecological transitions
- instability landscapes
- transitional ecological configurations
- context-dependent perturbation dynamics
rather than strict deterministic dominance hierarchies.
------------------------------------------------------------------------

## Analytical Workflow

| Analytical Stage |  Inputs |  Outputs |
|------------------|------------|--------------------|
| Genome Reconstruction | Shotgun metagenomes | Dereplicated MAG catalog |
| Functional Annotation | CARD + VFDB annotations | Harmonized ARG/VF profiles |
| Mechanism Aggregation | Gene-level TPM abundances | Mechanism-level abundance matrix |
| Batch Correction & CLR | Mechanism abundance matrix | Batch-corrected CLR matrix |
| Ecological Modeling | CLR matrix + metadata | Diversity metrics, PERMANOVA, CAP |
| Differential Enrichment | Mechanism-level profiles | Environment-associated enrichment patterns |
| Functional Integration | Mechanism abundances | Entropy metrics, FCR, Integration Index |
| Machine Learning | CLR features + integration metrics | Classification models and feature importance |
| Synthetic Perturbation | Environment centroids | Transition trajectories and entropy landscapes |

------------------------------------------------------------------------

## Repository Structure

    01_genome_processing
    02_batch_correction
    03_ecology_analysis
    04_machine_learning
    05_synthetic_ecology

Scripts follow the analytical progression described above.

------------------------------------------------------------------------

## Requirements

### Python ≥ 3.9

-   pandas
-   numpy
-   scipy
-   scikit-learn
-   matplotlib
-   seaborn

### R ≥ 4.3

-   tidyverse
-   vegan
-   limma
-   compositions
-   ComplexHeatmap
-   FSA
-   ggpubr

------------------------------------------------------------------------

## Reproducibility
All analyses use fixed random seeds.
Batch correction preserves environment-associated ecological signal.
Nested cross-validation prevents information leakage.
Supplementary robustness analyses evaluate taxonomic contribution to predictive performance.
Synthetic perturbation analyses explicitly model intermediate ecological states and non-linear transition dynamics.

------------------------------------------------------------------------

## Manuscript Status
This repository accompanies the manuscript:
### "Coordinated resistance–virulence organization across urban microbiomes reveals context-dependent ecological structuring"
currently 
### under peer review at npj Biofilms and Microbiomes.

------------------------------------------------------------------------
# Citation and Usage
If you use this repository, analytical framework, or associated workflows, please cite the corresponding manuscript once published.
Until the peer-review process is complete, we respectfully request that the repository not be redistributed, repackaged, or used in derivative publications without prior communication with the authors.
For collaborations, questions, or reuse requests, please contact:

Suleiman Aminu
(suleiman.aminu@um6p.ma)
or 
Rachid Daoud
(rachid.daoud@um6p.ma)
------------------------------------------------------------------------
## License

MIT License — free to use, adapt, and cite with attribution. 

---

 
