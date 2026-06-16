# Beta-Burst-Detection-Using-Nonlinear-Features

MATLAB code for detecting and classifying EEG beta bursts using nonlinear dynamic features and machine learning.

## Overview

This repository contains the MATLAB implementation developed during my Master's thesis and the associated published study:

Designing a model to detect beta burst in EEG using nonlinear dynamic features based on machine learning

The project investigates whether nonlinear dynamic features can characterize beta bursts across different motor tasks and burst durations, and whether these features can be used for automated beta-burst classification.

## Authors

- Armin Hakkak Moghadam Torbati*†
- Narges Davoudi†
- Giuseppe Longo

\* Corresponding author

† Co-first authors

## Publication

Hakkak Moghadam Torbati A†, Davoudi N†, Longo G.

*Designing a model to detect beta burst in EEG using nonlinear dynamic features based on machine learning.*

Biomedical Signal Processing and Control, 2026.

DOI: [(https://doi.org/10.1016/j.bspc.2026.110115)]

> Note: The final published PDF is not included in this repository due to publisher copyright restrictions.

## Dataset

- Participants: 26 healthy adults
- EEG system: 64-channel EEG
- Sampling rate: 1000 Hz
- Tasks:
  - Right-hand isometric contraction (10% MVC)
  - Right-hand dynamic movement
  - Bimanual dual task

## Processing Pipeline

1. EEG preprocessing
2. Beta-band filtering (13–30 Hz)
3. Hilbert envelope extraction
4. Beta-burst detection using a 75th percentile threshold
5. Burst duration categorization:
   - Short (< 100 ms)
   - Medium (100–149 ms)
   - Long (≥ 150 ms)
6. Nonlinear feature extraction
7. Machine learning classification

## Extracted Features

- Fractal Dimension (FD)
- Sample Entropy (SE)
- Wavelet Entropy (WE)
- Nonlinear Energy Operator (NEO)

## Machine Learning Models

The following classifiers were evaluated:

- Decision Tree (DT)
- Logistic Regression (LR)
- Linear Discriminant Analysis (LDA)
- Quadratic Discriminant Analysis (QDA)
- Naive Bayes (NB)
- k-Nearest Neighbors (KNN)
- Support Vector Machine (SVM)
- Random Forest (RF)
- Neural Network (NN)

## Main Findings

- Nonlinear features varied significantly across beta-burst durations and motor tasks.
- Fractal Dimension, Sample Entropy, and Wavelet Entropy showed the highest discriminative power.
- The proposed framework achieved up to 91.1% validation accuracy and 85.7% independent test accuracy.
  

  ## Requirements

- MATLAB R2024b or later
- FieldTrip toolbox

FieldTrip: https://www.fieldtriptoolbox.org/

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Torbati2026BetaBurst,
  author = {Hakkak Moghadam Torbati, Armin and Davoudi, Narges and Longo, Giuseppe},
  title = {Designing a model to detect beta burst in EEG using nonlinear dynamic features based on machine learning},
  journal = {Biomedical Signal Processing and Control},
  year = {2026},
  doi = {10.1016/j.bspc.2026.110115}
}
```
## Contact

Narges Davoudi

MSc in Data Science

University of Naples Federico II

Email: n.davoudi@studenti.unina.it

Alternative email: narges.davoudi@yahoo.com


## Repository Structure

├── code/
│   └── MATLAB source files

└── README.md
