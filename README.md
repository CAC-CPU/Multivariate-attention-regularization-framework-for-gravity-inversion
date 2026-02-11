# Multivariate-attention-regularization-framework-for-gravity-inversion
# Overview
This repository provides a MATLAB implementation of the multivariate attention-based regularization framework for gravity inversion.
The method extends the classical Tikhonov regularization by introducing cell-wise adaptive regularization weights, which are dynamically updated based on:
  Data sensitivity (Jacobian-based sensitivity measure)
  Prior-model credibility (data-consistency evaluation)
The framework enables spatially heterogeneous utilization of prior information while maintaining stability within a deterministic inversion scheme.

# Project Structure
├── config/          % Parameter configuration files
├── data/            % Input data (Excel format)
├── model/           % Model definition files
├── src/
│   ├── forward/     % Forward modeling operators
│   ├── inversion/   % Inversion algorithms
│   └── utils/       % Utility functions
├── main.m           % Main execution script
└── README.md
Folder Description

config/ – Contains inversion parameters and runtime configuration.
data/ – Stores observed data, prior models, and experiment settings (Excel format).
model/ – Defines synthetic or field model structures.
src/forward/ – Implements gravity forward modeling.
src/inversion/ – Implements Tikhonov and multivariate attention inversion schemes.
src/utils/ – Auxiliary numerical routines.
main.m – Entry point of the inversion workflow.

# Installation

Requirements
MATLAB (R2021a or later recommended)
Microsoft Excel
No additional MATLAB toolboxes are required beyond the standard installation.

Setup Instructions

1.Download or clone this repository.
2.Open MATLAB.
3.Set the project root directory as the working folder.
4.Run:
  main

# Usage

1.Modify model parameters in config/.
2.Edit observed data or prior models in the data/ Excel files.
3.Execute main.m.
4.Results (inversion models and error maps) will be generated automatically.

# Data Format

All input data are stored in Excel (.xlsx) format, including:
  Observed gravity data
  Prior density models
  Grid definitions
  Inversion parameters
MATLAB reads the Excel files automatically during execution.

# Reproducibility 

To reproduce the results in the manuscript:
1.Use the provided Excel datasets.
2.Keep default configuration parameters.
3.Execute main.m.
All figures in the manuscript can be regenerated using the provided scripts.

# License
This project is released for academic research purposes.
Commercial use requires permission from the authors.


