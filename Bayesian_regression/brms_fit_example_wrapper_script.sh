#!/bin/bash

### Wrapper to run R script for Bayesian Regression modelling of interactions in Swedish Littorina 
### saxatilis. This job was run on the server 'Rackham', part of the National Academic Infrastructure 
### for Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 2024-01-12

### Job parameters
#SBATCH -A naiss2023-22-881 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 3-00:00:00
#SBATCH -J brms_response_predictor
#SBATCH --output logfiles/brms_response_predictor.out
#SBATCH --error logfiles/brms_response_predictor.err

### Required software
module load R/4.2.1
module load R_packages/4.2.1

### Paramters
Pheno=...
Habitat=... # Set to 'Inv' for Pheno~Inv models
Iter=...
Chain=...

### Command
Rscript --no-restore --quite --no-save --vanilla brms_fit_MODEL.R $Pheno $Habitat $Iter $Chain
