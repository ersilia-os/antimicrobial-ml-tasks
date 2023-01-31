# antimicrobial-ml-tasks
Antimicrobial activity prediction with automated machine learning

# Process overview

![Process overview flow chart](doc/images/pipeline_overview.png)

# Installation

This process has been developed in Ubuntu Linux. For other environments, it may require adaptations.

The installation instructions assume you have Ubuntu Linux, and that the conda package manager is installed.

## Installation required for part 1 (create datasets)

1. Clone this repository (https://github.com/ersilia-os/antimicrobial_ml_tasks.git)

2. Not required but recommended: create a conda environment for this project. Activate it.

```
conda create -n antimicrobial python=3.7
conda activate antimicrobial
```

3. Install the package https://github.com/ersilia-os/chembl_ml_tools.git , following the instructions in that 
repository. This includes the instructions to install the ChEMBL database in PostgreSQL.


4. Create a directory "models" in your home. Your models and model data will be stored here.

```
mkdir ~/models
```
Note: If you prefer to use a different directory, just edit it in the variable `BASE_PATH` in the program `code/create_datasets.py`.

## Installation required for part 2 (build models)

1. Install the Ersilia Model Hub: https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/installation

2. Install ZairaChem by following the instructions in the repository: https://github.com/ersilia-os/zaira-chem

3. Copy the script `scripts/call_zairachem.sh` to the directory `~/bin` (create `~/bin` if it does not exist).

# Running part 1 - Create datasets

1. Make sure that the PostgreSQL server containing the ChEMBL database is running. In case of doubt, review step 3 of the installation.

By default the programs assume that PostgreSQL is running in the local computer, and that the database user `chembl_user` with
password `aaa` has read access to the tables of ChEMBL. This can be changed in program `code/create_datasets.py`.

2. Edit the file `config/pathogens.csv` to select the pathogens for which we need models.

This file has two columns:

- **pathogen_code**: Choose a short code to identify the pathogen, alphanumeric only, **without spaces**. Example: "efaecium".
- **search_text**: A search string, case insensitive, to search for the pathogen name in the ORGANISM field in the ChEMBL database. This may require some trial and error. Example: "Enterococcus Faecium".

3. Run the program `create_datasets.py`
```
python code/create_datasets.py
```

This will create the required directory structure under the base path (`~/models`). In the subdirectory of each model there is an `input` subdirectory. The input dataset for that model will be created there.

This will also generate the file `model_metadata/dataset.csv` containing a list of all the datasets and their counts.

# Running part 2 - Build models

1. Run the script to perform the train-test split (instructions pending)

2. Run the script to fit and assess the models (instructions pending)

# Results

The directory for each model (example: `~/models/saureus/saureus_organism_anytype`) will contain the following subdirectories:

- input: Contains the files:
  - input.csv: full input data
  - train.csv: input data for training
  - test.csv: input data for test
  - input_rejected.csv: cases that ZairaChem has rejected (typically because the molecule's SMILES is not valid)
  
- model: Contains the model definition, in the format used by ZairaChem

- test: Predictions for the test data and assessment reports of the model

- log: The log files resulting from the split, test and predict runs

