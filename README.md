# antimicrobial-ml-tasks
Antimicrobial activity prediction with automated machine learning

# Process overview

![Process overview flow chart](doc/images/pipeline_overview.png)

(this document is work in progress)

# Installation

This process has been developed in Linux. For other environments, it may require adaptations.

The installation instructions asume you have Linux, and that the conda package manager is installed.

1. Not required but recommended: create a conda environment for this project. Activate it.

```
conda create -n antimicrobial python=3.7
conda activate antimicrobial
```

2. Clone the repository https://github.com/ersilia-os/chembl_ml_tools.git
Follow the installation instructions in that repository. This will include installing the ChEMBL database.

3. Create a directory "models" in your home. Your models and model data will be stored here.

```
mkdir ~/models
```

Note: If you prefer to use a different directory, just edit it in the variable `BASE_PATH` in the program `code/create_datasets.py`.

4. By default, the program `create_datasets.py` will create the datasets required for the models of 6 pathogens known as *ESKAPE*. If you need different pathogens, edit the variables `LIST_PATHOGEN_CODES` and `LIST_PATHOGEN_SEARCH_TEXT` in the program `code/create_datasets.py`.


# Running part 1 - Create datasets

1. Run the program `create_datasets.py`

```
python code/create_datasets.py
```

This will create the required directory structure under the base path (`~/models`). In the subdirectory of each model there is an `input` subdirectory. The input dataset for that model will be there.
