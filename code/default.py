import pandas as pd
import os

abspath = os.path.dirname(os.path.abspath(__file__))

CHEMBL_USR = "chembl_user"
CHEMBL_PWD = "aaa"
UNITS_MASTER = pd.read_csv(os.path.join(abspath, "..", "config", "ucum.csv"))
DATAPATH = os.path.join(abspath, "..", "data")
CONFIGPATH = os.path.join(abspath, "..", "config")
MIN_SIZE_ASSAY_TASK = 1000 #Top assays with at least this data size will get a specific task
MIN_SIZE_PROTEIN_TASK = 250 #Top proteins with at least this data size will get a specific task
MIN_COUNT_POSITIVE_CASES = 30
TOP_ASSAYS = 3
TOP_TYPES = 3
TOP_PROTEINS = 3
DATASET_SIZE_LIMIT = 1e6
SPLIT_METHOD = 'similarity'
