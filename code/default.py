import pandas as pd
import os

abspath = os.path.dirname(os.path.abspath(__file__))

CHEMBL_USR = "chembl_user"
CHEMBL_PWD = "aaa"
UNITS_MASTER = pd.read_csv(os.path.join(abspath, "..", "config", "ucum.csv"))
DATAPATH = os.path.join(abspath, "..", "data")