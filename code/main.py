import sys
import os
import pandas as pd

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from pathogens import PathogenGetter
from utils import UnitStandardiser, RawCleaner
from default import DATAPATH

#pathogens = sys.argv[1] #csv file indicating the pathogens to process, see example in ../config/pathogens.csv

#Obtain desired pathogen data
"""
df_pathogens = pd.read_csv(pathogens)
list_pathogen_codes = df_pathogens.pathogen_code
list_pathogen_search_text = df_pathogens.search_text
for i, patho_code in enumerate(list_pathogen_codes):
    print('------------------------------------------------------------')
    print(f'Creating data for pathogen {patho_code} ({list_pathogen_search_text[i]})')
    print('------------------------------------------------------------')
    patho_getter = PathogenGetter(list_pathogen_search_text[i])
    df = patho_getter.create_datasets_pathogen()
    if not os.path.exists(os.path.join(DATAPATH, "pathogen_original")):
        os.makedirs(os.path.join(DATAPATH, "pathogen_original"))
    df.to_csv(os.path.join(DATAPATH, "pathogen_original", "{}.csv".format(patho_code)), index=False)


# List all files in the folder
file_list = os.listdir(os.path.join(DATAPATH, "pathogen_original"))

# Iterate through each file
for filename in file_list:
    df = pd.read_csv(os.path.join(DATAPATH, "pathogen_original", filename), low_memory=False)
"""
#Standardise units

df = pd.read_csv(os.path.join(DATAPATH, "pathogen_original", "efaecium.csv"))
rc = RawCleaner()
df = rc.run(df)

us = UnitStandardiser()
df = us.standardise(df)

df.to_csv(os.path.join(DATAPATH, "test.csv"), index=False)