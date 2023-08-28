import os
import sys
import pandas as pd
import numpy as np
from chemblmltools import chembl_activity_target


abspath = os.path.abspath(__file__)
sys.path.append(abspath)
from default import CHEMBL_PWD, CHEMBL_USR, DATAPATH

pathogens = sys.argv[1] #csv file indicating the pathogens to process, see example in ../config/pathogens.csv


class PathogenGetter():
    """
    This class obtains the full data related to the specified pathogens and creates the directories to save it
    """
    def __init__(self, patho_search):
        self.db_user = CHEMBL_USR
        self.db_pwd = CHEMBL_PWD
        self.patho_search = patho_search

    def create_datasets_pathogen(self):
        df = chembl_activity_target(
                db_user= self.db_user,
                db_password=self.db_pwd,
                organism_contains=self.patho_search,
                max_heavy_atoms=100)
        return df
    

#Obtain desired pathogen data  
df_pathogens = pd.read_csv(pathogens)
list_pathogen_codes = df_pathogens.pathogen_code
list_pathogen_search_text = df_pathogens.search_text
for i, patho_code in enumerate(list_pathogen_codes):
    print('------------------------------------------------------------')
    print(f'Creating data for pathogen {patho_code} ({list_pathogen_search_text[i]})')
    print('------------------------------------------------------------')
    patho_getter = PathogenGetter(list_pathogen_search_text[i])
    df = patho_getter.create_datasets_pathogen()

    if not os.path.exists(os.path.join(DATAPATH, patho_code)):
        os.makedirs(os.path.join(DATAPATH, patho_code))
    df.to_csv(os.path.join(DATAPATH, patho_code, "{}_original.csv".format(patho_code)), index=False)
