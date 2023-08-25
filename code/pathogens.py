import os
import sys
import pandas as pd
import numpy as np
from chemblmltools import chembl_activity_target

from default import CHEMBL_PWD, CHEMBL_USR


class PathogenGetter():
    """
    This class obtains the full data related to the specified pathogens
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