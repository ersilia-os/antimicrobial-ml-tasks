import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from default import UNITS_MASTER


class RawCleaner():
    def __init__(self):
        self.ucum = UNITS_MASTER
    
    def get_comment(self, df):
        """The variable comment_active may sometimes include an indication that the result 
        of the experiment is "Active" or "Not Active".
        We will use this information in the calculation of the target variables.
        The variable comment_active will be:
            - True if activity_comment = 'active'
            - False if activity_comment = 'not active'
            - NA otherwise
        """
        df['comment_active'] = df['activity_comment'].str.upper()\
                .apply(lambda x: 
                        1 if x=='ACTIVE' else
                        0 if x=='NOT ACTIVE'
                        else np.nan
                    ).astype('float')
        return df

    def eliminate_rows(self, df):
        pref = len(df)
        df = df[~df["canonical_smiles"].isna()]
        print("removing rows with empty smiles: {}".format(pref-len(df)))
        df['standard_units'] = df['standard_units'].fillna('N/A')    
        # Remove rows where both standard_value and comment_active are null 
        # (we can't know if they are active or not)
        rows_before_filter = len(df)
        df = df[~((df.standard_value.isnull()) & (df.comment_active.isnull()) )]
        rows_after_filter = len(df)
        print('\nRemoving rows where standard_value is null and activity_comment does not inform on activity.')
        print(f'Removed {rows_before_filter-rows_after_filter} rows. Remaining {rows_after_filter} rows.')
        return df
    
    def drop_unwanted_cols(self, df):
        cols_to_drop = ['doc_id', 'assay_id', 'activity_id', 'assay_type',
       'assay_confidence_score', 'assay_bao_format', 'pchembl_value', 'activity_comment',
       'target_tax_id', 'protein_accession_class', 'year',
       'pubmed_id', 'count_activity_rows', 'doc_id_all',
       'assay_id_all', 'activity_id_all', 'assay_description']
        cols_int = list(set(cols_to_drop).intersection(df.columns))
        df = df.drop(columns=cols_int, errors='ignore')
        return df

    def _mol_weight(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return Descriptors.MolWt(mol)

    def add_mol_weight(self, df):
        df["molecular_weight"] = df["canonical_smiles"].apply(self._mol_weight)
        #if molweight cannot be calculated, drop row
        df = df[~df["molecular_weight"].isna()]
        return df
    
    def add_ucum_units(self, df):
        units_to_val_units = dict(zip(self.ucum['units'], self.ucum['val_unit']))
        #df['val_units'] = [units_to_val_units[unit] for unit in df['standard_units']]
        df['val_units'] = [units_to_val_units.get(unit, unit) for unit in df['standard_units']]
        return df
    
    def run(self, df):
        df = self.get_comment(df)
        df = self.eliminate_rows(df)
        df = self.drop_unwanted_cols(df)
        df = self.add_mol_weight(df)
        df = self.add_ucum_units(df)
        return df

class UnitStandardiser():
    #the units file has been processed in UCUM. Selected units will be converted to standard:
    #Molar: umol
    #weight/volume: ug.ml-1 to uM
    #weight/weight: ug.mg-1
    #molar/weight: umol/mg

    def __init__(self):
        self.unit_col = "standard_units"
        self.value_col = "standard_value"
        self.mw_col = "molecular_weight"
        self.ucum = UNITS_MASTER
    
    def _parse_function(self, s):
        if 'standard_value' not in s:
            return None
        if "molecular_weight" in s:
            p = "lambda x,y: "
        else:
            p = "lambda x: "
        s = s.replace("molecular_weight", "y")
        s = s.replace("standard_value", "x")
        s = p + s
        return eval(s)

    def _umol_converter(self):
        converter_str = {}
        converter_frm = {}
        for i,r in self.ucum.iterrows():
            if r["final_unit"] == "umol":
                converter_str[r["units"]] = r["transformer"]
        for k,v in converter_str.items():
            f = self._parse_function(v)
            converter_frm[k]=f
        return converter_frm

    def _ugmg_converter(self):
        converter_str = {}
        converter_frm = {}
        for i,r in self.ucum.iterrows():
            if r["final_unit"] == "ug.mg-1":
                converter_str[r["units"]] = r["transformer"]
        for k,v in converter_str.items():
            f = self._parse_function(v)
            converter_frm[k]=f
        return converter_frm
    
    def _umolmg_converter(self):
        converter_str = {}
        converter_frm = {}
        for i,r in self.ucum.iterrows():
            if r["final_unit"] == "umol.mg-1":
                converter_str[r["units"]] = r["transformer"]
        for k,v in converter_str.items():
            f = self._parse_function(v)
            converter_frm[k]=f
        return converter_frm

    def standardise(self,df):
        final_units = []
        final_value = []
        umol_converter = self._umol_converter()
        ugmg_converter = self._ugmg_converter()
        umolmg_converter = self._umolmg_converter()

        for i,r in df.iterrows():
            if r[self.unit_col] in umol_converter.keys():
                final_units += ["umol"]
                if r["val_units"] in ["umol", "nmol", "pmol", "mmol", "mol"]:
                    print("here", r[self.unit_col])
                    final_value += [umol_converter[r[self.unit_col]](r[self.value_col])] 
                else:
                    print("not here", r[self.unit_col])
                    final_value += [umol_converter[r[self.unit_col]](r[self.value_col], r[self.mw_col])]
            elif r[self.unit_col] in umol_converter.keys():
                final_units += ["ug.mg-1"]
                final_value += [ugmg_converter[r[self.unit_col]](r[self.value_col], r[self.mw_col])]
            elif r[self.unit_col] in umol_converter.keys():
                final_units += ["umol.mg-1"]
                final_value += [umolmg_converter[r[self.unit_col]](r[self.value_col], r[self.mw_col])]
            else:
                final_units += [r[self.unit_col]]
                final_value += [r[self.value_col]]
        df["final_units"] = final_units
        df["final_value"] = final_value
        return df



