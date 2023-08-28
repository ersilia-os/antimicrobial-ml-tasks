import sys
import os
import pandas as pd

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from utils import UnitStandardiser, RawCleaner, Binarizer
from generate_datasets import AssayDatasets, ProteinDatasets, TypeDatasets
from default import DATAPATH
from split_datasets import Splitter


pathogen = "efaecium"
"""
#Standardise units
df = pd.read_csv(os.path.join(DATAPATH, pathogen, "{}_original.csv".format(pathogen)), low_memory=False)
rc = RawCleaner()
df = rc.run(df)

us = UnitStandardiser()
df = us.standardise(df)
df.to_csv(os.path.join(DATAPATH, pathogen, "{}_processed.csv".format(pathogen)), index=False)

bin = Binarizer()
df = bin.run(df)
df.to_csv(os.path.join(DATAPATH, pathogen, "{}_binary.csv".format(pathogen)), index=False)


td = TypeDatasets(pathogen)
lc, hc = td.run_any(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anytype_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anytype_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} assay types".format(pathogen))

types_data = td.run_type(df)
for k,v in types_data.items():
    try:
        v[0].to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_lc.csv".format(pathogen, k)), index=False)
    except:
        print("No Low Cut data for {} {} assay".format(pathogen, k))
    try:           
        v[1].to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_hc.csv".format(pathogen, k)), index=False)
    except:
        print("No High Cut data for {} {} assay".format(pathogen, k))

"""
s = Splitter(pathogen)
s.create_directoy_structure()
s.create_input_files()