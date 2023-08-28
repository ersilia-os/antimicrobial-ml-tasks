import os
import sys
import pandas as pd

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from default import DATASET_SIZE_LIMIT, DATAPATH

class Splitter():
    def __init__(self, pathogen):
        self.pathogen = pathogen
        self.base_path  = os.path.join(DATAPATH, pathogen)

    def _cap_dataset_size(self, df):
        """If the dataset is larger than a certain limit, get a sample. Otherwise get full dataset.
        This is done separately for each of the classes, to optimize the balance.
        """
        if len(df) <= DATASET_SIZE_LIMIT:
            return df
        else:
            class_size_limit = DATASET_SIZE_LIMIT // 2  # Size limit for each class (0/1)
            df_pos = df[df.activity==1]
            if len(df_pos) > class_size_limit:
                df_pos = df_pos.sample(n=class_size_limit)  # Randomly sample values if size too large
            df_neg = df[df.activity==0]
            if len(df_neg) > class_size_limit:
                df_neg = df_neg.sample(n=class_size_limit)  # Randomly sample values if size too large

            df_new = pd.concat([df_pos, df_neg])
            # Shuffle rows
            df_new = df_new.sample(frac=1).reset_index(drop=True)
            return df_new

    def _create_dir_if_not_exists(self, dirpath):
        """Create a directory if it does not exist"""
        if not os.path.isdir(dirpath):
            os.mkdir(dirpath)
            print('Directory created:', dirpath)
    
    def _get_available_tasks(self):
        task_files = [fn for fn in os.listdir(self.base_path) if os.path.isfile(os.path.join(self.base_path, fn)) and ("hc" in fn or "lc" in fn)]
        tasks = []
        for fn in task_files:
            t = fn.split(".")[0].split("_")[1:]
            t = "_".join(t)
            tasks += [t]
        return tasks

    def create_directoy_structure(self):
        """Create directory structure and input files for the models of current pathogen"""
        # To prevent errors, the base directory must exist previously.
        if not os.path.isdir(self.base_path):
            raise ValueError(f'The directory {self.base_path} does not exist.')

        # For each task, create directorate structure and input files
        for t in self._get_available_tasks():
            task_path = os.path.join(self.base_path, t)
            print(task_path)
            # Create directories
            self._create_dir_if_not_exists(task_path)
            self._create_dir_if_not_exists(os.path.join(task_path, 'input'))
            self._create_dir_if_not_exists(os.path.join(task_path, 'model'))
            self._create_dir_if_not_exists(os.path.join(task_path, 'test'))
            self._create_dir_if_not_exists(os.path.join(task_path, 'log'))

    def create_input_files(self):
        #task_files = [fn for fn in os.listdir(self.base_path) if "hc" in fn or "fc" in fn]
        task_files = [fn for fn in os.listdir(self.base_path) if os.path.isfile(os.path.join(self.base_path, fn)) and ("hc" in fn or "lc" in fn)]
        for fn in task_files:
            t = fn.split(".")[0].split("_")[1:]
            t = "_".join(t)
            df = pd.read_csv(os.path.join(self.base_path, fn))
            df = self._cap_dataset_size(df[["smiles", "activity"]])
            df.to_csv(os.path.join(self.base_path,t, "input", "input.csv"), index=False)