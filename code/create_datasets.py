import pandas as pd
import numpy as np
import os
from chemblmltools import chembl_activity_target
from chemblmltools import ChemblMoleculeSampler

# Set BASE_PATH to the base directory where all the model data will be stored
BASE_PATH = os.path.expanduser('~/models')

# If the data is large, in a pilot modeling phase we may wish to build models with a sample of the data.
# To do that, change DATASET_SIZE_LIMIT to the desired maximum value, e.g. 2000 or 5000
# If we do not want to limit, we set it to 1e6
DATASET_SIZE_LIMIT = 1e6

# Split method (similarity, scaffold or random): will be used in the generated script split_all.sh
SPLIT_METHOD = 'similarity'


def create_dataset_for_task(moleculeSampler, df_input, target_var):
    """Create a dataset that will be used to train a model
    - Start with df_input
    - Rename as "activity" the indicated variable 
      (normally this will be 'activity_lc' or 'activity_hc')
    - Keep only the required output variables
    - Remove duplicate molecules. Keep the most frequent value of activity (1 in case of tie)

    Returns:
        Dataset
        (If very few postitive cases (less than MIN_COUNT_POSITIVE_CASES), return None)
    
    """
    # Tasks with very few positive cases will be discarded
    MIN_COUNT_POSITIVE_CASES = 30

    # Create dataset with required columns
    df = df_input.rename(columns={target_var:'activity'})[
        ['compound_chembl_id', 'canonical_smiles', 'activity']]
        
    # Remove duplicates. Within each group, activity takes the majority value (1 in case of tie)
    # To do this, we calculate the mean of activity, then assign 1 if the mean is 0.5 or higher    
    df = df.groupby(['compound_chembl_id', 'canonical_smiles'])['activity'].mean().reset_index()
    df['activity'] = df['activity'].apply(lambda x: 1 if x>=0.5 else 0)

    count_pos = len(df[df.activity==1])
    count_neg = len(df[df.activity==0])
    
    # If there are fewer negatives than positives, fill up negatives with a random sample
    if count_neg < count_pos:
        df_sample_neg = moleculeSampler.negative_sample(
                num_molecules=count_pos-count_neg,
                list_positive_molecules=list(df[df.activity==1].compound_chembl_id))
        df_sample_neg.rename(columns={'chembl_id':'compound_chembl_id'}, inplace=True)
        df_sample_neg['activity'] = 0
        df = pd.concat([df, df_sample_neg], ignore_index=True)
    
    # Provide dataset only if minimum number of positive cases achieved.
    if count_pos >= MIN_COUNT_POSITIVE_CASES:
        df.rename(columns={'canonical_smiles':'smiles'}, inplace=True)
        return df.copy()
    else:
        return None


def store_dataset(dict_task_dataset, dict_task_desc, dataset, task_code, task_desc):
    """Store the provided dataset in dictionary dict_task_dataset
       Store the task description in dictionary dict_task_desc
       If the dataset is empty, do nothing
       """
    if dataset is not None:
        dict_task_dataset[task_code] = dataset
        dict_task_desc[task_code] = task_desc
        print(f'[{task_code}] - Created dataset for task. Description:')
        print(f'{task_desc}\n')


def warning_for_missing_config_entries(df):
    """Give a warning if missing important (type,units) combination in config table"""

    # Top 10 combinations of (type, units) that appear at least 100 times
    top_values_type_unit = df[['standard_type', 'standard_units']].value_counts()[0:10]
    top_values_type_unit = top_values_type_unit[top_values_type_unit>100]
    # Detailed data for those (type, units)
    df_top = df.merge(top_values_type_unit.to_frame(), left_on=['standard_type', 'standard_units'], right_index=True)
    # Combinations not present in config table ()
    missing_top_type_units = df_top[df_top.is_in_config_table=='left_only']\
            [['standard_type', 'standard_units']].value_counts()
    if len(missing_top_type_units) > 0:
        print('\n--- WARNING - The following combinations of (standard_type, standard_unit) are'
              ' often present in the data but do not exist in the configuration table. Please'
              ' consider updating the configuration table:')
        print(missing_top_type_units)
        print()


def create_datasets_pathogen(pathogen_code, organism_contains):
    """Create datasets for given pathogen
    
    For selected pathogen, obtain and process required data, then create
    all the defined datasets (one for each model).
    """

    ###### Read data
    
    # Get activity_target data from chembl
    print(f'Current pathogen: {pathogen_code}. Searching string in Chembl "{organism_contains}".')
    df1 = chembl_activity_target(
            db_user='chembl_user',
            db_password='aaa',
            organism_contains=organism_contains,
            max_heavy_atoms=100)
    
    # Import standard_type_config
    df_standard_type_config = pd.read_csv(
            '../config/standard_type_config.csv',
            usecols=['standard_type', 'standard_units', 'active_direction', 'low_cut', 'high_cut'],
            keep_default_na=False, na_values=['']
            )
    print('\nstandard_type_config shape:', df_standard_type_config.shape)
    
    ###### Data cleaning
    
    # Clean values of standard_units
    df1['standard_units'] = np.where(df1.standard_units=='ug ml-1', 'ug.mL-1',
                                     df1.standard_units)
    df1['standard_units'] = np.where(df1.standard_units=='ppm', 'p.p.m.',
                                     df1.standard_units)
    df1['standard_units'] = np.where(df1.standard_units.str.lower()=='mg kg-1', 'mg.kg-1',
                                     df1.standard_units)

    # standard_units: convert NaN into the text 'N/A'
    df1['standard_units'] = df1['standard_units'].fillna('N/A')    
    
    ###### Generate calculated variables
    
    # Define comment_active based on the value of activity_comment (case insensitive)
    """The variable comment_active may sometimes include an indication that the result 
        of the experiment is "Active" or "Not Active".
        We will use this information in the calculation of the target variables.
        The variable comment_active will be:
            - True if activity_comment = 'active'
            - False if activity_comment = 'not active'
            - NA otherwise
    """
    df1['comment_active'] = df1['activity_comment'].str.upper()\
            .apply(lambda x: 
                       1 if x=='ACTIVE' else
                       0 if x=='NOT ACTIVE'
                       else np.nan
                  ).astype('float')

    # Remove rows where both standard_value and comment_active are null 
    # (we can't know if they are active or not)
    rows_before_filter = len(df1)
    df1 = df1[~( (df1.standard_value.isnull()) & (df1.comment_active.isnull()) )]
    rows_after_filter = len(df1)
    print('\nRemoving rows where standard_value is null and activity_comment does not inform on activity.')
    print(f'Removed {rows_before_filter-rows_after_filter} rows. Remaining {rows_after_filter} rows.')
    
    ###### Add variables from config table
    # We will do a left join to get the required variables from the configuration table:
    # - active_direction
    # - low_cut
    # - high_cut
    df = df1.merge(df_standard_type_config, how='left',
              on=['standard_type', 'standard_units'],
              indicator='is_in_config_table')
    
    # Issue warning if any frequent combinations of (standard_type, standard_units) does not exist
    # in the configuration data.
    warning_for_missing_config_entries(df)
        
    # Exclude cases with no active_direction
    rows_before_filter = len(df1)
    df = df[df.active_direction.notnull()]
    rows_after_filter = len(df1)
    print('\nRemoving rows where the combination of standard_type and standard_units is not in the config table.')
    print(f'Removed {rows_before_filter-rows_after_filter} rows. Remaining {rows_after_filter} rows.')
    
    # Generate the two target variables
    # - ind_active_lc: Active indicator low-confidence, based on low_cut
    # - ind_active_hc: Active indicator high-confidence, based on high_cut
    
    def calculate_active(row, cut):
        if np.isnan(row.standard_value):
            return row.comment_active
        elif row.active_direction == 1:  # Higher value is more active
            if row.standard_value >= cut:
                return 1
            else:
                return 0
        elif row.active_direction == -1:  # Lower value is more active
            if row.standard_value <= cut:
                return 1
            else:
                return 0
        
    df['activity_lc'] = df.apply(
            lambda row: calculate_active(row, row.low_cut), 
            axis=1).astype('float')

    df['activity_hc'] = df.apply(
            lambda row: calculate_active(row, row.high_cut), 
            axis=1).astype('float')
    
    # Remove cases where standard_direction is not consistent with our activity definition
    # Background:
    # If the value resulting of an experiment is beyond the range that can be measured, 
    # instead of reporting the value, it will be reported as ">x" or "<x".

    # The variable standard_direction contains "=" if the precise value is reported. It will contain "<", "<=", ">" or ">=" if the reported value is a lower or upper bound. For example, if the "real" value is 500 but only up to 100 can be measured, then standard_value=100 and standard_direction=">".

    # Taking this into account, it makes sense that, for results that we label as ACTIVE:
    # - If active_direction=1
    #     * standard_relation may be '>' (it indicates a "large" value)
    # - If active_direction=-1
    #     * standard_relation may be '<' (it indicates a "small" value)

    # For results that we label as NOT ACTIVE:
    # - If active_direction=1
    #     * standard_relation may be '<' (it indicates a "small" value)
    # - If active_direction=-1
    #     * standard_relation may be '>' (it indicates a "large" value)

    print()
    
    rows_to_drop = df[(df.comment_active.isnull()) &
       (df.activity_hc==0) & 
       (df.active_direction==1) & 
       (df.standard_relation.isin(['>','>=']))
      ].index
    df.drop(rows_to_drop, inplace=True)
    print(f'Removed {len(rows_to_drop)} cases with active direction +, relation ">", but labeled as not active')

    rows_to_drop = df[(df.comment_active.isnull()) &
       (df.activity_hc==0) & 
       (df.active_direction==-1) & 
       (df.standard_relation.isin(['<','<=']))
      ].index
    df.drop(rows_to_drop, inplace=True)
    print(f'Removed {len(rows_to_drop)} cases with active direction -, relation "<", but labeled as not active')

    rows_to_drop = df[(df.comment_active.isnull()) &
       (df.activity_hc==1) & 
       (df.active_direction==1) & 
       (df.standard_relation.isin(['<','<=']))
      ].index
    df.drop(rows_to_drop, inplace=True)
    print(f'Removed {len(rows_to_drop)} cases with active direction +, relation "<", but labeled as active')

    rows_to_drop = df[(df.comment_active.isnull()) &
       (df.activity_hc==1) & 
       (df.active_direction==-1) & 
       (df.standard_relation.isin(['>','>=']))
      ].index
    df.drop(rows_to_drop, inplace=True)
    print(f'Removed {len(rows_to_drop)} cases with active direction -, relation ">", but labeled as active')

    print('Cases remaining after filter:', len(df))

    ###### Generate a data set for each task
    
    df_organism = df[df.target_type == 'ORGANISM']
    df_protein = df[df.target_type.str.contains('PROTEIN')]
    
    # Top assays with at least this data size will get a specific task
    MIN_SIZE_ASSAY_TASK = 1000
    # Top proteins with at least this data size will get a specific task
    MIN_SIZE_PROTEIN_TASK = 250
    
    # Find top 3 assays targetting organism. Must have at least MIN_SIZE_ASSAY_TASK activities
    assay_counts = df_organism.assay_id.value_counts()[:3]
    assay_counts = assay_counts[assay_counts >= MIN_SIZE_ASSAY_TASK]
    list_top_assays = list(assay_counts.index)
    print(f'\n{pathogen_code} - selected top assays:', list_top_assays)

    # Find top 5 target proteins. Must have at least MIN_SIZE_PROTEIN_TASK activities
    protein_counts = df_protein.target_pref_name.value_counts()[:5]
    protein_counts = protein_counts[protein_counts >= MIN_SIZE_PROTEIN_TASK]
    list_top_proteins = list(protein_counts.index)
    print(f'\n{pathogen_code} - selected top proteins:', list_top_proteins)

    ###### Create a data set for each task

    # Initialize molecule sampler object. Will be used to sample random molecules
    # from Chembl when there are not enough negative cases.
    # The data path is just a directory to save the data extraction of all Chembl
    # molecules and avoid obtaining it more than once.
    moleculeSampler = ChemblMoleculeSampler(
            data_path='../tmp', db_user='chembl_user', db_password='aaa',
            max_heavy_atoms=100)

    dict_task_dataset = {}  # Create empty dictionary to store the datasets
    dict_task_desc = {}  # Create empty dictionary to store the task descriptions

    # Create datasets for main modeling tasks

    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism,
            target_var='activity_lc',)
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_anytype',
            task_desc='All experiments targeting organism. Using activity low-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism,
            target_var='activity_hc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_anytype_hc',
            task_desc='All experiments targeting organism. Using activity high-confidence cut')

    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type == 'MIC'],
            target_var='activity_lc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_mic',
            task_desc='Experiments of type MIC targeting organism. Using activity low-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type == 'MIC'],
            target_var='activity_hc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_mic_hc',
            task_desc='Experiments of type MIC targeting organism. Using activity high-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type == 'IZ'],
            target_var='activity_lc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_iz',
            task_desc='Experiments of type IZ targeting organism. Using activity low-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type == 'IZ'],
            target_var='activity_hc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_iz_hc',
            task_desc='Experiments of type IZ targeting organism. Using activity high-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type.isin(['Inhibition', 'Activity'])],
            target_var='activity_lc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_inhibition',
            task_desc='Experiments of type Activity or Inhibition, targeting organism. '
                    'Using activity low-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type.isin(['Inhibition', 'Activity'])],
            target_var='activity_hc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_inhibition_hc',
            task_desc='Experiments of type Activity or Inhibition, targeting organism. '
                    'Using activity high-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type == 'IC50'],
            target_var='activity_lc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_ic50',
            task_desc='Experiments of type IC50 targeting organism. Using activity low-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_organism[df_organism.standard_type == 'IC50'],
            target_var='activity_hc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='organism_ic50_hc',
            task_desc='Experiments of type IC50 targeting organism. Using activity high-confidence cut')

    df = create_dataset_for_task(moleculeSampler,
            df_input=df_protein,
            target_var='activity_lc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='protein_any',
            task_desc='All experiments targeting proteins. Using activity low-confidence cut')
    
    df = create_dataset_for_task(moleculeSampler,
            df_input=df_protein,
            target_var='activity_hc')
    store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
            task_code='protein_any_hc',
            task_desc='All experiments targeting proteins. Using activity high-confidence cut')

    # Create datasets for top assays
    for i, assay_id in enumerate(list_top_assays):
        df = create_dataset_for_task(moleculeSampler,
                df_input=df_organism[df_organism.assay_id == assay_id],
                target_var='activity_lc')
        store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
                task_code=f'organism_assay_top{i+1}',
                task_desc=f'Results of top {i+1} assay (id={assay_id}). Using activity low-confidence cut')
        
        df = create_dataset_for_task(moleculeSampler,
                df_input=df_organism[df_organism.assay_id == assay_id],
                target_var='activity_hc')
        store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
                task_code=f'organism_assay_top{i+1}_hc',
                task_desc=f'Results of top {i+1} assay (id={assay_id}). Using activity high-confidence cut')
    
    # Create datasets for top proteins
    for i, protein_name in enumerate(list_top_proteins):
        dataset_code = f'protein_top{i+1}'
        df = create_dataset_for_task(moleculeSampler,
                df_input=df_protein[df_protein.target_pref_name == protein_name],
                target_var='activity_lc')
        store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
                task_code=dataset_code,
                task_desc=f'Experiments targeting top {i+1} protein ({protein_name}). Using activity low-confidence cut')
        
        dataset_code = f'protein_top{i+1}_hc'
        df = create_dataset_for_task(moleculeSampler,
                df_input=df_protein[df_protein.target_pref_name == protein_name],
                target_var='activity_hc')
        store_dataset(dict_task_dataset=dict_task_dataset, dict_task_desc=dict_task_desc, dataset=df,
                task_code=dataset_code,
                task_desc=f'Experiments targeting top {i+1} protein ({protein_name}). Using activity high-confidence cut')
        
    return dict_task_dataset, dict_task_desc


def cap_dataset_size(df, dataset_size_limit):
    """If the dataset is larger than a certain limit, get a sample. Otherwise get full dataset.
       This is done separately for each of the classes, to optimize the balance.
    """
    
    if len(df) <= dataset_size_limit:
        return df
    
    class_size_limit = dataset_size_limit // 2  # Size limit for each class (0/1)
    
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


def create_directoy_structure(base_path, patho_code, dict_task_dataset):
    """Create directory structure and input files for the models of current pathogen"""

    # To prevent errors, the base directory must exist previously.
    if not os.path.isdir(base_path):
        raise ValueError(f'The directory {base_path} does not exist.')

    def create_dir_if_not_exists(path):
        """Create a directory if it does not exist"""
        if not os.path.isdir(path):
            os.mkdir(path)
            print('Directory created:', path)

    # Create directory for current pathogen
    patho_path = os.path.join(base_path, patho_code)
    create_dir_if_not_exists(patho_path)

    # For each task, create directorate structure and input files
    for task_code in (dict_task_dataset):
        patho_task_code = patho_code + '_' + task_code
        task_path = os.path.join(base_path, patho_code, patho_task_code)

        # Create directories
        create_dir_if_not_exists(task_path)
        create_dir_if_not_exists(os.path.join(task_path, 'input'))
        create_dir_if_not_exists(os.path.join(task_path, 'model'))
        create_dir_if_not_exists(os.path.join(task_path, 'test'))
        create_dir_if_not_exists(os.path.join(task_path, 'log'))

        # Create input files (full, train and test)
        df_current = cap_dataset_size(dict_task_dataset[task_code][['smiles', 'activity']],
                                      dataset_size_limit=DATASET_SIZE_LIMIT)
        df_current.to_csv(os.path.join(task_path, 'input', 'input.csv'), index=False)


def create_model_metadata(dict_pathogen_task_dataset, dict_pathogen_task_desc):
    """Create the metadata tables "dataset.csv" and "model.csv"
        - dataset.csv contains the dataset master table (list of datasets and their counts)
        - model.csv contains one row for each model. Useful for systematic analysis 
          of results (see examples in analysis directory)
    """

    # The "cap" field contains the value of DATASET_SIZE_LIMIT, or "full" if we
    # do not want to cap (that is, if we set DATASET_SIZE_LIMIT>=1e6)
    if DATASET_SIZE_LIMIT >= 1e6:
        cap_value = 'full'
    else:
        cap_value = str(DATASET_SIZE_LIMIT)

    dataset_list = []
    model_list = []
    # Loop over all pathogens
    for _, patho_code in enumerate(dict_pathogen_task_dataset.keys()):
        # With pathogen patho_code, loop over all tasks
        for task_code in (dict_pathogen_task_dataset[patho_code]):
            df_current = dict_pathogen_task_dataset[patho_code][task_code]
            total_cases = len(df_current)
            positive_cases = len(df_current[df_current.activity==1])
            task_desc = dict_pathogen_task_desc[patho_code][task_code]
            dataset_list.append([patho_code, task_code, total_cases, positive_cases, task_desc])
            model_list.append([cap_value, SPLIT_METHOD, 'zairachem', patho_code, task_code])

    # Create dataset.csv
    df_dataset = pd.DataFrame(dataset_list, 
            columns=['patho_code', 'task_code', 'total_cases', 'positive_cases', 'task_desc'])
    filename_dataset = '../model_metadata/dataset.csv'
    df_dataset.to_csv(filename_dataset, index=False)
    print(f'Created file {os.path.abspath(filename_dataset)}')

    # Create model.csv
    df_model = pd.DataFrame(model_list, 
            columns=['cap', 'split', 'is_lazy', 'patho_code', 'task_code'])
    filename_model = '../model_metadata/model.csv'
    df_model.to_csv(filename_model, index=False)
    print(f'Created file {os.path.abspath(filename_model)}')


def create_scripts(dict_pathogen_task_dataset, dict_pathogen_task_desc, destination_path):
    """Create the shell scripts to run ZairaChem once for each model
    Two scripts will be created: split_all.sh and fit_predict_all.sh

    Example of fit_split_all.sh contents:
        cd efaecium/efaecium_organism_anytype
        ../../bin/call_zairachem.sh fit
        ../../bin/call_zairachem.sh predict
        cd ../..
    """

    # Open the two files to be created: split_all and fit_predict_all shell scripts
    filename_split = os.path.join(destination_path, 'split_all.sh')
    filename_fit = os.path.join(destination_path, 'fit_predict_all.sh')
    f_split = open(filename_split, 'w')
    f_fit = open(filename_fit, 'w')
    
    # Loop over all pathogens
    for _, patho_code in enumerate(dict_pathogen_task_dataset.keys()):
        # With pathogen patho_code, loop over all tasks
        for task_code in (dict_pathogen_task_dataset[patho_code]):
            patho_task_code = patho_code + '_' + task_code
            task_rel_path = os.path.join(patho_code, patho_task_code)

            f_split.writelines([
                'cd ' + task_rel_path +'\n',
                f'../../bin/call_zairachem.sh split {SPLIT_METHOD}\n',
                'cd ../..\n\n'
                ])

            f_fit.writelines([
                'cd ' + task_rel_path + '\n',
                '../../bin/call_zairachem.sh fit\n',
                '../../bin/call_zairachem.sh predict\n',
                'cd ../..\n\n',
                ])

    f_split.close()
    f_fit.close()
    print(f'Created file {os.path.abspath(filename_split)}')
    print(f'Created file {os.path.abspath(filename_fit)}')


#################### MAIN ########################################


# Run dataset creation for all pathogens and tasks

# Read required pathogens from file pathogens.csv
df_pathogens = pd.read_csv('../config/pathogens.csv')
print(f'Number of required pathogens: {len(df_pathogens)} (configured in file config/pathogens.csv)')
list_pathogen_codes = df_pathogens.pathogen_code
list_pathogen_search_text = df_pathogens.search_text

# Dictionaries to contain, for each pathogen, all their tasks and datasets
dict_pathogen_task_dataset = {}
dict_pathogen_task_desc = {}

# Run for all pathogens
for i, patho_code in enumerate(list_pathogen_codes):
    print('------------------------------------------------------------')
    print(f'Creating data for pathogen {patho_code} ({list_pathogen_search_text[i]})')
    print('------------------------------------------------------------')
    
    dict_task_dataset, dict_task_desc = create_datasets_pathogen(
            pathogen_code=patho_code,
            organism_contains=list_pathogen_search_text[i])

    dict_pathogen_task_dataset[patho_code] = dict_task_dataset
    dict_pathogen_task_desc[patho_code] = dict_task_desc
    
    create_directoy_structure(base_path=BASE_PATH,
                              patho_code=patho_code,
                              dict_task_dataset=dict_task_dataset)    

# Create "dataset.csv" and "model.csv".
create_model_metadata(dict_pathogen_task_dataset, dict_pathogen_task_desc)

# Create the scripts "split_all.sh" and "fit_predict_all.sh"
create_scripts(dict_pathogen_task_dataset, dict_pathogen_task_desc, BASE_PATH)

print('DONE')
