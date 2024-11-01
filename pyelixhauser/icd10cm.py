

'''
pyelixhauser icd9cm module contains
comorbidity_from_string,
comorbity_from_array methods for encoding icd9cm codes to Elixhauser Comorbidity
sourced from



http://mchp-appserv.cpe.umanitoba.ca/concept/Elixhauser%20Comorbidities%20-%20Coding%20Algorithms%20for%20ICD-9-CM%20and%20ICD-10.pdf

this table was recreated as a csv to store data for the icd9cm Elixhauser map

Quan H, Sundararajan V, Halfon P, et al. Coding algorithms for defining Comorbidities in ICD-9-CM and
ICD-10 administrative data. Med Care. 2005 Nov; 43(11): 1130-9
example inputs

['Congestive heart failure ', 'Cardiac arrhythmias ', 'Valvular disease',
       'Pulmonary circulation Disorders', 'Peripheral vascular disorders',
       'Hypertension, uncomplicated', 'Hypertension, complicated', 'Paralysis',
       'Other neurological disorders ', 'Chronic pulmonary disease ',
       'Diabetes, uncomplicated ', 'Diabetes, complicated', 'Hypothyroidism ',
       'Renal failure', 'Liver disease',
       'Peptic ulcer disease excluding bleeding', 'AIDS/H1V ', 'Lymphoma ',
       'Metastatic cancer ', 'Solid tumor without Metastasis',
       'Rheumatoid arthritis/ collagen vascular diseases', 'Coagulopathy ',
       'Obesity ', 'Weight loss', 'Fluid and electrolyte Disorders',
       'Blood loss anemia', 'Deficiency anemia', 'Alcohol abuse ',
       'Drug abuse', 'Psychoses', 'Depression']



'''


import pkg_resources
import pandas as pd
import re
import logging
from pyelixhauser.setup_logger import logger
from pyelixhauser.utils import load_resource


_resource_name = "/resources/icd10_elixhauser.csv"
_manager = 'pyelixhauser'
_path = pkg_resources.resource_stream(_manager,_resource_name)

_icd10cm_lookup_df = pd.read_csv(_path)
_feature_names = _icd10cm_lookup_df.loc[:, 'Comorbidities']
_feature_names = [n.replace('\n', ' ' ).strip() for n in _feature_names]

def _parse_codes(s, pattern='[A-Z][0-9]+[.]?x?[0-9]*'
):
    return re.findall(pattern, s)



def _isin_range(s, min_code, max_code= None):
    s = str(s).lower().replace('\n', '').replace('\t', '').strip()
    min_code = str(min_code).lower().replace('\n', '').replace('\t', '').strip()
    if max_code is None:
        max_code = min_code
    else:
        max_code =str(max_code).lower().replace('\n', '').replace('\t', '').strip()
    logger.debug(F'checking inputs {s} is in range {min_code, max_code} ..')
    if all((s[0]==min_code[0], s[0]==max_code[0])):
        s = s[1:]
        min_code = min_code[1:]
        max_code = max_code[1:]
        try:
            v = float(s)
        except ValueError:
            logger.debug(F'input string  {s} to _isin_range in not in range {min_code, max_code}')
            return False

        if all(('x' in min_code, 'x' in max_code)):
            try:
                min_val = float(re.sub('x', '', min_code))
                max_val = float(re.sub('x', '', max_code)) + 1
                if all((v>=min_val, v< max_val)):
                    return True
                else:
                    return False
            except ValueError:
                logger.debug(F'input string  {s} to _isin_range in not in range {min_code, max_code}')
                return False

        if any(('x' in min_code, 'x' in max_code)):
            try:
                min_val = float(re.sub('x', '', min_code))
                max_val = float(re.sub('x', '', max_code)) + 1
                if all((v>=min_val, v< max_val)):
                    return True
                else:
                    return False
            except ValueError:
                logger.debug(F'input string  {s} to _isin_range in not in range {min_code, max_code}')
                return False
        else:
            try:
                min_val = float(min_code)
                max_val = float(max_code)
                if all((v>=float(min_val), v<= float(max_val))):
                    return True
                else:
                    return False
            except ValueError:
                return False
    else:
        return False

def _icd_str_look_up(input_str, ref_str):
    input_codes = _parse_codes(input_str)
    reference_codes = _parse_codes(ref_str)
    ## Single reference code
    if len(reference_codes ) == 1:
        if 'x' in reference_codes[0]:
            isin_results = [_isin_range(code,reference_codes[0],reference_codes[0]) for code in input_codes]
        else:
            isin_results = [code == reference_codes[0] for code in input_codes]
        return any(isin_results)

    elif len(reference_codes ) == 2:
        isin_results = [_isin_range(code,reference_codes[0],reference_codes[1]) for code in input_codes]
        return any(isin_results)
    elif len(reference_codes) == 0:
        return False
    else:
        raise ValueError(F" multipe reference code {reference_codes } ranges with 3 or more codes not supported")

def _icd_str_has_comorbity(input_str, ref_str):
    ref_str = ref_str.replace('\n', '')
    ref_code_ranges =[c.strip() for c in  ref_str.split(',')]
    results = list(map(lambda x: _icd_str_look_up(input_str, x),     ref_code_ranges))
    return any(results)



def comorbidity_from_string(s):
    '''
    comorbidity_from_string
    function to detect comorbities from a string if icd10cm codes
    param s: string (icd10cm codes shoulld be seperated with , or space )
    returns pandas series, index is comorbities

    example usage:
    comorbidity_from_string('K29.2 | K70.0')


    '''
    ref_str_array = _icd10cm_lookup_df.loc[:, "Enhanced ICD-9-CM"]
    results = list(map(lambda x: _icd_str_has_comorbity(s, x), ref_str_array))
    return pd.Series(results, index=_feature_names).replace({True:1, False:0})

def get_elix(s):
    '''
    get_elix(
    function to detect comorbities from a string if icd10cm codes
    param s: string (icd10cm codes shoulld be seperated with , or space )
    returns str

    example usage:
    get_elix('Z71.5')


    '''
    array = comorbidity_from_string(s)
    if array.sum() == 1:
        return array.loc[array == 1].index[0]
    elif array.sum() > 1:
        return  ' | '.join(array.loc[array == 1].index)
    else:
        return None

def comorbidity_from_array(array):
    '''
    comorbidity_from_string
    function to detect comorbities from a string if icd10cm codes
    param s: string (icd10cm codes shoulld be seperated with , or space )
    yield pandas series, index is comorbities from ICD9 cm codes

    example usage:
    comorbidity_from_string(' F34.1','K29.2 | K70.0')


    '''

    results = list(map(comorbidity_from_string, array))
    return pd.DataFrame(results, columns=_feature_names)
       
def cumulative_comorbidity_plot(patient_icd_array, event_dates):
    '''
    Plots the cumulative number of comorbidities for a specific patient over time, using ICD-9 codes.
    
    Parameters:
    - patient_icd_array (list of str): List of strings containing ICD-9 codes (one string per event).
    - event_dates (list of str or datetime-like): List of event dates (same length as patient_icd_array).
    
    Returns:
    - None: Displays a plot of the cumulative comorbidities over time.
    
    Example Usage:
    patient_icd_array = ['490', '506', 'E11.9']
    event_dates = ['2021-01-01', '2021-06-15', '2022-03-20']
    cumulative_comorbidity_plot(patient_icd_array, event_dates)
    '''
    
    # Convert event_dates to datetime if they are not already
    if not pd.api.types.is_datetime64_any_dtype(event_dates):
        event_dates = pd.to_datetime(event_dates, errors='coerce')

    assert len(patient_icd_array) == len(event_dates), "The number of ICD-9 codes must match the number of event dates."
    
    # Initialize a set to store unique comorbidities
    unique_comorbidities = set()
    cumulative_counts = []

    # Loop through each event
    for i, icd_str in enumerate(patient_icd_array):
        # Add new ICD codes to the set of unique comorbidities
        unique_comorbidities.add(icd_str)
        
        # Store the current cumulative count of comorbidities
        cumulative_counts.append(len(unique_comorbidities))
    
    plt.figure(figsize=(10, 6))
    plt.plot(event_dates, cumulative_counts, marker='o', linestyle='-', color='b')
    plt.title("Cumulative Comorbidities Over Time")
    plt.xlabel("Event Dates")
    plt.ylabel("Cumulative Comorbidities")
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def preprocess_dates(df, date_columns):
    ''' Converts and localizes date columns in the DataFrame. 
    
    Parameters:
    - df (pd.DataFrame): The DataFrame containing date columns.
    - date_columns (list of str): List of column names that need to be converted to datetime.
    
    Returns:
    - pd.DataFrame: DataFrame with specified date columns converted to datetime and localized.
    
    Example Usage:
    df = pd.DataFrame({
        'condition_start_datetime': ['2021-01-01', '2021-02-15', '2021-03-20'],
        'date_of_birth': ['1980-05-10', '1990-06-15', '2000-07-20']
    })
    date_columns = ['condition_start_datetime', 'date_of_birth']
    df = preprocess_dates(df, date_columns)
    '''
    for column in date_columns:
        df[column] = pd.to_datetime(df[column], errors='coerce')
        df[column] = df[column].dt.tz_localize(None)
    return df


def calculate_cumulative_comorbidities(person_id, data, comorbidity_function, interval=5):
    '''
    Calculates the cumulative comorbidities for a specific patient over defined intervals using ICD-9 codes.

    Parameters:
    - person_id (int or str): The unique identifier for the patient.
    - data (pd.DataFrame): A DataFrame containing patient data with date and ICD-9 code information.
      Must include columns: 'person_id', 'condition_start_datetime', and 'source_concept_code'.
    - comorbidity_function (function): A function that takes a string of ICD codes and returns a comorbidity index.
    - interval (int): The number of years for each time interval. Default is 5 years.

    Returns:
    - List of dicts: Each dict contains 'person_id', 'interval_start_year', and 'comorbidity_index' for each time interval.
    
    Example Usage:
    df = pd.DataFrame({
        'person_id': [1060119, 1060119, 1060119],
        'condition_start_datetime': ['2021-01-01', '2021-06-15', '2022-03-20'],
        'source_concept_code': ['I10', 'E11.9', 'K70.0']
    })

    def calculate_comorbidity_index(icd_codes_string):
        return len(icd_codes_string.split())

    results = calculate_cumulative_comorbidities(1060119, df, calculate_comorbidity_index, interval=5)
    for result in results:
        print(result)
    '''
    results = []
    accumulated_icd_codes = set()
    
    # Ensure date columns are preprocessed
    data = preprocess_dates(data, ['condition_start_datetime'])
    
    start_year = int(data['condition_start_datetime'].dt.year.min())
    end_year = int(data['condition_start_datetime'].dt.year.max())
    
    for year in range(start_year, end_year + 1, interval):
        start_date = pd.Timestamp(f'{year}-01-01')
        end_date = pd.Timestamp(f'{year + interval - 1}-12-31')
        
        interval_data = data[(data['person_id'] == person_id) & 
                             (data['condition_start_datetime'] >= start_date) & 
                             (data['condition_start_datetime'] <= end_date)]
        
        if interval_data.empty:
            continue
        
        # Accumulate ICD codes
        accumulated_icd_codes.update(interval_data['source_concept_code'].astype(str))
        
        # Join accumulated ICD codes into a single string
        icd_codes_string = ' '.join(accumulated_icd_codes)
        
        # Calculate the comorbidity index using the cumulative ICD codes
        comorbidity_index = comorbidity_function(icd_codes_string)
        
        # Ensure comorbidity_index is a scalar
        if isinstance(comorbidity_index, (pd.Series, list)):
            comorbidity_index = comorbidity_index.sum()
        
        results.append({'person_id': person_id, 'interval_start_year': year, 'comorbidity_index': comorbidity_index})
    
    return results

def plot_cumulative_comorbidities_for_person(person_id, data, comorbidity_function, interval=5):
    '''
    Plots the cumulative comorbidity index for a specific patient over defined intervals using ICD-9 codes.

    Parameters:
    - person_id (int or str): The unique identifier for the patient.
    - data (pd.DataFrame): A DataFrame containing patient data with date and ICD-9 code information.
      Must include columns: 'person_id', 'condition_start_datetime', and 'source_concept_code'.
    - comorbidity_function (function): A function that takes a string of ICD codes and returns a comorbidity index.
    - interval (int, optional): The number of years for each time interval. Default is 5 years.

    Returns:
    - None: Displays a plot of the cumulative comorbidity index over time.
    
    Example Usage:
    df = pd.DataFrame({
        'person_id': [1060119, 1060119, 1060119],
        'condition_start_datetime': ['2021-01-01', '2021-06-15', '2022-03-20'],
        'source_concept_code': ['I10', 'E11.9', 'K70.0']
    })

    def calculate_comorbidity_index(icd_codes_string):
        return len(icd_codes_string.split())

    plot_cumulative_comorbidities_for_person(1060119, df, calculate_comorbidity_index, interval=5)
    '''
    # Calculate cumulative comorbidities over time
    results = calculate_cumulative_comorbidities(person_id, data, comorbidity_function, interval)
    
    if not results:
        print(f"No data available for person_id {person_id}")
        return
    
    comorbidity_intervals_df = pd.DataFrame(results)

    # Plot the results
    plt.figure(figsize=(14, 8))
    sns.lineplot(data=comorbidity_intervals_df, x='interval_start_year', y='comorbidity_index', marker='o')
    plt.title(f'Cumulative Comorbidity Index Over Time for Person ID {person_id}')
    plt.xlabel('Year')
    plt.ylabel('Comorbidity Index')
    plt.show()
