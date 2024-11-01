

'''
pyelixhauser icd9cm module contains
comorbidity_from_string,
comorbity_from_array methods for encoding icd9cm codes to Elixhauser Comorbidity
sourced from



http://mchp-appserv.cpe.umanitoba.ca/concept/Elixhauser%20Comorbidities%20-%20Coding%20Algorithms%20for%20ICD-9-CM%20and%20ICD-10.pdf

this table was recreated as a csv to store data for the icd9cm Elixhauser map

Quan H, Sundararajan V, Halfon P, et al. Coding algorithms for defining Comorbidities in ICD-9-CM and
ICD-10 administrative data. Med Care. 2005 Nov; 43(11): 1130-9
example inputs ; 532.91, 533.41,

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


_resource_name = "/resources/icd9cm_elixhauser.csv"
_manager = 'pyelixhauser'
_path = pkg_resources.resource_stream(_manager,_resource_name)

_icd9cm_lookup_df = pd.read_csv(_path)
_feature_names = _icd9cm_lookup_df.loc[:, 'Comorbidities']
_feature_names = [n.replace('\n', ' ' ).strip() for n in _feature_names]



def _parse_codes(s, pattern='V?[0-9]+[.]?x?[0-9]*'):
    return re.findall(pattern, s)

assert _parse_codes('490.x-505.x') == ['490.x', '505.x']

def _isin_range(s, min_code, max_code):
    s = str(s).lower().replace('\n', '').replace('\t', '').strip()
    min_code = str(min_code).lower().replace('\n', '').replace('\t', '').strip()
    max_code =str(max_code).lower().replace('\n', '').replace('\t', '').strip()

    if all((s.startswith('v'), min_code.startswith('v'), max_code.startswith('v'))):
        min_code.replace('v', '')
        max_code.replace('v', '')
        s.replace('v', '')
    else:
        pass
    if all((s.startswith('v'), min_code.startswith('v'), max_code.startswith('v'))) == False:
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
        else:
            try:
                min_val = float(min_code)
                max_val = float(max_code)
                if all((v>=float(min_val), v< float(max_val))):
                    return True
                else:
                    return False
            except ValueError:
                return False
    else:
        return False

def _icd9_str_look_up(input_str, ref_str):
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

def _icd9_str_has_comorbity(input_str, ref_str):
    ref_str = ref_str.replace('\n', '')
    ref_code_ranges =[c.strip() for c in  ref_str.split(',')]
    results = list(map(lambda x: _icd9_str_look_up(input_str, x),     ref_code_ranges))
    return any(results)



def comorbidity_from_string(s):
    '''
    comorbidity_from_string
    function to detect comorbities from a string if icd10cm codes
    param s: string (icd9cm codes shoulld be seperated with , or space )
    returns pandas series, index is comorbities

    example usage:
    comorbidity_from_string('490.1 | 506.0')


    '''
    ref_str_array = _icd9cm_lookup_df.loc[:, "Enhanced ICD-9-CM"]
    results = list(map(lambda x: _icd9_str_has_comorbity(s, x), ref_str_array))
    return pd.Series(results, index=_feature_names).replace({True:1, False:0})

def get_elix(s):
    '''
    get_elix(
    function to detect comorbities from a string if icd10cm codes
    param s: string (icd9cm codes shoulld be seperated with , or space )
    returns str

    example usage:
    get_elix(('490.1')


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
    comorbidity_from_string('490.1 | 506.0', '175', 'V45.1')




    '''

    results = list(map(comorbidity_from_string, array))
    return pd.DataFrame(results, columns=_feature_names)
def cumulative_comorbidity_plot(patient_icd_array, event_dates):
    '''
    Plots the cumulative number of comorbidities for a specific patient over time.
    
    Parameters:
    - patient_icd_array (list of str): List of strings containing ICD-9 codes (one string per event).
    - event_dates (list of str or datetime-like): List of event dates (same length as patient_icd_array).
    
    Returns:
    - None: Displays a plot of the cumulative comorbidities over time.
    
    Example Usage:
    patient_icd_array = ['490.x', '506.0', 'E11.9']
    event_dates = ['2021-01-01', '2021-06-15', '2022-03-20']
    cumulative_comorbidity_plot(patient_icd_array, event_dates)

    
    '''
    
    # Convert event_dates to datetime if they are not already
    if not pd.api.types.is_datetime64_any_dtype(event_dates):
        event_dates = pd.to_datetime(event_dates, errors='coerce')

    assert len(patient_icd_array) == len(event_dates), "The number of ICD-9 code arrays must match the number of event dates."
    
    # Initialize a set to store unique comorbidities
    unique_comorbidities = set()
    cumulative_counts = []

    # Loop through each event
    for i, icd_str in enumerate(patient_icd_array):
        # Extract comorbidities for this event
        current_comorbidities = comorbidity_from_string(icd_str)
        
        # Add new comorbidities to the set of unique comorbidities
        new_comorbidities = current_comorbidities[current_comorbidities == 1].index
        unique_comorbidities.update(new_comorbidities)
        
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
    ''' to convert and localize date columns 
    param s: pandas df, list of strings
    yields pandas df
    example usage:
    df = pd.DataFrame({
    'condition_start_datetime': ['2021-01-01', '2021-02-15', '2021-03-20'],
    'date_of_birth': ['1980-05-10', '1990-06-15', '2000-07-20']
    })

    # Columns to preprocess
    date_columns = ['condition_start_datetime', 'date_of_birth']
    preprocess_dates(df, date_columns)

    
    '''

    for column in date_columns:
        df[column] = pd.to_datetime(df[column], errors='coerce')
        df[column] = df[column].dt.tz_localize(None)
    return df


def calculate_cumulative_comorbidities(person_id, data, comorbidity_function, interval=5):
    '''
    Calculates the cumulative comorbidities for a specific patient over defined intervals.

    Parameters:
    - person_id (int or str): The unique identifier for the patient whose comorbidities are being calculated.
    - data (pd.DataFrame): A DataFrame containing patient data with date and ICD code information. 
      Must include columns: 'person_id', 'condition_start_datetime', and 'source_concept_code'.
    - comorbidity_function (function): A function that takes a string of ICD codes and returns a comorbidity index. 
      This function should handle the conversion of ICD codes to a comorbidity score.
    - interval (int): The number of years for each time interval. Default is 5 years.

    Returns:
    - List of dicts: Each dict contains 'person_id', 'interval_start_year', and 'comorbidity_index' 
      for each time interval.
      
     Example Usage:
      # Sample DataFrame
df = pd.DataFrame({
    'person_id': [1060119, 1060119, 1060119],
    'condition_start_datetime': ['2021-01-01', '2021-06-15', '2022-03-20'],
    'source_concept_code': ['I10', 'E11.9', 'K70.0']
})

# Sample comorbidity function (needs to be defined elsewhere)
def calculate_comorbidity_index(icd_codes_string):
    # Placeholder function - replace with actual implementation
    return len(icd_codes_string.split())

# Calculate cumulative comorbidities for a specific person
results = calculate_cumulative_comorbidities(1060119, df, calculate_comorbidity_index, interval=5)

# Print the results
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

# Function to plot cumulative comorbidities for a specific person
def plot_cumulative_comorbidities_for_person(person_id, data, comorbidity_function, interval=5):
    '''
    Plots the cumulative comorbidity index for a specific patient over defined intervals.

    Parameters:
    - person_id (int or str): The unique identifier for the patient whose comorbidity index is to be plotted.
    - data (pd.DataFrame): A DataFrame containing patient data with date and ICD code information. 
      Must include columns: 'person_id', 'condition_start_datetime', and 'source_concept_code'.
    - comorbidity_function (function): A function that takes a string of ICD codes and returns a comorbidity index. 
      This function should handle the conversion of ICD codes to a comorbidity score.
    - interval (int, optional): The number of years for each time interval. Default is 5 years.

    Returns:
    - None: The function will display a plot of the cumulative comorbidity index over time.
    
    Example Usage: # Sample DataFrame
df = pd.DataFrame({
    'person_id': [1060119, 1060119, 1060119],
    'condition_start_datetime': ['2021-01-01', '2021-06-15', '2022-03-20'],
    'source_concept_code': ['I10', 'E11.9', 'K70.0']
})

# Sample comorbidity function (needs to be defined elsewhere)
def calculate_comorbidity_index(icd_codes_string):
    # Placeholder function - replace with actual implementation
    return len(icd_codes_string.split())

# Plot cumulative comorbidities for a specific person
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
