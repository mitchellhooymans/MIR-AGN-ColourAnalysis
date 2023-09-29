# This script is intended to provide a simple way to analyze the outputs of the MIR-Analysis of AGN Script and
# To output a table of the results, for each field based on a given sigma. The output of this script
# will be a table.

# First begin by importing the relevant libraries
# Begin by importing the required packages for the project
import numpy as np
from tabulate import tabulate
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import os



# Define the inport function
def import_AGN_data(field, all_data, selection_criteria, e_sigma, folder_path='diagnostic selections'):
    """
    Import AGN data based on selection criteria and field.

    This function imports AGN data from CSV files based on the specified selection criteria
    and field. It can retrieve either full data or specific selection criteria data.

    Args:
        field (str): The name of the field for which to import AGN data.
        all_data (bool): If True, import the full AGN data; if False, use the specified selection_criteria.
        selection_criteria (str): The selection criteria used for data retrieval ('Lacy', 'Messias', 'Donley', 'Combined').
        e_sigma (float): The sigma value associated with the data selection.
        folder_path (str, optional): The folder path where CSV files are located (default is 'diagnostic selections').

    Returns:
        pd.DataFrame: A DataFrame containing the imported AGN data based on the specified criteria.

    Note:
        - The function constructs the file name based on the input arguments.
        - It checks if 'all_data' is True to determine whether to import full data or criteria-specific data.
        - If 'all_data' is False, 'selection_criteria' is used to determine which data to import.
        - The function returns 0 and prints a message if an invalid 'selection_criteria' is provided.
    """
    if all_data == True:
        file_name = field+'_full_data_'+str(e_sigma)+'_sigma.csv'
        file_path = os.path.join(folder_path, file_name)
        selection_df = pd.read_csv(file_path)
    else:
        if selection_criteria == 'Lacy':
            file_name = field+'_lacy_selection_'+str(e_sigma)+'_sigma.csv'
            file_path = os.path.join(folder_path, file_name)
            selection_df = pd.read_csv(file_path)
        elif selection_criteria == 'Messias':
            file_name = field+'_messias_selection_'+str(e_sigma)+'_sigma.csv'
            file_path = os.path.join(folder_path, file_name)
            selection_df = pd.read_csv(file_path)
        elif selection_criteria == 'Donley':
            file_name = field+'_donley_selection_'+str(e_sigma)+'_sigma.csv'
            file_path = os.path.join(folder_path, file_name)
            selection_df = pd.read_csv(file_path)
        elif selection_criteria == 'Combined':
            file_name = field+'_combined_selection_'+str(e_sigma)+'_sigma.csv'
            file_path = os.path.join(folder_path, file_name)
            selection_df = pd.read_csv(file_path)
        else:
            print('Please enter a valid selection criteria')
            return 0
    return selection_df


def import_truth_sample(field):
    """
    Import and preprocess truth samples for a given field.

    This function reads in truth samples from FITS files, creates a DataFrame,
    and performs data preprocessing to filter and calculate AGN-related properties.

    Args:
        field (str): The name of the field for which to import truth samples.

    Returns:
        pd.DataFrame: A DataFrame containing preprocessed truth samples for the field.

    Note:
        - The function assumes the presence of FITS files with names following the pattern '{field}_truth.fits'.
        - It calculates AGN contribution and adds a 'Known AGN' column based on predefined thresholds.
    """
    # Read in the truth samples and then create dataframes from the fits files
    file_path =os.path.join("./truth sample/"+field+'_truth.fits')
    truth_sample = fits.open(file_path)
    truth_df = pd.DataFrame(np.array(truth_sample[1].data).byteswap(
    ).newbyteorder())  # Byteswap so that Pandas can read it
    # Rename the ID column so that it matches the other dataframes
    truth_df.rename(columns={'id_1': 'id'}, inplace=True)
    # bayes.agn.total_dust_luminosity <- AGN contribution
    truth_df['bayes.agn.total_dust_luminosity']

    # bayes.dust.luminosity <- stellar contribution
    truth_df['bayes.dust.luminosity']

    # filter the truth data based on error
    truth_df = truth_df[truth_df['bayes.agn.total_dust_luminosity'] > 0]
    truth_df = truth_df[truth_df['bayes.dust.luminosity'] > 0]

    # This will be the AGN luminosity contribution
    truth_df['agn contribution'] = truth_df['bayes.agn.total_dust_luminosity'] / \
        (truth_df['bayes.agn.total_dust_luminosity'] +
         truth_df['bayes.dust.luminosity'])

    # Add a new column for known AGN
    truth_df['Known AGN'] = np.where(truth_df['agn contribution'] > 0.5, 1, 0)

    # This will be the AGN luminosity contribution
    # num_true_AGN = len(truth_df[truth_df['Known AGN'] == 1])
    return truth_df

def calculateCompleteness(df, diagnostic):
    # Calculate the completeness of the selection
    # Completeness = Positive Diagnostic / Known AGN
    # Positive Diagnostic = AGN that are selected by CIGALE, and as AGN by the selection diagnostic
    # Known AGN = AGN selected by CIGALE
    positive_selection = df['Positive '+ diagnostic + " Selection"].value_counts()[1]
    known_AGN = df['Known AGN'].value_counts()[1]
    return positive_selection/known_AGN

def calculateReliability(df, diagnostic):
    # Calculate the reliability of the selection
    # Reliability = Positive Diagnostic / Positive Selection
    # Positive Diagnostic = AGN that are selected by CIGALE, and as AGN by the selection diagnostic
    # Positive Selection = AGN that are selected by the selection diagnostic
    positive_selection = df['Positive '+ diagnostic + " Selection"].value_counts()[1]
    diagnostic_selection = df[diagnostic + ' Selection'].value_counts()[1]
    return positive_selection/diagnostic_selection




# Import the truth sample for each field here to save time in the loop below
cdfs_truth = import_truth_sample('CDFS')
cosmos_truth = import_truth_sample('COSMOS')
uds_truth = import_truth_sample('UDS')


fields = {"CDFS": cdfs_truth, "COSMOS": cosmos_truth, "UDS": uds_truth}
# Read in the relevant data



# As we are trying to compare each field in a given sigma range, we need to read in data based on a sigma value
for sigma in range(1, 6):
    results = []
    # Read in the data for each field, and then proccess
    for field in fields:
        # Create a new dictionary to store results for this field
        field_results = {'Field': field}
        
        # Import the data 
        selections_df = import_AGN_data(
            field, False, 'Combined', sigma, 'diagnostic selections')
        selections_df = selections_df.join(
            fields[field].set_index('id'), on='id')
        
        # Create a new column for positive diagnostic, this will be set to zero initally
        selections_df['Positive Lacy Selection'] = 0
        selections_df['Positice Messias Selection'] = 0
        selections_df['Positive Donley Selection'] = 0

        # We now need to test against the AGN to see if we have a positive selection
        # First for Lacy, Donley, then Messias
        selections_df['Positive Lacy Selection'] = np.where(
            (selections_df['Known AGN'] == 1) & (selections_df['Lacy Selection'] == 1), 1, 0)
        selections_df['Positive Donley Selection'] = np.where(
            (selections_df['Known AGN'] == 1) & (selections_df['Donley Selection'] == 1), 1, 0)
        selections_df['Positive Messias Selection'] = np.where(
            (selections_df['Known AGN'] == 1) & (selections_df['Messias Selection'] == 1), 1, 0)
        
        # Calculate the completeness of each selection
        lacy_completeness = calculateCompleteness(selections_df, 'Lacy')
        messias_completeness = calculateCompleteness(selections_df, 'Messias')
        donley_completeness = calculateCompleteness(selections_df, 'Donley')
        
        # Calculate the reliability of each selection
        lacy_reliability = calculateReliability(selections_df, 'Lacy')
        messias_reliability = calculateReliability(selections_df, 'Messias')
        donley_reliability = calculateReliability(selections_df, 'Donley')
        
        # Add completeness and reliability results to the dictionary
        field_results['Lacy Completeness (%)'] = (round(lacy_completeness*100, 2))
        field_results['Lacy Reliability (%)'] = (round(lacy_reliability*100, 2))
        field_results['Donley Completeness (%)'] = (round(donley_completeness*100, 2))
        field_results['Donley Reliability (%)'] = (round(donley_reliability*100, 2))
        field_results['Messias Completeness (%)'] = (round(messias_completeness*100, 2))
        field_results['Messias Reliability (%)'] = (round(messias_reliability*100, 2))
        
        
        
        # This is where the magic will need to happen
        # In addition to this, if we are in the CDFS field we can also explore the Szokoly XRay AGN selection criteria
        if field == "CDFS":
            print("XRAY WORKING AREA\n\n\n\n\n")
            # Import the data 
            xfolder_path='diagnostic selections'
            xfile_name = 'CDFS_xagn_selection_'+str(sigma)+'_sigma.csv'
            xfile_path = os.path.join(xfolder_path, xfile_name)
            xray_df = pd.read_csv(xfile_path)
            
            #print("Elements of xray df: "+str(xray_df['id'].head()))
            #print("Elements of cdfs df: "+str(fields[field]['id'].head()))
        
            print("Xray df shape: "+str(xray_df['id']))
            print("CDFS df shape: "+str(fields[field]['id']))
            xray_df = xray_df.join(fields[field].set_index('id'), on='id')
            
            
            


            #print("Known AGNs in the CDFS field:" + str(xray_df['Known AGN'].value_counts()))
            
            # AGN Diagnostics using the truth sample
            print("Known AGNs in the CDFS field: " + str(xray_df['Known AGN'].value_counts()[1]))
            print("Known unknown in the CDFS field: " + str(xray_df['Known AGN'].value_counts()[0]))
            # Selected by Szokoly
            print("Selected by Szokoly: " + str(xray_df['Szokoly Selection'].value_counts()[1]))   
            print("Not selected by Szokoly: " + str(xray_df['Szokoly Selection'].value_counts()[0]))        
            
            
            # # Create a new column for positive diagnostic, this will be set to zero initally
            xray_df['Positive Szokoly Selection'] = 0
            
            # Need to make a positive selection for when their are both known AGN (as per the truth sample) and a positive Szokoly selection
            xray_df.loc[(xray_df['Known AGN'] == 1) & (xray_df['Szokoly Selection'] == 1), 'Positive Szokoly Selection'] = 1
            
            
            
            # Present in both Szoloky and the truth sample
            #print("Present in both Szoloky and the truth sample: " + str(xray_df['Positive Szokoly Selection'].value_counts()[1]))
            
            # # # Calculate the completeness of each selection
            # szokoly_completeness = calculateCompleteness(xray_df, 'Szokoly')
            
            # # # Calculate the reliability of each selection
            # szokoly_reliability = calculateReliability(xray_df, 'Szokoly')
            
            # # # Add completeness and reliability results to the dictionary
            # field_results['Szokoly Completeness (%)'] = (round(szokoly_completeness*100, 2))
            # field_results['Szokoly Reliability (%)'] = (round(szokoly_reliability*100, 2))

        # Create a dataframe from the results list
        results_df = pd.DataFrame(results)
        print(str(sigma)+ "sigma")
        print(results_df)
        # Save the results to a CSV file
        folder_path = 'diagnostic comparisons'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        results.append(field_results)
    
    
    # Seperately we also have XRay AGN data. Unfortuantely this only exsists for the CDFS field
    # as such we can edit the diagnostic comparisons above for each sigma value
    # and add an extra column for the Szokoly XRay AGN Selection and then calculate the completeness and reliability
    # remembering that the CDFS field is the only field with XRay AGN data. 
    # the file selected has the name CDFS_xagn_selection_n_sigma.csv where the n is the sigma number
    print(results)
    
    
    
    
    
    
    
        
    # Create a dataframe from the results list
    results_df = pd.DataFrame(results)
    print(str(sigma)+ "sigma")
    print(results_df)
    # Save the results to a CSV file
    folder_path = 'diagnostic comparisons'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = 'diagnostic_comparisons_'+str(sigma)+'_sigma.csv'
    file_path = os.path.join(folder_path, file_name)
    results_df.to_csv(file_path, index=False)
        
