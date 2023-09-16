## This script is intended to provide a simple way to analyze the outputs of the MIR-Analysis of AGN Script and 
## To output a table of the results, for each field based on a given sigma. The output of this script
## will be a table.

# First begin by importing the relevant libraries
# Begin by importing the required packages for the project
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import os


# Define the inport function 
def import_AGN_data(field, all_data, selection_criteria, e_sigma, folder_path = 'diagnostic selections'):
    
    if all_data == True:
        file_name =  field+'_full_data_'+str(e_sigma)+'_sigma.csv'
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


fields = ["CDFS", "COSMOS", "UDS"]
# Read in the relevant data


## As we are trying to compare each field in a given sigma range, we need to read in data based on a sigma value
for sigma in [1, 2, 3, 4, 5]:
    for field in fields:

        # Read in the truth samples and then create dataframes from the fits files
        truth_sample = fits.open(field+'_truth.fits')
        truth_df=pd.DataFrame(np.array(truth_sample[1].data).byteswap().newbyteorder()) # Byteswap so that Pandas can read it
        truth_df.rename(columns={'id_1':'id'}, inplace=True) # Rename the ID column so that it matches the other dataframes
        # bayes.agn.total_dust_luminosity <- AGN contribution
        truth_df['bayes.agn.total_dust_luminosity']

        # bayes.dust.luminosity <- stellar contribution
        truth_df['bayes.dust.luminosity']

        # filter the truth data based on error
        truth_df = truth_df[truth_df['bayes.agn.total_dust_luminosity'] > 0]
        truth_df = truth_df[truth_df['bayes.dust.luminosity'] > 0]

        # This will be the AGN luminosity contribution
        truth_df['agn contribution'] = truth_df['bayes.agn.total_dust_luminosity']/(truth_df['bayes.agn.total_dust_luminosity'] + truth_df['bayes.dust.luminosity'])
    
        # Add a new column for known AGN
        truth_df['Known AGN'] = np.where(truth_df['agn contribution'] > 0.5, 1, 0)

        # This will be the AGN luminosity contribution
        num_true_AGN = len(truth_df[truth_df['Known AGN'] == 1])


