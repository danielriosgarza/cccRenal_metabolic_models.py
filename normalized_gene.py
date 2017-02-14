# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 05:46:22 2017

@author: edwina
"""

import pandas as pd
import cobra
import numpy as np
import quantile_normalization as qn

'''
This script:
1. Remove all genes not in the model
2. Remove missing values
3. Quantile normalize the gene expressions
4. Take median of expressin values for each reaction
'''
# Open cell line DB and add blank reactionID
Celldb = pd.read_excel('1_Cell_line_db.xlsx', sheetname =0)

for i in range(len(Celldb)):
    if pd.isnull(Celldb['Reaction ID'][i]):
        Celldb['Reaction ID'][i] = Celldb['Reaction ID'][i-1]
        
# Make a set of all reactions in original model
model = cobra.io.read_sbml_model('0_iRenalCancer1410.xml')

reaction_in_model = set()
for reaction in model.reactions:
    reaction_in_model.add(str(reaction))
    
# Subset reaction only in the model
Celldb_model = Celldb[Celldb['Reaction ID'].isin(reaction_in_model)]

# Remove missing values
Celldb_model = Celldb_model[pd.notnull(Celldb_model['SKRC-7'])].reset_index()

# Make numpy array of SKRC7 and VHL7 expressions
skrc7 = Celldb_model['SKRC-7'].tolist()
vhl7 = Celldb_model['VHL-7'].tolist()
exp_array = np.array([skrc7, vhl7])

# Perform quantile normalization
n_array = qn.quantile_normalization(exp_array)

#Make dataframe with normalized values
dataframe = pd.DataFrame(n_array.T, columns = ['SKRC-7', 'VHL-7'])
dataframe = pd.concat([Celldb_model['Reaction ID'], dataframe], axis = 1)

# Take median values for each reaction
grouped = dataframe.groupby('Reaction ID', as_index = False)
normalized = grouped.median()

# Save to file\
normalized.to_csv("normalized_expression.csv", sep = '\t', index = False)






