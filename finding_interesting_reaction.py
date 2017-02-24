#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 05:30:40 2017

@author: edwina
"""

import cobra
import re
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import quantile_normalization as qn
from pylab import *
import scipy.stats as sts


## Make function to update upper and lower bound
def update_bound(model, exp_dict, exp_key, cell_line):
    '''
    Update model upper and lower bound based on expression values.
    Input:
        model    : cobra model
        exp_dict : dict
        exp_key  : list
        cell_line: int ; skrc7 = 0, vhl7 = 1 
    Output:
        model: cobra model
    '''
    md = model.copy()
    for key in exp_key:
        b = exp_dict[key][cell_line]
        #if b>1000:
        #    b=1000
        if re.search('<=>', model.reactions.get_by_id(key).reaction):
            md.reactions.get_by_id(key).upper_bound = b
            md.reactions.get_by_id(key).lower_bound = -b
        elif re.search('-->', model.reactions.get_by_id(key).reaction):
            md.reactions.get_by_id(key).upper_bound = b
    return md



'''          MODEL USING NORMALIZED MEDIAN EXPRESSION '''
    
# Import the normalized expression values for skrc7 nad vhl7
exp = open("14_normalized_expression.csv", "r")
header = exp.readline().strip().split('\t')

# Make dictionary of normalized expression value (in list [skrc7, vhl7])
exp_dict = {}
for line in exp:
    line = line.strip().split('\t')
    exp_dict[line[0]] = [float(line[1]), float(line[2])]
    
exp_key = exp_dict.keys()

v1 = [exp_dict[i][0] for i in exp_key]
v2 = [exp_dict[i][1] for i in exp_key]

mean1 = mean(v1)
mean2 = mean(v2)



    
#  Import the model
model_orig = cobra.io.read_sbml_model("0_iRenalCancer1410.xml")

'''Model 1
Set the upper and lower bound to normalized expressions'''
result_dict = {}



'''MODEL 1
Set upper and lower bound to expression values'''
m1_skrc7 = model_orig.copy()
m1_vhl7  = model_orig.copy()
# Update model upper and lower bound 
m1_skrc7=update_bound(m1_skrc7, exp_dict, exp_key, 0)
m1_vhl7 = update_bound(m1_vhl7, exp_dict, exp_key, 1)
# Set objective function and optimize the model
m1_skrc7.objective = ['CancerBiomass_OF']
m1_vhl7.objective = ['CancerBiomass_OF']

m1_skrc7.optimize()
m1_vhl7.optimize()


reaction_list = m1_skrc7.solution.x_dict.keys()

for i in reaction_list:
    if i not in exp_key:
        if '<=>' in model_orig.reactions.get_by_id(i).reaction:
            m1_skrc7.reactions.get_by_id(i).upper_bound = mean1
            m1_skrc7.reactions.get_by_id(i).lower_bound = -mean1
            m1_vhl7.reactions.get_by_id(i).upper_bound = mean2
            m1_vhl7.reactions.get_by_id(i).lower_bound = -mean2
        elif '-->' in model_orig.reactions.get_by_id(i).reaction:
            m1_skrc7.reactions.get_by_id(i).upper_bound = mean2
            m1_vhl7.reactions.get_by_id(i).upper_bound = mean2


zero_keys_1=[]
zero_keys_2=[]

for i in exp_key:
    if exp_dict[i][0]==0:
        if '<=>' in model_orig.reactions.get_by_id(i).reaction:
            zero_keys_1.append(i)
        else:
            zero_keys_2.append(i)


m1_skrc7.optimize()
m1_vhl7.optimize()
###############################################################################
'''                    FINDING INTERESTING REACTIONS
MODEL 1'''
# Make list of all reactions in model
reaction_list = m1_skrc7.solution.x_dict.keys()
# Make list of all flux in model 1
flux1_skrc7 = [m1_skrc7.solution.x_dict[key] for key in reaction_list]
flux1_vhl7  = [m1_vhl7.solution.x_dict[key] for key in reaction_list]



# Take ratio

m1_s_newlist = np.array(flux1_skrc7)

m1_v_newlist = np.array(flux1_vhl7)


m1_s_newlist = m1_s_newlist/np.linalg.norm(m1_s_newlist) 
          
m1_v_newlist = m1_v_newlist/np.linalg.norm(m1_v_newlist) 
            
ratio_m1_new = abs(m1_s_newlist - m1_v_newlist)

ratio_m1_new = ratio_m1_new

#ind3 = ratio_m1_new<200 

#ratio_m1_new = ratio_m1_new[ind3]

b = 1*std(ratio_m1_new)
ind = ratio_m1_new>b
reaction_array = np.array(reaction_list)
#reaction_array = reaction_array[ind3]
m1_interesting = reaction_array[ind]

# Save the interseting reaction to file
int_file = 'Reaction ID\n'
for i in xrange(len(m1_interesting)):
    int_file = int_file + str(m1_interesting[i]) + '\n'
f_int = open("15_interesting_reactions_2.csv", "w")
f_int.write(int_file)
f_int.close()

