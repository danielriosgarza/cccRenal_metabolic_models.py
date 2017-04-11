import cobra
import numpy as np
import gurobipy as gu
import scipy.stats as sts 
import math


model1 = cobra.io.read_sbml_model("/home/danielg/PhD/Students/edwin/model_SKRC7.xml")
model2 = cobra.io.read_sbml_model("/home/danielg/PhD/Students/edwin/model_VHL7.xml")

def parse_reaction(model):
    '''
    Parse reaction to dictionary
    input:
        cobra model
    output:
        dict ; dict['reaction'] = ({met: coeff}, (lower bound, upper bound))
    '''    
    result = {}
    r = [key.id for key in model.reactions]; r.sort()
    
    for key in r:
        string = model.reactions.get_by_id(key).reaction
        if '-->' in string:
            bound = (0, 1000)
        elif '<=>' in string:
            bound = (-1000, 1000)
        else:
            print string
        
        met_dict = {}
        m = model.reactions.get_by_id(key).metabolites
        for met in m:
            met_dict[met.id] = m[met]
            
        result[key] = (met_dict, bound)
    return result

def change_signs(dict_for_a_reactions):
    k = dict_for_a_reactions.keys()
    d = dict_for_a_reactions.copy()
    for i in k:
        d[i]=dict_for_a_reactions[i]*-1
    return d

def create_reverse_reactions(reaction_dict):
    d = reaction_dict.copy()
    for i in reaction_dict:
        if reaction_dict[i][1] == (-1000, 1000):
            d[i+'_f']= (reaction_dict[i][0], (0,1000))
            b = change_signs(reaction_dict[i][0])
            d[i+'_r'] = (b, (0,1000))
            d.pop(i)
    return d   
    


def build_stoichiometric_dict(present_reaction_dict):
    '''
    reaction_dict maps reactions to upper and lower bounds
    metab_dict is a dictionary of dictionaries mapping metabolites
    to reactions were they are a part and reactions to the stoichiometric
    coefficient of the metabolite'''
    reaction_dict={}
    metab_dict={}
    for i in present_reaction_dict:
        reaction_dict[i] = present_reaction_dict[i][1]
        metabolites = present_reaction_dict[i][0]
        for z in metabolites.keys():
            if z not in metab_dict:
                metab_dict[z]={}
        for z in metabolites:
            metab_dict[z][i]=metabolites[z] 
    return reaction_dict, metab_dict



def build_gurobi_model(reaction_dict, metabolite_dict, objective_name):
    #define the model
    m= gu.Model('mipl')    
    #define variables
    for i in reaction_dict:
        if i==objective_name:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=1., ub = reaction_dict[i][1], lb = reaction_dict[i][0])
        else:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=0., ub = reaction_dict[i][1], lb = reaction_dict[i][0])
    m.update()
    #set the objective function
    var = m.getVars()
    coef = [1. if i.VarName==objective_name else 0. for i in var]
    m.setObjective(gu.LinExpr(coef, var), gu.GRB.MAXIMIZE)
    m.update()    
    #set the stoichiometric constraints for metabolites    
    for i in metabolite_dict:
        var = metabolite_dict[i].keys()
        var = [m.getVarByName(z) for z in var]
        coef = metabolite_dict[i].values()
        m.addConstr(gu.LinExpr(coef, var), 'E', 0, i)    
    m.update()
    m.optimize()
    return m

def gapfill_gurobi_model(reaction_dict, metabolite_dict, objective_name, delta):
    #define the model
    m= gu.Model('mipl')    
    #define variables
    for i in reaction_dict:
        if i==objective_name:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=delta, ub = reaction_dict[i][1], lb = reaction_dict[i][0])
        else: 
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=-1, ub = reaction_dict[i][1], lb = reaction_dict[i][0])
        
    m.update()
    #set the objective function
    var = m.getVars()
    coef=[]
    for i in var:
        if i.VarName==objective_name:
            coef.append(delta)
        else:
            coef.append(-1)
        
    m.setObjective(gu.LinExpr(coef, var), gu.GRB.MAXIMIZE)
    m.update()    
    #set the stoichiometric constraints for metabolites    
    for i in metabolite_dict:
        var = metabolite_dict[i].keys()
        var = [m.getVarByName(z) for z in var]
        coef = metabolite_dict[i].values()
        m.addConstr(gu.LinExpr(coef, var), 'E', 0, i)    
    m.update()
    m.optimize()
    return m

    
def friday_algorithm(M, objective_name, delta):
    d= delta
    kbase = {i:(M[i][0],(0,1)) for i in M}
    current_sum = 0
    previous_sum = -1
    r, m = build_stoichiometric_dict(kbase)
    modelx = gapfill_gurobi_model(r, m, objective_name, kbase, d)
    sol = modelx.ObjVal
    if sol<0.000001:
        return 'no previous flux--algorithm failled'
    while abs(current_sum - previous_sum)>.0001:
        previous_sum = current_sum
        r, m = build_stoichiometric_dict(kbase)
        modelx = gapfill_gurobi_model(r, m, objective_name, kbase, d)
        sol = modelx.ObjVal
        while modelx.ObjVal>(0.1):
            k = kbase.keys()
            modely = build_gurobi_model(r, m, objective_name) 
            for i in k:
                if modely.getVarByName(i).X==0:
                    kbase.pop(i,None)
            r, m = build_stoichiometric_dict(kbase)
            
            #modely = build_gurobi_model(r, m, objective_name)
            #k = kbase.keys()
            #for i in k:
            #    if modely.getVarByName(i).X==0:
            #        kbase.pop(i,None)
            #r, m = build_stoichiometric_dict(kbase)
            d-=0.1
            print d, '\n\n\n\n\n\n'
            modelx = gapfill_gurobi_model(r, m, objective_name, kbase, d)
         
        
        
        
        
        z = [modelx.getVarByName(i).X for i in k]
    
        
                
        print 'size', '\t', len(kbase), '\n\n'
            
        
        #d = 2*sum(z)+1
            
        current_sum = sum(z)
                        
    return kbase



def tuesday_algorithm(reaction_dict, objective_name):
    '''
    Parameters
    ---------
    reactions_dict: python-dict
    output from parse_reaction(model)
    objective_name: string
    any reaction for which you would like to determine its context specific
    neighbours.
    
    '''
    
    score_dict={i:0 for i in reaction_dict}
    d = sum([reaction_dict[i][1][1] for i in reaction_dict])
    r, m = build_stoichiometric_dict(reaction_dict)
    modelx = gapfill_gurobi_model(r, m, objective_name, d)
    z = sum([modelx.getVarByName(i).X for i in score_dict])
    print z,'\n\n\n'
    for i in score_dict:
        score_dict[i]+=math.log(modelx.getVarByName(i).X+1) -math.log(z+1)
    sol = modelx.ObjVal
    if sol<0.000001:
        return 'no previous flux--algorithm failled'
    while modelx.ObjVal>(0.1):
        d-=0.2*modelx.ObjVal
        modelx = gapfill_gurobi_model(r, m, objective_name, d)
        z = sum([modelx.getVarByName(i).X for i in score_dict])
        for i in score_dict:
            score_dict[i]+=math.log(modelx.getVarByName(i).X+1) -math.log(z+1)
              
    return score_dict


def fix_names(l):
    g=[]
    for i in l:
        if '_f' in i:
            g.append(i[0:-2])
        elif '_r' in i:
            g.append(i[0:-2])
        else:
            g.append(i)
    return list(set(g))

a1 = parse_reaction(model1)
a2 = parse_reaction(model2)


objective = 'HMR_4006'
#a.pop(objective, None)

ac1 = a1.copy()
ac2 = a2.copy()

ac1.pop(objective, None)
ac2.pop(objective, None)


a_r1 = create_reverse_reactions(ac1)
a_r2 = create_reverse_reactions(ac2)


a_r1[objective] = a1[objective]
a_r2[objective] = a2[objective]



a_r1[objective] =(a_r1[objective][0],(0,1))  
a_r2[objective] =(a_r2[objective][0],(0,1))

k1 = a_r1.keys()
k2 = a_r2.keys()
for i in k1:
    if i not in a_r2:
        a_r2[i] = a_r1[i]


for i in k2:
    if i not in a_r1:
        a_r1[i] = a_r2[i]
#
#for i in k2:
#    if i not in a_r1:
#        a_r2.pop(i,None)

y1 = tuesday_algorithm(a_r1, objective)
y2 = tuesday_algorithm(a_r2, objective)


x,y,common=[],[],[]

for i in y1:
    if i in y2:
        x.append(y1[i])
        y.append(y2[i])
        common.append(i)

x,y,common = np.array(x), np.array(y), np.array(common)

from pylab import *
z =abs(x-y)
scatter(x,y);scatter(x[z>0.1], y[z>0.1], c='r')

show()
#model = cobra.io.read_sbml_model("/home/danielg/PhD/Students/Veerle/iRenalCancer1410.xml")
#
#model.change_objective(model.reactions.get_by_id(objective))
#
#model.optimize()
#
#
#k = np.array(y.keys())
#
#l = fix_names(k)
#
#v = np.array(y.values())
#ind = v>-460
#
#t = k[ind]
#
#l = fix_names(t)
#
#for i in model.reactions:
#    if i.id not in l:
#        model.reactions.get_by_id(i.id).upper_bound=0
#        model.reactions.get_by_id(i.id).lower_bound=0
##
#model.optimize()
##
#s = model.solution.f
##
#print 'model_so_far', s
##
#for i in l:
#    mc = model.copy()
#    mc.reactions.get_by_id(i).upper_bound=0
#    mc.reactions.get_by_id(i).lower_bound=0
##
#    mc.optimize()
#    print i, 'original: ', '\t', s, '\t', 'new:', '\t', mc.solution.f, '\t', 'diff:', '\t', s-mc.solution.f
