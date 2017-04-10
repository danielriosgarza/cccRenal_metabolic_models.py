import cobra
import numpy as np
import gurobipy as gu
import scipy.stats as sts 



model = cobra.io.read_sbml_model("/home/danielg/PhD/Students/Veerle/iRenalCancer1410.xml")

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

def gapfill_gurobi_model(reaction_dict, metabolite_dict, objective_name, candidate_set_dict, delta):
    #define the model
    m= gu.Model('mipl')    
    #define variables
    for i in reaction_dict:
        if i==objective_name:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=delta, ub = reaction_dict[i][1], lb = reaction_dict[i][0])
        elif i in candidate_set_dict:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=-1, ub = reaction_dict[i][1], lb = reaction_dict[i][0])
        else:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, obj=0, ub = reaction_dict[i][1], lb = reaction_dict[i][0])
    m.update()
    #set the objective function
    var = m.getVars()
    coef=[]
    for i in var:
        if i.VarName==objective_name:
            coef.append(delta)
        elif i.VarName in candidate_set_dict:
            coef.append(-1)
        else:
            coef.append(0)

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

    
def friday_algorithm(kbase_dict, M, objective_name, delta):
    d= delta
    kbase = {i:(kbase_dict[i][0],(0,1)) for i in kbase_dict}
    current_sum = 0
    previous_sum = -1
    while abs(current_sum - previous_sum)>.0001:
        previous_sum = current_sum
        r, m = build_stoichiometric_dict(kbase)
        modelx = gapfill_gurobi_model(r, m, objective_name, M, d)
        
            
               
        if modelx.getVarByName(objective_name).X<0.1:
            pass
        else:
            print len(kbase), '\n\n'
            z = [modelx.getVarByName(i).X for i in M]
            d = 2*sum(z)
            
            current_sum = sum(z)
            k = M.keys() 
                
            for i in k:
                if modelx.getVarByName(i).X==0:
                    kbase.pop(i,None)
                    M.pop(i,None)
            print len(kbase), '\n\n'
            
    return kbase


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

a = parse_reaction(model)
objective = 'HMR_2771'
#a.pop(objective, None)
ac = a.copy()
ac.pop(objective, None)

a_r = create_reverse_reactions(ac)

a_r[objective] = a[objective]

a_r[objective] =(a_r[objective][0],(0,1))  

r, m = build_stoichiometric_dict(a_r)
modelx = gapfill_gurobi_model(r, m, objective, a_r, 1000)


c=[]
for i in a_r:
    if modelx.getVarByName(i).X!=0:
        c.append(i)

n=[]
for i in a_r:
    if modelx.getVarByName(i).X==0:
        n.append(i)

model = cobra.io.read_sbml_model("/home/danielg/PhD/Students/Veerle/iRenalCancer1410.xml")
model.change_objective(model.reactions.HMR_2771)
model.optimize()

g = fix_names(c)

mc = model.copy()

l =[] 

for i in a_r:
    if i not in c:
        l.append(i)

l_f = fix_names(l)

for i in l_f:
    model.reactions.get_by_id(i).upper_bound=0
    model.reactions.get_by_id(i).lower_bound=0
model.optimize()

print model.solution.f
