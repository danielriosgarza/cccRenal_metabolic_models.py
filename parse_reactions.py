
import cobra

model = cobra.io.read_sbml_model('0_iRenalCancer1410.xml')


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
    
a = parse_reaction(model)
