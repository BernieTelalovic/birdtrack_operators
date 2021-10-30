#import sympy as sp
#import gravipy as gp
import numpy as np
import itertools
import yaml
import fractions as fr
from cycle import Cycle
import copy as cp

def read_yaml(conf_path):
    
    my_dict = None
    with open(conf_path, 'r') as stream:
        try:
            my_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return my_dict


def permutation_dict(partition):
     return {str(part): len([pr for pr in partition if pr == part]) for part in partition}

def num_elements_in_conjugacy_class(partition, N):
    
    if np.sum(partition) <= N:
        
        parts = {str(part): len([pr for pr in partition if pr == part]) for part in partition}
        denom = 1

        for pr in parts.keys():
            part = int(pr)
            denom = denom*(part**parts[pr])*np.math.factorial(parts[pr])

        num_elements = fr.Fraction(np.math.factorial(N)/denom)
    
    else:
        num_elements = None
        print('The partition numbers exceed the group N; returning None')
    return num_elements

def cycle_convert(perm):
    
    cyc = str([list(np.array(pr)+1) for pr in perm.cyclic_form]).replace('[]', '')\
                         .replace('[[', '(').replace(']]', ')')\
                         .replace('[','(').replace(']', ')').replace(', ', ',')
    return Cycle(cyc)

def elements_of_S(N):
    
    import sympy.combinatorics as spcom
    perms = spcom.generators.symmetric(N)
    
    return [cycle_convert(perm) for perm in perms]
    

def conjugacy_class(N,partition):
    
    num_elements = num_elements_in_conjugacy_class(partition, N)
    
    elements = []
    if num_elements is not None:
        
        perms = elements_of_S(N)
        for perm in perms:
            pr = [per.replace('(', '').replace(')', '').split(' ') for per in cp.deepcopy(perm.as_str()).split(')(')]
            perm_lengths = [len(per) for per in pr]
            
            if permutation_dict(perm_lengths) == permutation_dict([part for part in partition if part > 1]):
                elements.append(perm)
    return elements


def get_type(strin):
    
    return strin.replace("'>", '').split('.')[-1]



def multiply_cycles(a,b):

    #a o b
    cycle = {'image': [], 'domain': []}

    for entry_b in b['domain']:

        b_ind = b['domain'].index(entry_b)
        b_perm = b['image'][b_ind]

        if b_perm in a['domain']:
            a_ind = a['domain'].index(b_perm)
            a_perm = a['image'][a_ind]

            cycle['image'].append(a_perm)
            cycle['domain'].append(entry_b)

        else:
            cycle['image'].append(b_perm)
            cycle['domain'].append(entry_b)


    for entry_a in a['domain']:

        if not (entry_a in cycle['domain']):
            cycle['domain'].append(entry_a)
            perm_a = a['image'][a['domain'].index(entry_a)]
            cycle['image'].append(perm_a)

    cycle['image'] = [x for y,x in sorted(zip(cycle['domain'],cycle['image']))]
    cycle['domain'] = np.sort(cycle['domain'])

    return remove_single_cycles([cycle])[0]



def remove_single_cycles(ops):
    
    active_cycles = []
    
    for op_ind in range(len(ops)):
        
        left = ops[op_ind]['image']
        right = ops[op_ind]['domain']
        
        active_left = [left[ind] for ind in range(len(right)) if left[ind]!=right[ind]]
        active_right = [right[ind] for ind in range(len(right)) if right[ind]!=left[ind]]
        
        active_cycles.append({'image': active_left, 'domain': active_right,
                              'prefactor': ops[op_ind].get('prefactor', 0)})

    return active_cycles


def simplify_operator(ops):
    
    simplified = []
    
    for op in ops:

        appended = False
        for sim_ind in range(len(simplified)):
                
            sim = simplified[sim_ind]
                
            if (op['domain'] == sim['domain']) and (op['image'] == sim['image']):
                    
                simplified[sim_ind]['prefactor'] += op['prefactor']
                appended = True
                    
        if not appended:
                simplified.append(op)
                
    return simplified


def cycle_num(perm):
    
    left = perm['image']
    right = perm['domain']
    cycles_complete = 0
    
    if len(right) > 0:
        ind = 0
        ri_start = right[ind]
        ri = ri_start
        count = 0
        
        cycle = []

        indices = right.copy()

        while count < len(right):
    
            lf = left[ind]
            indices.remove(lf)
            cycle.append(lf)
        
            if lf == ri_start:
                
                if len(cycle) > 2:
                    cycles_complete += len(cycle) - 1
                else:
                    cycles_complete +=1
                ri_start = 0
                cycle = []
            
                if len(indices) > 0:
                    ri_start = indices[0]
                    ind = right.index(ri_start)

                else: 
                    count = len(right)
        
            else:
                ind = right.index(lf)
                ri = right[ind]

            count += 1
    
    return cycles_complete
            

def permutation_order(perm):
    
    perm_num = cycle_num(perm)
    
    if perm_num > 0:
        perm_num += 0
        
    return perm_num


def split_by_cycle(lis):
    
    rtn_str = ''
    
    if len(lis)> 0:
        list_sort = sorted(lis, key=len)
        len_fir = len(list_sort[0])
        
        for ent in list_sort:
            if len(ent) > len_fir:
                len_fir = len(ent)
                rtn_str = rtn_str + '\n ' + ent + ' + '
            else:
                rtn_str = rtn_str + ent + ' + '
    
    return rtn_str
    
    

def split_by_factor(strin):
    
    sum_str = strin.replace(' ', ',')
    sum_str = sum_str.replace('+', ' +').replace('-', ' -')
    
    sum_str = sum_str.split(' ')
    
    cycles = {pref.split('(')[0]: [] for pref in sum_str}
    
    for cyc in sum_str:
        
        split_cyc = cyc.split('(')
        cyc_str = ''
        for ind in range(len(split_cyc)):
            if ind > 0:
                cyc_str = cyc_str + '('+ split_cyc[ind]
        
        cycles[cyc.split('(')[0]].append(cyc_str)
        
    rtn_str = ''
    for ky in cycles.keys():
        rtn_str += ky + '[' + split_by_cycle(cycles[ky]) + ']' + '\n\n'
    return rtn_str
