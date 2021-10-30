#import sympy as sp
#import gravipy as gp
import numpy as np
import itertools
import tableau_utils as tu
from cycle import Cycle
import fractions as fr
from cycle_sum import CycleSum

class Symmetriser(CycleSum):

    def __init__(self, num_list, antisym = False, normalise = False):

        perms = list(itertools.permutations(num_list))
        self.domain = num_list

        odd = 0
        if antisym:
            odd = 1

        operator = []
        prefactors = []
        num_list.sort()
        for perm_ind in range(len(perms)):

            op = Cycle()
            op.set_cycle({'image': [x for x in perms[perm_ind]], 'domain': num_list})

            op.remove_single_cycles()
            perm_order = op.permutation_order()

            del_i = int((-1)**(odd*perm_order))
            prefactors.append(del_i)
            operator.append(op)
            

        self.operator = CycleSum(operator, prefactors, normalise = False)

        self.antisym = antisym
        
        if normalise:
            self.normalise()
        
        
        
    def get_total_domain(self):
        
        return self.operator.get_total_domain()
        
        
    def get(self, prpty):
        
        if prpty is 'domain':
            return self.get_total_domain()
        else: 
            return self.operator.get(prpty)
        
    def get_sum(self):
        
        return self.operator.get_sum()
    
    def get_operator(self):
        
        return self.operator
    
    
    def is_antisym(self):
        
        return self.antisym
    
    def is_normalised(self):
        
        return self.operator.is_normalised()

    def set_normalised(self, val):
        
        self.operator.set_normalised(val)
    
    
    def set_operator(self, op):
        
        self.operator = op
        
        
        
    def set_symmetriser_properties(self, op, antisym = False):
        
        self.operator = op
        self.antisym = antisym
        
        
    def reverse(self):
        
        self.operator = self.operator.reverse()
        return self
        

    def normalise(self):
    
        if not self.is_normalised():
            
            n = len(self.domain)
            self.operator = self.operator.multiply_by_constant(fr.Fraction(1,np.prod(range(1, n+1))))
        
            self.set_normalised(True)
    

    
    def act_on(self, op):
        
        op = op.get_operator() # cyc sum is always saved under 'operator'
        
        return self.operator.act_on(op)
    
    
    def is_subsym_of(self, op):
        
        me_antisym = self.antisym
        op_antisym = op.is_antisym()
        
        my_cycles = [self.operator[ind]['cycle'] for ind in range(len(self.operator))]
        op_cycles = []
        if 'Sym' in str(type(op)):
            op_cycles = op.get_operator().get('cycle')
        else:
            op_cycles = op.get('cycle')
        
        subsym = []
        is_subsym = False
        
        if me_antisym is op_antisym:
            
            for my_cyc in my_cycles:
                
                cyc_sub = False
                for op_cyc in op_cycles:
                    
                    if op_cyc.equivalent(my_cyc):
                        cyc_sub = True
                        
                if cyc_sub:
                    subsym.append(1)
                    
        if len(subsym) == len(my_cycles):
            is_subsym = True
            
        return is_subsym
        
    
    
    def write_as_deltas(self, left_label = 'x', right_label = 'y'):
    
        return self.operator.write_as_deltas(left_label = 'x', right_label = 'y')



    def write_as_cycles(self):
        
        return self.operator.write_as_cycles()
    
    def trace(self, N = 1, num_lines = None):
        
        return self.operator.trace(N = 1, num_lines = None)
