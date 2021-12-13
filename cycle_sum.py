import sympy as sp
import gravipy as gp
import numpy as np
import itertools
from cycle import Cycle
import copy as cp
import fractions as fr

class CycleSum():
    
    def __init__(self, cycle_list, prefactors = [], str_prefactors = [], normalise = False):
    
        cyc_lis = []
        prefacs = []
        self.normalised = False
        for cycle_ind in range(len(cycle_list)):
            
            if (len(cycle_list) == len(prefactors)):
                prefacs.append(prefactors[cycle_ind])
                if prefactors[cycle_ind] != 0:
                    if len(cycle_list) == len(str_prefactors):
                        cyc_lis.append({'cycle': cycle_list[cycle_ind], 'prefactor': fr.Fraction(str(prefactors[cycle_ind])),
                                    'prefactor_str': str_prefactors[cycle_ind]})
                    else:
                        cyc_lis.append({'cycle': cycle_list[cycle_ind], 'prefactor': fr.Fraction(str(prefactors[cycle_ind])),
                                    'prefactor_str': str(prefactors[cycle_ind])})
            else:
                cyc_lis.append({'cycle': cycle_list[cycle_ind], 'prefactor': fr.Fraction(1,1),
                                'prefactor_str': str(1)})
                prefacs.append(1)
                
        self.cycle_sum = cyc_lis
        self.cycle_list = cycle_list
        if normalise:
            normalised = self.normalise()
            prefacs = normalised.prefactors
            self.cycle_sum = normalised.get_sum()
        
        self.prefactors = prefacs
        
        
    def CycleSumFromStr(cyc_str):
        
        cyc_str = cyc_str.replace('+', '`+')
        cyc_str = cyc_str.replace('-', '`-')
        cyc_strlis = cyc_str.split('`')
        start = '('
        end = ')'
        prefs = []
        cycles = []

        for s in cyc_strlis:
            if len(s) > 0:
                cyc = '('+ s[s.find(start)+len(start):s.rfind(end)] + ')'
                pref = s.split(cyc)[0]
                if cyc in '()':
                    cyc = ''
                cycles.append(Cycle(cyc))
                prefs.append(pref)
        cycsum = CycleSum(cycles, prefactors = prefs)
        
        return cycsum
        
    def is_normalised(self):
        
        return self.normalised
    
        
    def get_operator(self):
        
        return self
        
        
    def get_sum(self):
        
        return self.cycle_sum
    
    
    def get_total_domain(self):
    
        all_cycs = self.get('cycle')
        doms = []
        for cyc in all_cycs:
            dom = cyc.get('domain')
            for ent in dom:
                doms.append(ent)
        
        return list(set(doms))
    
    
    def get_identity_factor(self):
        
        id_fac = None
        for cyc in self.cycle_sum:
            if cyc['cycle'].as_str() in '()':
                id_fac = cyc['prefactor']
                
        if id_fac is None:
            print('WARNING: id not found in this element, returning the overall factor instead.')
            id_fac = self.get_overall_factor()
            
        return id_fac
    
    
    def get_overall_factor(self):
        
        overall_fac = None
        
        prefs = self.get('prefactor')
        denoms = [pref.denominator for pref in prefs]
        numers = []
        denom_fac = None
        numer_fac = None
        if len(denoms) > 0:
            denom_lcm = np.lcm.reduce(denoms)
           # print(denom_lcm)
            numers = [int(prefs[ind].numerator*(denom_lcm/prefs[ind].denominator)) for ind in range(len(prefs))]
           # print(numers)
           # print(denoms)
           # print()
        
        if len(numers) > 0:
            numer_fac = np.gcd.reduce(numers)
            denom_fac = denom_lcm
            overall_fac = fr.Fraction(numer_fac,denom_lcm)
            
        return overall_fac
        
        
    def get(self, prpty):
        
        if prpty is 'domain':
            return self.get_total_domain()
        else:
            return [cyc.get(prpty, None) for cyc in self.cycle_sum]
    
    
    def set_cycle_sum(self, val):
        
        self.cycle_sum = val
        
    def set_normalised(self, val):
        
        self.normalised = val


    def sum_with(self, cycles, factor = 1):

        cycles_fac = cycles.multiply_by_constant(factor)

        all_cycs = self.cycle_sum + cycles_fac.cycle_sum

        cy = CycleSum([Cycle()])
        cy.set_cycle_sum(all_cycs)

        sums = cy.simplify()

        return sums

            
    
    def sum_with_depricated(self, cycles, factor = 1):
        
        added = cp.deepcopy(self).get('cycle')
        prefactors = cp.deepcopy(self).get('prefactor')

        added += cycles.get('cycle')
        prefactors += [fr.Fraction(str(factor))*fr.Fraction(str(cy)) for cy in cycles.get('prefactor')]

        added_cycs = CycleSum(added, prefactors = prefactors).simplify()
        
        return added_cycs
    
    
    def sum_with_list(self, cycles, factor = 1):
        
        added = cp.deepcopy(self)
        
        for cyc in cycles:
            
            added = added.sum_with(cyc, factor = factor)
            
        return added

    
    def multiply_by_constant(self, const):
        
        cycles = self.get('cycle')
        prefacs = self.get('prefactor')
        
        prefactors = [fr.Fraction(str(const))*fr.Fraction(str(pref)) for pref in prefacs]
        
        return CycleSum(cycles, prefactors = prefactors)
    
    
    def reverse(self):

        for cyc_ind in range(len(self.cycle_sum)):
            self.cycle_sum[cyc_ind]['cycle'] = self.cycle_sum[cyc_ind]['cycle'].reverse()            
        return self
    
    
        
    def normalise(self):
        
        rtn = self
        
        if not self.normalised:
            idem, squared = self.is_idempotent()
            if idem:
    
                norm_fac = self.get_normalisation_factor(squared = squared)
                if norm_fac is not None:
                    rtn = self.multiply_by_constant(norm_fac)
                    rtn.set_normalised(True)
                else:
                    print('Error: operator cannot be normalised: operator is vanishing.')
            else:
                print('Error: operator cannot be normalised: is not idempotent')
                
        return rtn
    
    
    def get_normalisation_factor(self, squared = None, verbose = False):
        
        norm_fac = None
        if not squared:
            idem, squared = self.is_idempotent()
            if idem:

                norm_fac = fr.Fraction(1,1)



                if len(self.cycle_sum) > 0:

                    #norm_fac = self.get_overall_factor()/squared.get_overall_factor()

                    first_term = self.cycle_sum[0]['cycle']
                    first_pref = self.cycle_sum[0]['prefactor']
                    sq_pref = None
                    ind = 0
                    while sq_pref is None:
                        sq_term = squared.get_sum()[ind]['cycle']

                        if sq_term.equivalent(first_term):
                            sq_pref = squared.get_sum()[ind]['prefactor']

                            norm_fac = fr.Fraction(str(first_pref/sq_pref))
                        ind += 1
                else:
                    norm_fac = None
            elif verbose:
                print('Cannot normalise, operator is not idempotent.')
        else:
            norm_fac = fr.Fraction(1,1)

            if len(self.cycle_sum) > 0:

                #norm_fac = self.get_overall_factor()/squared.get_overall_factor()

                first_term = self.cycle_sum[0]['cycle']
                first_pref = self.cycle_sum[0]['prefactor']
                sq_pref = None
                ind = 0
                while sq_pref is None:
                    sq_term = squared.get_sum()[ind]['cycle']

                    if sq_term.equivalent(first_term):
                        sq_pref = squared.get_sum()[ind]['prefactor']

                        norm_fac = fr.Fraction(str(first_pref/sq_pref))
                    ind += 1
            else:
                norm_fac = None
                if verbose:
                    print('Operator squared has no terms in common with operator.')
                
        return norm_fac
        

    
    
    def is_idempotent(self):
        
        idem = False
        slf = cp.deepcopy(self)
        squared = slf.act_on(slf)
        
        if len(squared.get_sum()) == len(slf.get_sum()):
            if squared.is_equivalent_to(slf, up_to_factor = True):
                idem = True
        
        return idem, squared
    
    
    def is_equivalent_to(self, cycsum, up_to_factor = False):
        
        equiv = False
        cyc_sum = cycsum
        slf = cp.deepcopy(self)

        if (not self.write_as_cycles() in '0') and (not cycsum.write_as_cycles() in '0'):
            pref_slf = 1
            pref_cyc = 1
            if up_to_factor:
                pref_slf = slf.get_identity_factor()
                pref_cyc = cyc_sum.get_identity_factor()
                #slf = slf.multiply_by_constant(pref_cyc)
                #cyc_sum = cyc_sum.multiply_by_constant(pref_slf)
                

            summed = slf.sum_with(cyc_sum, factor = -fr.Fraction(pref_slf,pref_cyc))
            #print(summed.write_as_cycles())
            
            if summed.write_as_cycles() in '0':
                equiv = True
        
        return equiv


    def simplify(self):

        apps = self.cycle_sum
        dic_lis = [cyc.update(cyc['cycle'].cycle) for cyc in apps]

        imgs = [apps[ind]['image'] for ind in range(len(apps))]

        k = sorted(imgs)
        imgs = list(k for k, _ in itertools.groupby(k))
        doms = [sorted(lis) for lis in imgs]

        prefs = [sum([apps[ind]['prefactor'] for ind in range(len(apps)) if apps[ind]['image'] == img]) for img in imgs]

        cycs = [Cycle(domain = doms[ind], image = imgs[ind]) 
                for ind in range(len(imgs)) if prefs[ind].numerator != 0]
        
        prefs = [pref for pref in prefs if pref.numerator != 0]

        new_sum = CycleSum(cycs, prefactors = prefs)

        return new_sum
    
    
    def simplify_depricated(self):
    
        simplified = []
        prefactors = []
        prefactors_str = []
        
        cycles = self.get('cycle')
        prefacs = self.get('prefactor')
        prefacs_str = self.get('prefactor_str')

        for op_ind in range(len(cycles)):

            op = cycles[op_ind]
            pref = prefacs[op_ind]
            pref_str = prefacs_str[op_ind]
            
            appended = False
            for sim_ind in range(len(simplified)):

                sim = simplified[sim_ind]

                if (op.get('domain') == sim.get('domain')) and (op.get('image') == sim.get('image')):

                    if str(float(prefactors[sim_ind]) + float(pref)) in '0.0000000000000000000000000000000000000000':
                    #if str(prefactors[sim_ind]/pref) in '-1.0':
                        prefactors[sim_ind] = fr.Fraction(0,1)
                    else:
                        prefactors[sim_ind] = fr.Fraction(str(prefactors[sim_ind])) + fr.Fraction(str(pref))

                    appended = True
                    
                elif (len(op.get('domain')) == 0) and (len(sim.get('domain')) == 0):
                    if str(float(prefactors[sim_ind]) + float(pref)) in '0.0000000000000000000000000000000000000000':
                    #if str(prefactors[sim_ind]/pref) in '-1.0':
                        prefactors[sim_ind] = fr.Fraction(0,1)
                    else:
                        prefactors[sim_ind] = fr.Fraction(str(prefactors[sim_ind])) + fr.Fraction(str(pref))
                    appended = True
                    
                elif op.as_str() == sim.as_str():
                    if str(prefactors[sim_ind]/pref) in '-1.0':
                        prefactors[sim_ind] = fr.Fraction(0,1)
                    else:
                        prefactors[sim_ind] = fr.Fraction(str(prefactors[sim_ind])) + fr.Fraction(str(pref))

                    appended = True

            if not appended:
                    simplified.append(op)
                    prefactors.append(pref)
                    prefactors_str.append(pref_str)

        return CycleSum(simplified, prefactors = prefactors, str_prefactors = prefactors_str)
    
    
    
    
    def act_on(self, operator):
        
        op1 = self.cycle_sum.copy()
        op2 = operator.get_sum().copy()
    
        final_ops = []
        prefactors = []
        prefacs_str = []

        for op1_ind in range(len(op1)):

            prefac1 = op1[op1_ind]['prefactor']
            prefac1_str = op1[op1_ind]['prefactor_str']

            for op2_ind in range(len(op2)):

                prefac2 = op2[op2_ind]['prefactor']
                prefac2_str = op2[op2_ind]['prefactor_str']
                cycle = op1[op1_ind]['cycle'].act_on(op2[op2_ind]['cycle'])

                #op = Cycle()
                prefactors.append(fr.Fraction(str(prefac1*prefac2)))
                prefacs_str.append('('+prefac1_str+')' + r"\times" + '('+prefac2_str+')')
                #op.set_cycle({'image': cycle.get('image'), 'domain': cycle.get('domain')})

                final_ops.append(cycle)
                
        cycle_list = CycleSum(final_ops, prefactors = prefactors, str_prefactors = prefacs_str)

        return cycle_list.simplify()


    def act_on_depricated(self, operator):
        
        op1 = self.cycle_sum.copy()
        op2 = operator.get_sum().copy()
    
        final_ops = []
        prefactors = []
        prefacs_str = []

        for op1_ind in range(len(op1)):

            prefac1 = op1[op1_ind]['prefactor']
            prefac1_str = op1[op1_ind]['prefactor_str']

            for op2_ind in range(len(op2)):

                prefac2 = op2[op2_ind]['prefactor']
                prefac2_str = op2[op2_ind]['prefactor_str']
                cycle = op1[op1_ind]['cycle'].act_on(op2[op2_ind]['cycle'])

                op = Cycle()
                prefactors.append(fr.Fraction(str(prefac1*prefac2)))
                prefacs_str.append('('+prefac1_str+')' + r"\times" + '('+prefac2_str+')')
                op.set_cycle({'image': cycle.get('image'), 'domain': cycle.get('domain')})

                final_ops.append(op)
                
        cycle_list = CycleSum(final_ops, prefactors = prefactors, str_prefactors = prefacs_str)

        return cycle_list.simplify_depricated()
    
    
    
    
    def write_as_deltas(self, left_label = 'x', right_label = 'y'):
    
        deltas = 0
        for op in self.cycle_sum:

            entries_left = op['cycle'].get('image')
            entries_right = op['cycle'].get('domain')
            prefactor = op['prefactor']
            delta = int(prefactor)

            for ind in range(len(entries_left)):

                delta = delta*sp.KroneckerDelta(sp.symbols(left_label + str(entries_left[ind])),
                                                sp.symbols(right_label + str(entries_right[ind])))
            deltas += delta
        return deltas
    
    

    def sub_cycles_as_str(self, op):
        
        left = op['cycle'].get('image')
        right = op['cycle'].get('domain')

        cycles = '('

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

                    ri_start = 0
                    cycles = cycles + str(cycle).replace('[', '').replace(']',')')

                    cycle = []

                    if len(indices) > 0:
                        ri_start = indices[0]
                        ind = right.index(ri_start)

                        cycles = cycles + '('

                    else: 
                        count = len(right)

                else:
                    ind = right.index(lf)

                count += 1
        if not (cycles[-1] is ')'):
            cycles = cycles + ')'

        return cycles



    def write_as_cycles(self):
        
        ops = self.cycle_list
        pref = self.prefactors

        ops = self.get('cycle')
        pref = self.get('prefactor')

        cycles = ''
        #for op in ops:
        for ind in range(len(ops)):

            prefac_str = ''
            prefac = pref[ind]#op['prefactor']
            #prefac_string = op['prefactor_str']
            if prefac > 0:
                if prefac is 1:
                    prefac_str = '+'
                else:
                    prefac_str = '+' + str(prefac)

            else:
                if prefac is -1:
                    prefac_str = '-'
                else:
                    prefac_str = str(prefac)

            if prefac != 0:
                cycles += prefac_str + ops[ind].as_str()#op['cycle'].as_str()
                
        if len(cycles) == 0:
            cycles = '0'

        return cycles
    
    
    def trace(self, N = 1, num_lines = None):
        
        tr = 0
        
        for cyc in self.cycle_sum:
            
            pref = cyc['prefactor']
            tr_cyc = cyc['cycle'].trace(N = N, num_lines = num_lines)
            
            tr = tr + pref*tr_cyc
            
        return tr