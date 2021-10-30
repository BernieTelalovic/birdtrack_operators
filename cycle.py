#import sympy as sp
#import gravipy as gp
import numpy as np
import itertools

class Cycle():
    
    def __init__(self, cyc = ''):
        
        cyc = cyc.replace(',', ' ').replace(') (', ')(')
        cycles_dict = []
        final_cycle = None
        if len(cyc) > 0:

            cycle_cut = cyc[1:-1]
            cycles_str = [[int(c) for c in cy.split(' ')] for cy in cycle_cut.split(')(')]

            for cy in cycles_str:

                cycle_dict = {'domain': [], 'image': []}

                for ind in range(len(cy)):
                    cycle_dict['domain'].append(cy[ind])
                    next_ind = ind + 1

                    if next_ind >= len(cy):
                        cycle_dict['image'].append(cy[0])

                    else:
                        cycle_dict['image'].append(cy[next_ind])
                
                cycles_dict.append(cycle_dict)
                
            final_cycle = cycles_dict[0]
            if len(cycles_dict) > 1:
 
                final_cyc = Cycle()
                final_cyc.set_cycle(final_cycle)
                
                for cyc_ind in range(len(cycles_dict)):
                    
                    if cyc_ind >= 1:
                        
                        next_cyc = cycles_dict[cyc_ind]
                        next_cycle = Cycle()
                        next_cycle.set_cycle(next_cyc)
                        
                        final_cyc = final_cyc.act_on(next_cycle)
                        
                final_cycle =  final_cyc.get_cycle()
                
        else:
            final_cycle = {'domain': [], 'image': []}
            
        self.cycle = final_cycle
        
        
    def get(self, prpty):
        
        return self.cycle.get(prpty, None)
        
        
    def get_cycle(self):
        
        return self.cycle
    
    
    def set_cycle(self, cycle):
        
        self.cycle = cycle
        
        
    def trace(self, N = 1, num_lines = None):
        
        domain = self.cycle['domain']
        image = self.cycle['image']
        
        tr = 0
        nums_in_cycle = []
        
        if not num_lines:
            if len(self.cycle['domain']) > 0:
                num_lines = max(self.cycle['domain'])
            else:
                num_lines = 1

        extra_lines = max([0, abs(num_lines - len(domain))])
        
        for ind in range(len(domain)):
            
            dom = domain[ind]
            img = image[ind]
            
            if dom in nums_in_cycle:
                tr = tr+1
                nums_in_cycle = []
            else:
                nums_in_cycle.append(img)
                
        return N**(tr+extra_lines)
        
        
    def reverse(self):
        
        img = self.cycle['domain'].copy()
        dom = self.cycle['image'].copy()
        
        image = [x for _,x in sorted(zip(dom,img))]
        domain = [x for x,_ in sorted(zip(dom,img))]
        
        self.set_cycle({'domain': domain, 'image': image})
        
        return self
        
        
    def act_on(self, cyc):

        #a=self o b= cyc
        a = self.cycle.copy()
        b = cyc.cycle.copy()
        
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
        
        final_cyc = Cycle()
        final_cyc.set_cycle(cycle)

        return final_cyc.remove_single_cycles()
    
    
    def equivalent(self, cyc):
        
        my_img = self.cycle['image']
        cyc_img = cyc.get('image')
        my_dom = self.cycle['domain']
        cyc_dom = cyc.cycle['domain']
        equivalent = False
        
        if (my_img == cyc_img) and (my_dom == cyc_dom):
        #(my_img is [img for img in my_img if img in cyc_img]) and (cyc_img is [img for img in cyc_img if img in my_img]):
            equivalent = True
        elif len(my_img) + len(cyc_img) + len(my_dom) + len(cyc_dom) < 1:
            equivalent = True
            
        return equivalent
    
    
    
    def remove_single_cycles(self):
        
        active_cycle = self
        if len(self.cycle['domain']) > 0:
            ops = self.cycle.copy()


            left = ops['image']
            right = ops['domain']

            active_left = [left[ind] for ind in range(len(right)) if left[ind]!=right[ind]]
            active_right = [right[ind] for ind in range(len(right)) if right[ind]!=left[ind]]

            active_cycle.set_cycle({'image': active_left, 'domain': active_right})

        return active_cycle
    
    
    
    def cycle_num(self):

        perm = self.cycle
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
    
    
    def permutation_order(self):
    
        perm_num = self.cycle_num()

        if perm_num > 0:
            perm_num += 0

        return perm_num
    
    
    
    def as_str(self):
        
        op = self.cycle
        left = op.get('image',[])
        right = op.get('domain', [])

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
                    cycles = cycles + str(np.roll(cycle,1)).replace('[', '').replace(']',')')

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
