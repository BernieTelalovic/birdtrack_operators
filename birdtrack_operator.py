from cycle import Cycle
from cycle_sum import CycleSum
from symmetriser import Symmetriser
from diagram import Diagram
import numpy as np
import copy as cp
import fractions as fr
import tableau_utils as tu
import math

class BirdtrackOperator():
    
    def __init__(self, ops, prefactor = 1, prefactor_str = '', order = 'r->l', normalise = True,):
        
        self.domain = list(set(sum([op.get_total_domain() for op in ops if len(ops) > 0], [])))
        
        self.operator_seq = ops
        
        self.order = order
        self.prefactor = fr.Fraction(prefactor)
        self.prefactor_str = prefactor_str
        self.normalised = False
        self.idempotent = None
        self.collapsed = None
        if normalise:
            self.idempotent = self.is_idempotent()
            if self.idempotent:
                if normalise:
                    self.normalise()
                else:
                    self.normalised = False


    def indices_seperate_by_domain(self):

        op_lis = self.operator_seq
        lis_domains = [op.get('domain') for op in op_lis]
        
        seperated_ops_inds = []
        total_domain = []
        app_inds = []
        for ind in range(len(op_lis)):
            appended = True

            if len(list(set(total_domain) & set(lis_domains[ind]))) == 0:
                app_inds.append(ind)
                total_domain = list(set(total_domain+lis_domains[ind]))
            else:
                total_domain = lis_domains[ind]
                seperated_ops_inds.append(app_inds)
                app_inds = [ind]

            if ind == len(op_lis)-1:
                total_domain = []
                seperated_ops_inds.append(app_inds)
                app_inds = []
                
        return seperated_ops_inds


    def transitions_to(self, op):

        op_seperate_indices = op.indices_seperate_by_domain()
        self_seperate_indices = self.indices_seperate_by_domain()

        immediately_zero = False
        self_right_most_ops = [[self.operator_seq[ind].domain, self.operator_seq[ind].antisym] 
                                for ind in self_seperate_indices[-1]]
        op_left_most_ops = [[op.operator_seq[ind].domain, op.operator_seq[ind].antisym] 
                                for ind in op_seperate_indices[0]]
        #print(self_right_most_ops)
        #print(op_left_most_ops)

        for slfop in self_right_most_ops:

            slf_dom = slfop[0]

            for opop in op_left_most_ops:
                    
                if not (slfop[1] is opop[1]):
                        
                    op_dom = opop[0]
                    
                    if len(list(set(slf_dom) & set(op_dom))) > 1:
                        immediately_zero = True
                        #print('immediately zero')


        op_mid_ind = int(len(op_seperate_indices)/2)

        self_mid_ind = int(len(self_seperate_indices)/2)
        self_start_ind = 0
        if len(self_seperate_indices)/2 > 2.5:
            self_start_ind = min(self_seperate_indices[self_mid_ind+2])

        self_side_ops_types = [self.operator_seq[ind].antisym for ind in self_seperate_indices[min([self_mid_ind+1, len(self_seperate_indices)-1])]]
        self_side_typ = max(set(self_side_ops_types), key = self_side_ops_types.count)
        #print(self_side_typ)
        
        op_side_ops_types = [op.operator_seq[ind].antisym for ind in op_seperate_indices[op_mid_ind-1]]
        #print('op off-doms', [op.operator_seq[ind].domain for ind in op_seperate_indices[op_mid_ind-1]])
        op_side_typ = max(set(op_side_ops_types), key = op_side_ops_types.count)
        #print(op_side_typ)
        
        min_ind = 0
        if not (op_side_typ is self_side_typ):
            min_ind = 1

        op_start_ind = 0
        if len(op_seperate_indices)/2 > 2.5:
            #print('here!')
            op_start_ind = max(min(op_seperate_indices[op_mid_ind-min_ind]),1)
            #print(op_start_ind)

        self_right_ops = self
        if len(self_seperate_indices)/2 > 2.5:
            self_ops = self.operator_seq[self_start_ind:]
            #print([op.domain for op in self_ops])
            self_right_ops = BirdtrackOperator(self_ops, normalise = False)

        op_left_ops = op
        if len(op_seperate_indices)/2 > 2.5:
            op_ops = op.operator_seq[:op_start_ind]
            #print([op.domain for op in op_ops])
            op_left_ops = BirdtrackOperator(op_ops, normalise = False)

        collapsed_mid = self_right_ops.deepact_on(op_left_ops)

        cyc_lists = collapsed_mid.get('cycle')
        cyc_prefs = collapsed_mid.get('prefactor')
        
        leftover_ops = self.operator_seq[:self_mid_ind+1]# + op.operator_seq[op_mid_ind-min_ind+1:]
        
        self_op = [[self.operator_seq[ind].domain, self.operator_seq[ind].antisym] 
                   for ind in self_seperate_indices[min([self_mid_ind+1, len(self_seperate_indices)-1])]]
        self_op_total_dom = []
        for p in self_op:
            self_op_total_dom += p[0]
        self_op_total_dom = list(set(self_op_total_dom))

        self_mid_op = [[self.operator_seq[ind].domain, self.operator_seq[ind].antisym] 
                   for ind in self_seperate_indices[self_mid_ind]]

        self_mid_op = [ent for ent in self_mid_op if len(list(set(ent[0]) & set(self_op_total_dom))) <= 1]

        op_op = [[op.operator_seq[ind].domain, op.operator_seq[ind].antisym]
                 for ind in op_seperate_indices[op_mid_ind-min_ind]]
        op_op_total_dom = []
        for p in op_op:
            op_op_total_dom += p[0]
        op_op_total_dom = list(set(op_op_total_dom))

        op_side_op = [[op.operator_seq[ind].domain, op.operator_seq[ind].antisym]
                 for ind in op_seperate_indices[op_mid_ind-min_ind+1]]

        op_side_op = [ent for ent in op_side_op if len(list(set(ent[0]) & set(op_op_total_dom))) <= 1]


        #self_left_most_op = [[self.operator_seq[ind].domain, self.operator_seq[ind].antisym] 
        #           for ind in self_seperate_indices[self_mid_ind]]
        #print('left left op', self_mid_op)
        #print('left op', self_op)
        #op_right_most_op = [[op.operator_seq[ind].domain, op.operator_seq[ind].antisym]
        #         for ind in op_seperate_indices[op_mid_ind-min_ind+2]]
        #print('right op', op_op)
        #print('right op', op_side_op)

        total_prefs = []

        if not immediately_zero:

            for ind in range(len(cyc_lists)):
                
                pref = cyc_prefs[ind]
                
                
                for slfop in self_op+self_mid_op:

                    slf_dom = slfop[0]
                    cyc = cyc_lists[ind].reverse()
                    expanded_cyc_dom = sorted(list(set(cyc.get('domain')+slf_dom)))
                    cyc.expand_domain(expanded_cyc_dom)
                            
                    slf_dom_cycled = []

                    for opop in op_op:
                        
                        if not (slfop[1] is opop[1]):
                            
                            op_dom = opop[0]
                            
                            for p in slf_dom:
                                ind = cyc.cycle['domain'].index(p)
                                slf_dom_cycled.append(cyc.cycle['image'][ind])
                            
                            if len(list(set(slf_dom_cycled) & set(op_dom))) > 1:
                                pref = None

                            
                if not pref is None:
                    total_prefs.append(pref)
                    
        return total_prefs



        
        
    def is_idempotent(self):
        
        idem = False
        if self.collapsed:
            if self.collapsed.is_idempotent():
                idem = True
        else:
            self_cp = cp.deepcopy(self)
            
            collapsed = self_cp.collapse()
            self.collapsed = collapsed
            if collapsed.is_idempotent():
                idem = True
            
        return idem
    
    def is_normalised(self):
        return self.normalised
        
    
    def set_normalised(self, val):
        
        self.normalised = val
        
        
    def set_prefactor(self, val):
        
        self.prefactor = val
        

    def get(self, prpty):
        
        if prpty is 'prefactor':
            return self.prefactor
        elif prpty is 'prefactor_str':
            return self.prefactor_str
        elif prpty is 'symmetrisers':
            return Operator([sym for sym in self.operator_seq if not sym.is_antisym()])
        elif prpty is 'antisymmetrisers':
            return Operator([asym for asym in self.operator_seq if sym.is_antisym()])
        elif prpty is 'domain':
            return self.domain
        elif prpty is 'operator':
            return self.operator_seq
        elif prpty is 'order':
            return self.order


    def symmetriser_lengths(self):
        lis = [math.factorial(len(op.domain)) for op in self.operator_seq]
        return np.product(lis)
        
    
        
    def normalise(self, verbose = False):
        
        norm_fac = None
        normalised = cp.deepcopy(self)
        sym_fac = 1

        if self.operator_seq[0].operator.normalised is False:
            sym_fac = self.symmetriser_lengths()
        
        if self.is_idempotent():
            if (not normalised.is_normalised()):
            
                self_cp = cp.deepcopy(self)

                norm_fac = self_cp.collapse().get_normalisation_factor()
                if not norm_fac is None:
                    self.set_normalised(True)
                    self.set_prefactor(norm_fac*sym_fac)
                    if verbose:
                        print('Operator has been normalised.')
                elif verbose:
                    print('Operator cannot be normalised.')
        
    

    def deepact_on(self, op):
        
        operator = op
        op_order = 'r->l'
        if 'Operator' in str(type(op)):
            op_order = op.get('order')
            operator = op.collapse()

        if operator is not None:
            
            my_order = self.order
            
            if not my_order is op_order:
                operator = operator.reverse()
                
            my_collapsed = self.collapse()

            deepact = my_collapsed.act_on(operator)
            
        return deepact

            
    def collapse(self):
        
        op_seq = self.operator_seq
        prefactor = self.prefactor
        
        collapsed = cp.deepcopy(op_seq[0]).multiply_by_constant(prefactor)

        for ind in range(len(op_seq)):
            
            if ind > 0:
                collapsed = collapsed.act_on(op_seq[ind])
                
        return collapsed
            

    
    def as_tableau(self):
        
        syms = self.operator['symmetrisers']
        asyms = self.operator['antisymmetrisers']
        
        prefac = self.operator['prefactor']
        prefac_str = self.operator['prefactor_str']
        
        
        sym_domains = self.sym_multiplets
        sdom_lengths = [len(dom) for dom in sym_domains]
        asym_domains = self.asym_multiplets
        adom_lengths = [len(dom) for dom in asym_domains]
        
        sorted_sdoms = [x for _,x in sorted(zip(sdom_lengths,sym_domains))]
        sorted_adoms = [x for _,x in sorted(zip(adom_lengths,asym_domains))]
        
        tab_sym = sorted_sdoms
        tab_asym = []
        
        for ind in range(len(sorted_adoms[0])):
            tab_asym[ind] = []
            for adom in sorted_adoms:
                if ind < len(adom):
                    tab_asym[ind].append(adom[ind])
                    
        print(tab_sym)
        print(sorted_adoms)
        
        
    def flip_vertically(self):
        
        op_seq = self.operator_seq
        
        new_seq = []
        
        for ind in range(len(op_seq)):
            
            new_seq.append(op_seq[-ind-1])
            
        return Operator(new_seq, prefactor = self.prefactor, prefactor_str = self.prefactor_str, 
                        order = self.order, normalised = self.normalised)
        
        
    def flip_horizontally(self):
        
        op_seq = self.operator_seq
        
        new_seq = []
        for ind in range(len(op_seq)):
            
            oper = op_seq[ind]
            oper_dom = oper.get_total_domain()
            oper_typ = str(type(oper))
            new_op = None
            if 'ymmetriser' in oper_typ:
                is_antisym = oper.is_antisym()
                
                new_dom = [max(oper_dom) -ent +1 for ent in oper_dom]
                new_op = Symmetriser(new_dom, antisym = is_antisym)
            else:
                print('ERROR: non-symmetriser types not supported yet in flip_horizontally()!!!')
                
            if not new_op is None:
                new_seq.append(new_op)
                
        return Operator(new_seq, prefactor = self.prefactor, prefactor_str = self.prefactor_str, 
                        order = self.order, normalised = self.normalised)
        
        
    def draw_symmetriser(self, lines, draw_cfg):
        
        
        return None
    
    
    
    def categorise_operators(self):
        
        op_seq = self.operator_seq
        op_types = [tu.get_type(str(type(op))) for op in op_seq]
        print(op_types)
        
        categories = []
        
        if len(op_seq) > 0:
            category = [op_seq[0]]
            typ = op_types[0]
            for op_ind in range(len(op_seq)):

                if op_ind > 0:
                    op = op_seq[op_ind]
                    op_typ = op_types[op_ind]
                    
                    if typ is op_typ:
                        category.append(op)
                    else:
                        print(category)
                        categories.append(category)
                        category = [op]
                        typ = op_typ
                        
        return categories 
            
            
            
    def draw_invertical_operators(self, draw_cfg, strin, lines):
    
        line_width = draw_cfg['line width']
        vert_sep = draw_cfg['vertical separation']
        horiz_sep = draw_cfg['horizontal separation']
        cycle_len = draw_cfg['cycle length']
        sym_width = draw_cfg['symmetriser width']
        perm_dist = draw_cfg['permutation distance']
        
        
    def draw_cycle(self, draw_cfg, strin, lines, cycle):
    
        return None
        
        
        
    def write_latex1(self, scale = 1, write_prefac = False):
        
        draw_cfg = tu.read_yaml('latex_config.yaml')

        strin = draw_cfg['start string']
        line_width = draw_cfg['line width']
        vert_sep = draw_cfg['vertical separation']
        horiz_sep = draw_cfg['horizontal separation']
        cycle_len = draw_cfg['cycle length']
        sym_width = draw_cfg['symmetriser width']
        perm_dist = draw_cfg['permutation distance']
        
        domain = np.arange(min(self.domain), max(self.domain)+1)
        lines = {str(num): [(num)*vert_sep*scale, 0, (len(domain) - num + 1)*vert_sep*scale] \
                 for num in domain}

        op_ind = 0
        print(len(self.operator_seq))
        while op_ind < len(self.operator_seq):
            print('fir', op_ind)
            op = self.operator_seq[op_ind]
            typ = str(type(op))
            cycle = None
            dom = op.get("domain")
            print('domain', dom)
            perm_str = ''
            extra_draw = ''
            
            next_op = None
            
            if "ymmetriser" in typ:
                asym = asym = op.is_antisym()
                print('----' + str(asym) + '---- at ind: '+str(op_ind))
                next_asym = asym
                op_ind_x = op_ind
                op = self.operator_seq[op_ind_x]
                dom = op.get("domain")
                prev_top = min(dom)
                prev_bottom = min(dom) - 1
                
                exit = False
                cycle = {'domain': domain, 'image': list(np.zeros(len(domain)))}
                
                while (op_ind_x < len(self.operator_seq)) and (not exit):

                    op = self.operator_seq[op_ind_x]
                    
                    if ("ymmetriser" in str(type(self.operator_seq[op_ind_x]))):
                        
                        next_asym = op.is_antisym()
                        dom = op.get("domain")
                        image = dom
                        print('next op type: '+ str(next_asym) + 'at ind: '+ str(op_ind_x))

                        if next_asym == asym:

                            vert = vert_sep*0.4
                            horiz = abs((horiz_sep - sym_width)/2)
                            asym = op.is_antisym()
                            sym_str = 'symmetrizer'

                            start = np.round(lines[list(lines.keys())[-1]][1] - horiz, 3)
                            end = np.round(lines[list(lines.keys())[-1]][1] - horiz - sym_width, 3)
                            
                            top_line = max([min(dom), prev_bottom+1])
                            prev_top = top_line
                            top = np.round(lines[str(top_line)][2] + vert, 4)
                            
                            bottom_line = top_line + len(dom) - 1
                            bottom = np.round(lines[str(bottom_line)][2] - vert,4)
                            prev_bottom = top_line + len(dom) - 1

                            if asym:
                                sym_str = 'anti'+sym_str
                                
                            extra_draw += "\draw["+ sym_str + "](" + str(end) + "," + str(top) + ")rectangle(" +\
                                         str(start) + "," + str(bottom) +");" + "\n"
                            
                                                        
                            inner_lines = list(np.arange(top_line, bottom_line+1))
                            for ind in range(len(cycle['domain'])):
                                
                                entry = cycle['domain'][ind]
                                if (cycle['image'][ind]-0.5) < 0:
                                    if entry in dom:
                                    
                                        for dom_ent_ind in range(len(dom)):
                                        
                                            dom_ent = dom[dom_ent_ind]
                                            if entry == dom_ent:
                                                cycle['image'][ind] = inner_lines[dom_ent_ind]
                            
                        else:
                            exit = True
                            op_ind = op_ind_x
                    
                    else:
                        exit = True
                        op_ind = op_ind_x
                        
                    print('op_ind',op_ind)
                    
                    op_ind_x += 1
                    if (op_ind_x >= len(self.operator_seq)) and (not exit):
                        op_ind = op_ind_x
                        
                    print('new x', op_ind_x)

                missing_lines = [cycle['domain'][ind] for ind in range(len(cycle['domain'])) \
                                 if not cycle['domain'][ind] in cycle['image']]
                print('miss lin', missing_lines)
                for lin in missing_lines:
                    ex = False
                    for ind in range(len(cycle['image'])):
                        
                        
                        if cycle['image'][ind] - 0.5 < 0 and (not ex):
                            cycle['image'][ind] = lin
                            ex = True
                print('final cyc', cycle)
            
            elif "Cycle" in str(type(op)):
                image = op.get("image")
                
                
            print(extra_draw)  
            print('cycle', cycle)    
            doma = domain
            imag = []
            for dm in doma:
                
                if dm in dom:
                    for ind in range(len(dom)):
                        if dm == dom[ind]:
                            imag.append(image[ind])
                else:
                    imag.append(dm)
            
            for ind in range(len(doma)):
                line = lines[str(doma[ind])][0]
                start = np.round(lines[str(doma[ind])][1], 3)
                target_line = lines[str(imag[ind])][0]
                end = np.round(start - horiz_sep*scale, 3)
                
                to_str = ''
                if line != target_line:
                    to_str = "[wavy]"
                    
                perm_str += "\draw[line width =" + str(line_width) + "pt](" + str(end) +\
                            "," + str(target_line) + ")to"+ to_str +"("+ str(start) + "," + str(line) +");" +"\n"
                
                lines[str(doma[ind])][1] = end
                
            perm_str = perm_str + extra_draw
            
            strin += perm_str
            
        return strin + "\end{tikzpicture}"
    
    
    def type_of_op(op):
    
        if op.is_antisym():
            return 'antisymmetriser'
        else:
            return 'symmetriser'
    
    def separate_operators(self):

        operators = self.operator_seq
        lis_sep = [[operators[0]]]
        lis_sep_ind = 0

        last_ent = type_ent(operators[0])

        for ind in range(len(operators)):
            if ind > 0:
                entry = operators[ind]

                if type_ent(entry) is last_ent:

                    lis_sep[lis_sep_ind].append(entry)
                else:
                    lis_sep.append([entry])
                    lis_sep_ind = lis_sep_ind + 1

                last_ent = type_ent(entry)
                
        return lis_sep
    
    def sort_operators(ops, domain):
        
        domains = [op.domain() for op in ops]
        scores = [sum(dom) for dom in domains]
        
        sorted_ops = [x for _, x in sorted(zip(scores, ops))]
        cyc_images = [x for _, x in sorted(zip(scores, domains))]
        
        cyc_img = []
        for dom in cyc_images:
            for ent in dom:
                cyc_img.append(ent)
                
        for ind in domain:
            fuck = ''
        
        return cyc_img, sorted_ops
    
    
    def draw_operators(ops, domain, lines, config):
        
        strin = ''
        
        
    
    def draw_cycle(cyc_dom, cyc_img, domain, lines, config):
    
        line_width = draw_cfg['line width']
        vert_sep = draw_cfg['vertical separation']
        horiz_sep = draw_cfg['horizontal separation']
        cycle_len = draw_cfg['cycle length']
        perm_dist = draw_cfg['permutation distance']
        strin = ''
        for ind in range(len(cyc_dom)):
            
            line_dom = cyc_dom[ind]
            start_x = lines[str(line_dom)][1]
            start_y = lines[str(line_dom)][2]
            
            line_img = cyc_img[ind]
            end_x = lines[str(line_img)][1] + cycle_len
            lines[str(line_img)][1] += cycle_len
            end_y = lines[str(line_img)][2]
            
            wavy = ''
            
            if line_dom != line_img:
                wavy = '[wavy]'
                
            strin = strin + "\draw[line width =" + line_width +"](" + start_x + "," + start_y +")to" +\
                    wavy + "(" + end_x + "," + end_y + ");"
            
        return strin
    
    def write_latex(self, scale = 1, write_prefac = False):
        
        draw_cfg = tu.read_yaml('latex_config.yaml')

        strin = draw_cfg['start string']
        line_width = draw_cfg['line width']
        vert_sep = draw_cfg['vertical separation']
        horiz_sep = draw_cfg['horizontal separation']
        cycle_len = draw_cfg['cycle length']
        sym_width = draw_cfg['symmetriser width']
        perm_dist = draw_cfg['permutation distance']
        
        domain = np.arange(min(self.domain), max(self.domain)+1)
        lines = {str(num): [(num)*vert_sep*scale, 0, (len(domain) - num + 1)*vert_sep*scale] \
                 for num in domain}
        
        operators_list = self.separate_operators()
        
        cyc_dom = domain
        
        for operators in operators_list:
            
            cyc_img, ops_sorted = sort_operators(operators)
            
            cyc_strin = draw_cycle(cyc_dom, cyc_img, domain, lines, draw_cfg)
            op_strin = draw_operators(ops_sorted, domain, lines, draw_cfg)
            
            cyc_dom = cyc_img
            strin = strin + cyc_strin + op_strin
            
        return strin
            
                         