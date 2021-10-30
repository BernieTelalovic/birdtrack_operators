import tableau_utils as tu
import numpy as np
import copy as cp
import fractions as fr


class Diagram():
    
    def __init__(self, size_array, multiplicity = 1):
    
        if 0 in size_array:
            size_array.remove(0)
        diag = [[[]]*ind for ind in size_array]

        self.diagram = {'diagram': diag, 'prefactor': fr.Fraction(multiplicity)}

    
    def complex_conjugate(self, N, multiplicity = 1):
        
        max_col = max(self.size_array())
        sizes = []
        for ind in range(N):

            row_len = 0
            if ind < len(self.diagram['diagram']):
                row_len = len(self.diagram['diagram'][-ind])
            sizes.append(max_col - row_len)
        sizes.reverse()
        cc = Diagram(sizes, multiplicity = multiplicity)

        return cc


    def numerator(self, N):

        num = 1

        for ind_r in range(len(self.diagram['diagram'])):

            for ind_c in range(len(self.diagram['diagram'][ind_r])):

                num = num*(N-ind_r+ind_c)

        return num



    def hook_length(self):

        length = 1

        for ind_r in range(len(self.diagram['diagram'])):

            for ind_c in range(len(self.diagram['diagram'][ind_r])):

                row_entries = self.diagram['diagram'][ind_r][ind_c:]
                column_entries = [[] for ind in self.diagram['diagram'][ind_r:] if len(ind) > ind_c]

                num = len(row_entries)+len(column_entries)-1

                length = length*num

        return length


    def dimension(self, N):

        return self.numerator(N)/self.hook_length()


    def simplify_diagram_list(diags):

        simplified = []

        for diag in diags:

            appended = False
            for sim_ind in range(len(simplified)):

                sim = simplified[sim_ind]
                
                if (diag.size_array() == sim.size_array()):

                    simplified[sim_ind].set_property('prefactor', fr.Fraction(sim.get('prefactor')) + fr.Fraction(diag.get('prefactor')))
                    appended = True

            if not appended:
                    simplified.append(diag)

        return simplified


    def direct_sum(diags):

        diagrams = Diagram.simplify_diagram_list(diags)

        return diagrams


    def size_array(self):

        return [len(row) for row in self.diagram['diagram']]
    


    def is_diagram(self):

        diag_sizes = self.size_array()
        sz = list(sorted(diag_sizes.copy(),reverse=1))

        return max([abs(diag_sizes[ind] - sz[ind]) for ind in range(len(sz))]) == 0


    def get_column(self, ind):

        return [row[ind] for row in self.diagram['diagram'] if ind < len(row)]

    def get_row(self, ind):

        return self.diagram['diagram'][ind]
    
    
    def get(self, prpty):
        
        return self.diagram.get(prpty, None)
    
    
    def set_property(self, prpty, val):
        
        if type(val) is type(self.diagram[prpty]):
        
            self.diagram[prpty] = val
        else:
            typ_prop = str(type(self.diagram[prpty]))
            
            print('Error: value type '+str(type(val)) + ' is not '+ prpty + 'type ' + typ_prop)


    def column_admissable(self):

        diag = self.diagram['diagram']
        admissable = True
        num_colmn = np.max([len(row) for row in diag])

        for col_ind in range(num_colmn):

            column = self.get_column(col_ind)
            entries = [[en[0]] for en in column if len(en)>0]

            for ent in entries:

                if len(ent) > 0:

                    entry = ent[0]
                    eq_entries = []
                    for en in entries:
                        if len(en) > 0:
                            if en[0] == entry:
                                eq_entries.append(entry)

                    if len(eq_entries) > 1:
                        return False

        return admissable


    def row_admissable(self):

        diag = self.diagram['diagram']
        admissable = True

        for row_ind in range(len(diag)):

            row = self.get_row(row_ind)

            for ind in range(len(row)):

                ind_col = len(row)-ind-1
                entry = row[ind_col]

                if len(entry) > 0:
                    read_from_entry = []
                    for r_idx in range(len(diag)):

                        if r_idx < row_ind:

                            prev_row = cp.deepcopy([ent[0] for ent in diag[r_idx] if len(ent)>0])
                            if len(prev_row) > 0:
                                prev_row.reverse()

                                for ent in prev_row:
                                    read_from_entry.append(ent)
                        elif r_idx == row_ind:
                            for c_idx in range(len(row)):

                                col_index = len(row) - c_idx - 1

                                if col_index >= ind_col:

                                    ent_this_row = row[col_index]
                                    if len(ent_this_row) > 0:

                                        read_from_entry.append(ent_this_row[0])

                    if len(read_from_entry) > 0:

                        for num_idx in range(len(read_from_entry)):

                            series = []

                            series = read_from_entry[0:num_idx+1]

                            for num in series:

                                nums = [x for x in series if x == num]

                                for les_num_ind in range(num):

                                    les_num = les_num_ind + 1
                                    less_nums = [y for y in series if y == les_num]

                                    if les_num < num:
                                        if len(nums) > len(less_nums):

                                            return False



        return admissable



    def remove_unadmissable_diags(self,diags):

        admissable = []

        for diag in diags:

            if diag.column_admissable():
                if diag.row_admissable():

                    admissable.append(cp.deepcopy(diag))

        return admissable


    def is_equivalent_to(self,diag2):

        diag1 = self
        di1 = diag1.get('diagram')
        di2 = diag2.get('diagram')

        row_equiv = []
        if len(di1) == len(di2):

            rows1 = diag1.size_array()
            rows2 = diag2.size_array()

            if max(abs(rows1[ind] - rows2[ind]) for ind in range(len(rows1))) == 0:

                for r_ind in range(len(rows1)):
                    equiv = False

                    d1_row = diag1.get_row(r_ind)
                    d1_entries = [ent[0] for ent in d1_row if len(ent) > 0]

                    d2_row = diag2.get_row(r_ind)
                    d2_entries = [tr[0] for tr in d2_row if len(tr) > 0]

                    if len(d1_entries) == len(d2_entries):

                        diffs = [abs(d1_entries[ind] - d2_entries[ind]) for ind in range(len(d1_entries))]
                        diffs.append(0)

                        if max(diffs) == 0:
                            equiv = True

                    row_equiv.append(equiv)
            else:
                row_equiv.append(False)
        else:
            row_equiv.append(False)

        equivalent = True
        if False in row_equiv:
            equivalent =  False

        return equivalent


    def remove_duplicates(self,diags):

        no_duplicates = []

        for diag in diags:
            diag_dupd = False

            for no_dups in no_duplicates:

                if diag.is_equivalent_to(no_dups):
                    diag_dupd = True

            if not diag_dupd:
                no_duplicates.append(diag)

        return no_duplicates





    def remove_indices(self,diags):

        new_diags = []
        for diag in diags:

            size = cp.deepcopy(diag).size_array()
            new_diag = Diagram(size, multiplicity = fr.Fraction(diag.get('prefactor')))

            new_diags.append(new_diag)
        return new_diags

            
    


    def direct_multiple(self, diag2, verbose = False, collect_terms = True):

        diag1 = self
        diag2_sizes = diag2.size_array()
        diag2_enum = [[rw+1]*diag2_sizes[rw] for rw in range(len(diag2_sizes))]

        diag2_enum_cp = diag2_enum.copy()

        diags = []

        new_diagrams = [cp.deepcopy(diag1)]

        for d2_r in diag2_enum:
            for d2_ind in range(len(d2_r)):

                new_diags = []
                entry = d2_r[d2_ind]

                for di_ind in range(len(new_diagrams)):

                    di = cp.deepcopy(new_diagrams[di_ind])
                    di_cp = cp.deepcopy(di)

                    for pl in range(len(di.get('diagram'))):

                        diag = cp.deepcopy(di_cp.get('diagram'))
                        diag[pl].append([entry])

                        candidate = Diagram([1], multiplicity = fr.Fraction(di.get('prefactor')))
                        candidate.set_property('diagram', diag)

                        if candidate.is_diagram():
                            new_diags.append(candidate)

                    diag2 = cp.deepcopy(di.get('diagram'))

                    diag2.append([[entry]]) # append also to a new row at the bottom
                    
                    candidate = Diagram([1], multiplicity = fr.Fraction(di.get('prefactor')))
                    candidate.set_property('diagram', diag2)

                    if candidate.is_diagram():
                        new_diags.append(candidate)

                new_diagrams = self.remove_unadmissable_diags(self.remove_duplicates(new_diags))

        rtn_diags = []
        if not verbose:
            rtn_diags = self.remove_indices(new_diagrams)
        else:
            rtn_diags = new_diagrams
            
        if collect_terms:
            rtn_diags = Diagram.direct_sum(rtn_diags)

        return rtn_diags




    def direct_multiple_list(diags, verbose = False, collect_terms = True):

        diags_cp = cp.deepcopy(diags)
        new_diagrams = []

        if len(diags) > 0:
            new_diagrams = [diags[0]]
            diags_cp = diags_cp[1:]

        while len(diags_cp) > 0:

            new_diags = cp.deepcopy(new_diagrams)
            nds = []

            for nd_ind in range(len(new_diags)):

                nd = cp.deepcopy(new_diags[nd_ind])

                multiples = nd.direct_multiple(cp.deepcopy(diags_cp[0]), verbose = verbose, collect_terms = collect_terms)

                for mult in multiples:
                    nds.append(mult)

            diags_cp.remove(diags_cp[0])
            new_diagrams = cp.deepcopy(nds)

        return new_diagrams
    
    
    def first_occurence(self, N):
        
        sizes = self.size_array()
        max_col = max(sizes)
        arr = []
        for ind in range(N):
            if ind < len(sizes):
                arr_row = [1]*(sizes[ind])
                other = [-1]*(max_col - (sizes[ind]))
                row = arr_row + other
                arr.append(row)
            else:
                arr.append([-1]*max_col)
        first_occ = None
        for col in range(max_col):

            mins = 0
            ons = 0

            for ind in range(col+1):
                column_bef = [row[ind] for row in arr if ind < len(row)]
                mins += len([ent for  ent in column_bef if ent == -1])

            for ind in range(max_col - col):
                ind_a = max_col - ind
                column_aft = [row_a[ind_a] for row_a in arr if ind_a < len(row)]
                ons += len([ent for ent in column_aft if ent == 1])

            if (mins == ons) and mins > 0:
                first_occ = mins

        if max_col ==1 and len(sizes) == N:
            first_occ = 0

        return first_occ
    
    

    def write_diagram(self):

        strin = ''

        di = self.get('diagram')
        weight = self.get('prefactor')

        strin = strin + str(weight)+'\n'+str(di[0])[1:-1].replace(',','')+'\n'

        for row_ind in range(len(di)):
            if row_ind+1 < len(di):
                strin = strin + str(di[row_ind+1])[1:-1].replace(',','')+'\n'
        return strin

    
    def write_diagrams(diags):

        strin = ''
        for diag in diags:

            di = diag.get('diagram')
            weight = diag.get('prefactor')

            strin = strin + diag.write_diagram() #str(weight)+'\n'+str(di[0])[1:-1].replace(',','')+'\n'

           # for row_ind in range(len(di)):
           #     if row_ind+1 < len(di):
           #         strin = strin + str(di[row_ind+1])[1:-1].replace(',','')+'\n'
        return strin