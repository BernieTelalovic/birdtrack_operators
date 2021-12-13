from cycle import Cycle
#from cycle_sum import CycleSum

import numpy as np
import copy as cp
import fractions as fr
from birdtrack_operator import BirdtrackOperator
from cycle_sum import CycleSum
from symmetriser import Symmetriser
from diagram import Diagram


class Tableau(Diagram):
    
    
    def __init__(self, diagram, prefactor = 1, prefactor_str = ''):
    
        pruned_diagram = []
        for row in diagram:
            if len(row) > 0:
                pruned_diagram.append(row)
                
    
        sizes = [len(row) for row in pruned_diagram]
        self.org_diag = pruned_diagram
        diag = Diagram(sizes)
        
        full_diagram = []
        for row in pruned_diagram:
            full_row = []
            for ent in pruned_diagram:
                full_row.append([ent])
            full_diagram.append(full_row)

        diag.set_property('diagram', pruned_diagram)
    
        self.diagram = diag
        self.prefactor = fr.Fraction(prefactor)
        self.prefactor_str = prefactor_str
        if len(prefactor_str) == 0:
            self.prefactor_str = str(prefactor)
            
            
    def get(self, prpty):
        
        if prpty is 'diagram':
            return self.diagram
        elif prpty is 'original_diagram':
            return self. org_diag
        elif prpty is 'prefactor':
            return self.prefactor
        elif prpty is 'prefactor_str':
            return self.prefactor_str
        else:
            return None
        
    def get_row(self, ind):
        return self.diagram.get_row(ind)
    
    def get_column(self, ind):
        return self.diagram.get_column(ind)
        
        
    def set_diagram(self, diag):
        self.diagram = diag
    
    def set_prefactor(self, pref):
        self.prefactor = pref
        
    def set_prefactor_str(self, pref):
        self.prefctor_str = str(pref)
        
        
    def complex_conjugate(self):
        
        cc = self.diagram.complex_conjugate()
        return cc
    
    def dimension(self):
        return self.diagram.dimension()
    
    def hook_length(self):
        return self.diagram.hook_length()
    
    def first_occurence(self):
        return self.diagram.first_occurence()
    
    def size_array(self):
        return self.diagram.size_array()
    
    
    def row_wise_entries(self):
        
        ents = []
        
        diag = self.diagram.get('diagram')
        for row in diag:
            for col in row:
                ents.append(col)
                
        return ents
    
    def column_wise_entries(self):
        
        ents = []
        diag = self.diagram.get('diagram')
        sizes = self.diagram.size_array()
        for col_ind in range(max(sizes)):
            
            column = self.get_column(col_ind)
            for col in column:
                ents.append(col)
                
        return ents
    
    
    def list_of_columns(self):
        
        sizes = self.size_array()
        max_col = max(sizes)
        
        columns = []
        for col_ind in range(max_col):
            columns.append(self.diagram.get_column(col_ind))
            
        return columns
    
    def list_of_rows(self):
        
        sizes = self.size_array()
        
        rows = []
        for row_ind in range(len(sizes)):
            rows.append(self.diagram.get_row(row_ind))
            
        return rows
    
        
    def is_standard(self):
        
        standard = False
        rows = self.list_of_rows()
        columns = self.list_of_columns()
        col_std = True
        row_std = True
        
        for row in rows:
            row_strict = all(i < j for i, j in zip(row, row[1:]))
            
            if not row_strict:
                row_std = False
        
        for column in columns:
            column_strict = all(i < j for i, j in zip(column, column[1:]))
            
            if not column_strict:
                col_std = False
        
        if col_std and row_std:
            standard = True
            
        return standard
        
        
        
    def is_semistandard(self):
        
        standard = False
        rows = self.list_of_rows()
        columns = self.list_of_columns()
        col_std = True
        row_std = True
        
        for row in rows:
            row_strict = all(i <= j for i, j in zip(row, row[1:]))
            
            if not row_strict:
                row_std = False
        
        for column in columns:
            column_strict = all(i < j for i, j in zip(column, column[1:]))
            
            if not column_strict:
                col_std = False
        
        if col_std and row_std:
            standard = True
            
        return standard
    
    
    def parent(self):
        
        diag = cp.deepcopy(self.org_diag)
        maxval = np.max([np.max(row) for row in diag])
        
        parent_diag = [[val for val in row if val != maxval] for row in diag]
        tab = Tableau(parent_diag, prefactor = self.prefactor, prefactor_str = self.prefactor_str)
        
        return tab
    
    def ancestor(self, n):
        
        tab = self
        if n < len(sum(self.org_diag, [])):
            for ind in range(n):
                tab = tab.parent()
        else:
            print('Error: '+str(n)+"-th ancestor does not exist for tableau of "+ str(len(sum(self.org_diag, []))) + 'boxes.')
        return tab
    
    
    def is_ancestor_of(self, tab):
        
        len_self = len(sum(self.org_diag, []))
        len_tab = len(sum(tab.get('original_diagram'), []))
        
        tab_cp = cp.deepcopy(tab)
        ancestor = False
        while len_self <= len_tab:
            
            if tab_cp.get('original_diagram') == tab.get('original_diagram'):
                ancestor = True
            else:
                tab_cp = tab_cp.parent()
                len_tab = len(sum(tab_cp.get('original_diagram'), []))
        return ancestor
        
            
    def remove_indices(tableau):

        sizes = Diagram.size_array(tableau)

        return Diagram(sizes, multiplicity = tableau.get_property('prefactor'))


    def remove_brackets(ls):

        removed = []

        for rw in ls:

            entries = []

            for ent in rw:
                if len(ent) > 0:
                    entries.append(ent[0])
                else:
                    entries.append(0)

            removed.append(entries)

        return removed
    
    
    def write_diagram(self):

        strin = ''

        di = self.diagram.get('diagram')
        weight = self.diagram.get('prefactor')

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

            strin = strin + diag.write_diagram()
        return strin
    
    
    def as_operator(self, order = 'r->l'):
        
        rows = self.list_of_rows()
        columns = self.list_of_columns()
        
        ops = []
        
        for column in columns:
            if len(column) > 1:
                ops.append(Symmetriser(column, antisym = True))
        for row in rows:
            if len(row) > 1:
                ops.append(Symmetriser(row))

        
        prefac = self.prefactor
        prefac_str = self.prefactor_str
        
        return BirdtrackOperator(ops, prefactor = prefac, prefactor_str = prefac_str, order = order)
    
    
    def row_word(self):
        
        word = []
        
        for row in self.org_diag:
            for ent in row:
                word.append(ent)
            
        return word
    
    
    def column_word(self):
        
        word = []
        
        row_len = [len(row) for row in self.org_diag]
        
        for ind in range(np.max(row_len)):
            column = self.get_column(ind)
            
            for ent in column:
                word.append(ent)
                
        return(word)
    
    
    def row_ordered(self):
        
        row_word = self.row_word()
        return all(x<y for x, y in zip(row_word, row_word[1:]))
    
    
    def column_ordered(self):
        
        column_word = self.column_word()
        return all(x<y for x, y in zip(column_word, column_word[1:]))
    
    def row_mold(self):
        
        mold = 0
        tab = cp.deepcopy(self)
        while not tab.row_ordered():
            mold = mold+1
            tab = tab.parent()
            
        return mold
            
            
    def column_mold(self):
    
        mold = 0
        tab = cp.deepcopy(self)
        while not tab.column_ordered():
            mold = mold+1
            tab = tab.parent()
            
        return mold
    
    
    def as_mold_operator(self, normalise = True):
        
        tableau = {'rows': [row for row in self.list_of_rows() if len(row) > 1], 
                   'columns': [col for col in self.list_of_columns() if len(col) > 1]}
        
        row_mold = self.row_mold()
        column_mold = self.column_mold()
        
        mold = []
        
        if row_mold <= column_mold:
            
            if not row_mold%2: # if it's even
                
                mold = [Symmetriser(row) for row in tableau['rows']] +\
                       [Symmetriser(column, antisym = True) for column in tableau['columns']] +\
                       [Symmetriser(row) for row in tableau['rows']]
                
                next_op = {'type': 'columns', 'antisym': True} # antisyms next
                next_op1 = {'type': 'rows', 'antisym': False} # antisyms next
                
                anc = cp.deepcopy(self)

                for ind in range(row_mold):
                    
                    anc = anc.parent()
                    tableau = {'rows': [row for row in anc.list_of_rows() if len(row) > 1], 
                               'columns': [col for col in anc.list_of_columns() if len(col) > 1]}
                    
                    mold = [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]] +\
                           mold +\
                           [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]]
                    
                    next_op2 = next_op # next type of operator
                    next_op = next_op1
                    next_op1 = next_op2
                
            else: # if it's odd
                
                mold = [Symmetriser(column, antisym = True) for column in tableau['columns']] +\
                       [Symmetriser(row) for row in tableau['rows']] +\
                       [Symmetriser(column, antisym = True) for column in tableau['columns']]
                
                next_op = {'type': 'rows', 'antisym': False} # antisyms next
                next_op1 = {'type': 'columns', 'antisym': True} # antisyms next
                
                anc = cp.deepcopy(self)

                for ind in range(row_mold):
                    
                    anc = anc.parent()
                    tableau = {'rows': [row for row in anc.list_of_rows() if len(row) > 1], 
                               'columns': [col for col in anc.list_of_columns() if len(col) > 1]}
                    
                    mold = [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]] +\
                           mold +\
                           [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]]
                    
                    next_op2 = next_op # next type of operator
                    next_op = next_op1
                    next_op1 = next_op2
        
        
        else:
            
            if not column_mold%2: # if it's even
                
                mold = [Symmetriser(column, antisym = True) for column in tableau['columns']] +\
                       [Symmetriser(row) for row in tableau['rows']] +\
                       [Symmetriser(column, antisym = True) for column in tableau['columns']]
                
                next_op = {'type': 'rows', 'antisym': False} # antisyms next
                next_op1 = {'type': 'columns', 'antisym': True} # antisyms next
                
                anc = cp.deepcopy(self)

                for ind in range(row_mold):
                    
                    anc = anc.parent()
                    tableau = {'rows': [row for row in anc.list_of_rows() if len(row) > 1], 
                               'columns': [col for col in anc.list_of_columns() if len(col) > 1]}
                    
                    mold = [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]] +\
                           mold +\
                           [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]]
                    
                    next_op2 = next_op # next type of operator
                    next_op = next_op1
                    next_op1 = next_op2
                
            else: # if it's odd
                
                mold = [Symmetriser(row) for row in tableau['rows']] +\
                       [Symmetriser(column, antisym = True) for column in tableau['columns']] +\
                       [Symmetriser(row) for row in tableau['rows']]
                
                next_op = {'type': 'columns', 'antisym': True} # antisyms next
                next_op1 = {'type': 'rows', 'antisym': False} # antisyms next
                
                anc = cp.deepcopy(self)

                for ind in range(column_mold):
                    
                    anc = anc.parent()
                    tableau = {'rows': [row for row in anc.list_of_rows() if len(row) > 1], 
                               'columns': [col for col in anc.list_of_columns() if len(col) > 1]}
                    
                    mold = [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]] +\
                           mold +\
                           [Symmetriser(dom, antisym = next_op['antisym']) for dom in tableau[next_op['type']]]
                    
                    next_op2 = next_op # next type of operator
                    next_op = next_op1
                    next_op1 = next_op2
        
        prefac = self.prefactor
        prefac_str = self.prefactor_str
        
        mold_op = BirdtrackOperator(mold, prefactor = prefac, prefactor_str = prefac_str, normalise = normalise)
            
        return mold_op