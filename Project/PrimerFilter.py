from Bio.SeqUtils import MeltingTemp as mt
from statistics import mean
from itertools import product


class PrimerFilter():
    def __init__(self, primers, config):
        super(PrimerFilter, self).__init__()
        self.primers = primers
        self.config = config

    
    def get_tm_avg(self, seq):
        tm_methods = ['Tm Wallace', 'Tm NN', ]
        tm_methods.extend([f'Tm NN DNA_NN{i} Table' for i in range(1, 5)])
        tm_methods.extend([f'Tm GC valueset {i}' for i in range(1, 9)])

        results = []
        for method in tm_methods:
            if method == 'Tm Wallace':
                tm = round(mt.Tm_Wallace(seq, strict=False))
            elif method == 'Tm NN':
                tm = round(mt.Tm_NN(seq, strict=False), 2)
            elif 'DNA_NN' in method:
                nn_table = getattr(mt, method.split(' ')[2])
                tm = round(mt.Tm_NN(seq, nn_table=nn_table, strict=False), 2)
            elif 'GC' in method:
                valueset = int(method[-1])
                tm = round(mt.Tm_GC(seq, valueset=valueset, strict=False), 2)

            results.append(tm)

        tm = mean(results)

        return tm


    def calc_primers_prop(self, seq):
        seq = seq.upper()
        size = len(seq)
        tm = self.get_tm_avg(seq)
        return size, tm


    def get_primers_prop(self):
        self.primers['Size'], self.primers['Tm'] = zip(*self.primers['nuc'].apply(self.calc_primers_prop))

        return self.primers

    def get_combinations(self, seq):
        ambiguity_table = {
            'A': ['A'],
            'C': ['C'],
            'G': ['G'],
            'T': ['T'],
            'W': ['A', 'T'],
            'S': ['C', 'G'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T'],
        }

        combinations = 0

        for i in product(*[ambiguity_table[j] for j in seq]):
            combinations += 1

        return combinations


    def get_valid_amb_comb(self):
        accepted_primers = []
        for primer_id, primer_seq in zip(self.filtered_primers['id'], self.filtered_primers['nuc']):
            combinations = self.get_combinations(primer_seq)
            if combinations <= self.config['CRITERIA']['MAX_AMB_COMB']:
                accepted_primers.append(primer_id)
            else:
                print(primer_seq)

        return accepted_primers


    def filter_primers(self):
        self.primers = self.get_primers_prop()

        min_len = self.config['CRITERIA']['LENGTH']['MIN']
        max_len = self.config['CRITERIA']['LENGTH']['MAX']
        tm_min = self.config['CRITERIA']['TM']['MIN']
        tm_max = self.config['CRITERIA']['TM']['MAX']
        
        self.filtered_primers = self.primers[
            (self.primers['Size'].between(min_len, max_len)) & (self.primers['Tm'].between(tm_min, tm_max))
        ]

        if self.config['CRITERIA']['INOSINE'] is False:
            self.filtered_primers = self.filtered_primers.loc[~self.filtered_primers['nuc'].str.contains('I')]

        if self.config['CRITERIA']['MAX_AMB_COMB'] is not None:
            accepted_primers = self.get_valid_amb_comb()
            self.filtered_primers = self.filtered_primers.loc[self.filtered_primers['id'].isin(accepted_primers)]

        return self.filtered_primers
