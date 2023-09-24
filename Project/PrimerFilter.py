from Bio.SeqUtils import MeltingTemp as mt
from statistics import mean
from itertools import product


class PrimerFilter():
    def __init__(self, primers, config):
        super(PrimerFilter, self).__init__()
        self.primers = primers
        self.config = config


    def extend_ambiguous_dna(self, seq):
        degenerated_table = {
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

        r = []

        # Aplica produto cartesiano nos conjuntos (conjunto = possíveis bases para cada letra)
        for i in product(*[degenerated_table[j] for j in seq]):
            r.append("".join(i))

        return r


    def get_tm_avg(self, seq):
        combinations = self.extend_ambiguous_dna(seq)
        tms = []

        for i in combinations:

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
            tms.append(tm)

        return mean(tms)


    def calc_primers_prop(self, seq):
        seq = seq.upper()
        size = len(seq)
        tm = self.get_tm_avg(seq)
        return size, tm


    def get_primers_prop(self, primers):
        primers['Size'], primers['Tm'] = zip(*primers['nuc'].apply(self.calc_primers_prop))

        return primers

    def get_number_of_combinations(self, seq):
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
            combinations = self.get_number_of_combinations(primer_seq)
            if combinations <= self.config['CRITERIA']['MAX_AMB_COMB']:
                accepted_primers.append(primer_id)

        return accepted_primers


    def filter_primers(self):

        min_len = self.config['CRITERIA']['LENGTH']['MIN']
        max_len = self.config['CRITERIA']['LENGTH']['MAX']
        tm_min = self.config['CRITERIA']['TM']['MIN']
        tm_max = self.config['CRITERIA']['TM']['MAX']
        
        if self.config['CRITERIA']['INOSINE'] is False:
            self.filtered_primers = self.primers.loc[~self.primers['nuc'].str.contains('I')]
            
        if self.config['CRITERIA']['MAX_AMB_COMB'] is not None:
            accepted_primers = self.get_valid_amb_comb()
            self.filtered_primers = self.filtered_primers.loc[self.filtered_primers['id'].isin(accepted_primers)]

        self.filtered_primers = self.get_primers_prop(self.filtered_primers)
        size = self.filtered_primers[self.filtered_primers['Size'].between(min_len, max_len)]
        tm = size[size['Tm'].between(tm_min, tm_max)]


        self.filtered_primers = self.filtered_primers[
            (self.filtered_primers['Size'].between(min_len, max_len)) & (self.filtered_primers['Tm'].between(tm_min, tm_max))
        ]


        return self.filtered_primers
