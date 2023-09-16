import pandas as pd
import primer3
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from statistics import mean


df = pd.read_csv('primers.csv', sep=';')

def get_tm_avg(seq):
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

# Função para calcular as propriedades dos primers
def calcular_propriedades_primers(sequencia):
    sequencia = sequencia.upper()  # Converter para letras maiúsculas
    tamanho = len(sequencia)
    gc_percent = GC(sequencia)
    tm = get_tm_avg(sequencia)
    return tamanho, gc_percent, tm

# Aplicando as funções para calcular as colunas de tamanho, % GC e Tm
df['tamanho'], df['% GC'], df['Tm'] = zip(*df['nuc'].apply(calcular_propriedades_primers))

# Definindo os critérios de filtro
tamanho_min, tamanho_max = 19, 30
gc_min, gc_max = 45, 55
tm_min, tm_max = 55, 68

# Aplicando os filtros, retirei o filtro de GC
filtered_df = df[
    (df['tamanho'].between(tamanho_min, tamanho_max)) &
    (df['Tm'].between(tm_min, tm_max))
]

not_passed_df = df[~df.index.isin(filtered_df.index)]

print(list(filtered_df['id']))




