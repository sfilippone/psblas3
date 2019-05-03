import pandas as pd
import numpy as np

filename = 'cz308.output'

f = pd.read_csv(filename, sep=';')

swap_di = {16:'persistant' , 32: 'nonpersistant'}

q = f[['np','ave_time','swap_mode']]
q['swap_mode'].replace(swap_di, inplace=True)
