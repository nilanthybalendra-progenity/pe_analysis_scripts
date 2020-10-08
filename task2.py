import numpy as np
import pandas as pd

from pathlib import Path
from scipy import stats

# function to calculate Cohen's d for independent samples
def cohend(d1, d2):
	# calculate the size of samples
	n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
	s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
	# calculate the pooled standard deviation
	s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
	u1, u2 = np.mean(d1), np.mean(d2)
	# calculate the effect size
	return (u1 - u2) / s


data_dir = Path('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/pe_analysis/data')
data = pd.read_csv(data_dir / 'data_reorg.tsv', sep='\t', header=0) # created in task 1

# calculate pearson correlation
print(data.columns)
markers = ['CD274', 'DCN', 'ENG', 'FGF21', 'KIM1', 'PlGFdiss', 'PlGFfree', 'sFlt.1', 'PlGF(%f)']
d_out_opt = []
d_out_ver = []
for m  in markers:
    opt_nonpe = data.loc[(data['Phase-Class'] == 'Opt-NonPE') & (data['Markers'] == m), 'Conc']
    ver_nonpe = data.loc[(data['Phase-Class'] == 'PRO-129-NonPE')  & (data['Markers'] == m), 'Conc']
    print(ver_nonpe)

    opt_pe    = data.loc[(data['Phase-Class'] == 'Opt-PE')  & (data['Markers'] == m), 'Conc']
    ver_pe    = data.loc[(data['Phase-Class'] == 'PRO-129-PE')  & (data['Markers'] == m), 'Conc']
    d_out_opt.append(np.abs(cohend(opt_pe, opt_nonpe)))
    d_out_ver.append(np.abs(cohend(ver_pe, ver_nonpe)))


result = pd.DataFrame()
result['Marker'] = markers
result['Cohen_d_Opt'] = d_out_opt
result['Cohen_d_Ver'] = d_out_ver

result.to_csv('cohen.tsv', sep='\t', index=None)
