import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

from pathlib import Path
from scipy.stats import ks_2samp

import plotnine
from plotnine import *


# read in data
data_dir = Path('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/pe_analysis/data')
output_dir = Path('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/pe_analysis/plots')

# read in optimization data
opt_concs = pd.read_csv(data_dir / 'optimization_data_reformatted_20200925.csv', header=0)
opt_concs.rename(columns={'Sample_ID': 'SID', 'Dataset': 'Phase', 'Perc_free_PlGF' : 'PlGF(%f)'}, inplace=True)
opt_concs = opt_concs.loc[opt_concs['Visit'] == 'Baseline']

# read in PRO-129 clinical verification data
ver_concs = pd.read_csv(data_dir / 'PRO-129_clinical_verification_data_merged_20200925.csv', header=0)
ver_concs.rename(columns={'Dataset': 'Phase', 'Perc_free_PlGF' : 'PlGF(%f)'}, inplace=True)
ver_concs = ver_concs.loc[ver_concs['Visit'] == 'Baseline']

# summary statistics
print(opt_concs.describe())
print(ver_concs.describe())

# merge into one data set, reformat for easy plotting
other = ['SID', 'Phase', 'Class']
markers = ['CD274', 'DCN', 'ENG', 'FGF21', 'KIM1', 'PlGFdiss', 'PlGFfree', 'sFlt.1', 'PlGF(%f)']

data_set = pd.concat([opt_concs[markers + other], ver_concs[markers + other]], axis=0, join='inner')
data = pd.melt(data_set, id_vars=other, value_vars=markers, var_name='Markers', value_name='Conc')
data['Phase-Class'] = data['Phase'] + '-' + data['Class']

#data.to_csv('test.tsv', sep='\t', index=None)

di = {
    'Optimization-Non-PE': 'Opt-NonPE', 
    'Optimization-PE': 'Opt-PE', 
    'PRO-129-PE': 'PRO-129-PE', 
    'PRO-129-Non-PE': 'PRO-129-NonPE'
    }
data = data.replace({'Phase-Class': di})
cat_type = CategoricalDtype(categories=['Opt-NonPE', 'PRO-129-NonPE', 'Opt-PE', 'PRO-129-PE'], ordered=True)
data['Phase-Class'] = data['Phase-Class'].astype(cat_type)

#colors
class_cols = {'Non-PE': '#3182bd', 'PE': '#e6550d'}

# box plots
markers2 = ['CD274', 'DCN', 'ENG', 'FGF21', 'KIM1', 'PlGFdiss', 'PlGFfree', 'sFlt.1',]

p1 = ggplot(data.loc[data['Markers'].isin(markers2)], aes(x='Phase-Class', y='Conc', color='Class'))\
    + geom_boxplot(outlier_shape='') \
    + labs(x='Phase', y='Concentration (pg/mL)', title='Protein Marker Concentrations') \
    + geom_point(alpha=0.2, position=position_jitter(height=0, width=0.1), size=0.5)\
    + facet_wrap('Markers', scales='free_y') \
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=4), axis_text_y=element_text(size=4), subplots_adjust={'wspace': 0.5}) \
    + plotnine.scales.scale_y_log10() \
    + scale_color_manual(values=class_cols) \
    + scale_x_discrete(labels=['Opt', 'PRO-129', 'Opt', 'PRO-129'])

p1.save(output_dir / 'boxplot_markers.png', format='png', dpi=500)

p2 = ggplot(data.loc[data['Markers'] == 'PlGF(%f)'], aes(x='Phase-Class', y='Conc', color='Class'))\
    + geom_boxplot(outlier_shape='') \
    + labs(x='Phase', y='Percentage', title='PlGF(%f)') \
    + geom_point(alpha=0.3, position=position_jitter(height=0, width=0.1))\
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=6), axis_text_y=element_text(size=6)) \
    + scale_color_manual(values=class_cols) \
    + scale_x_discrete(labels=['Opt', 'PRO-129', 'Opt', 'PRO-129'])

p2.save(output_dir / 'boxplot_plgf.png', format='png', dpi=500)

p3 = ggplot(data.loc[data['Markers'] == 'PlGF(%f)'], aes(x='Phase-Class', y='Conc', color='Class'))\
    + geom_boxplot(outlier_shape='') \
    + labs(x='Phase', y='Percentage', title='PlGF(%f)') \
    + geom_point(alpha=0.3, position=position_jitter(height=0, width=0.1))\
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=6), axis_text_y=element_text(size=6)) \
    + plotnine.scales.scale_y_log10() \
    + scale_color_manual(values=class_cols) \
    + scale_x_discrete(labels=['Opt', 'PRO-129', 'Opt', 'PRO-129'])

p3.save(output_dir / 'boxplot_plgf_log.png', format='png', dpi=500)

# only Pro 129 (clinical verification)
p4 = ggplot(data.loc[data['Markers'].isin(markers2) & (data['Phase'] == 'PRO-129')], aes(x='Phase-Class', y='Conc', color='Class'))\
    + geom_boxplot(outlier_shape='') \
    + labs(x='Class', y='Concentration (pg/mL)', title='Clinical Verification (PRO-129): Protein Marker Concentrations') \
    + geom_point(alpha=0.2, position=position_jitter(height=0, width=0.1), size=0.5)\
    + facet_wrap('Markers', scales='free_y') \
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=4), axis_text_y=element_text(size=4), subplots_adjust={'wspace': 0.5}) \
    + plotnine.scales.scale_y_log10() \
    + scale_color_manual(values=class_cols) \
    + scale_x_discrete(labels=['', ''])

p4.save(output_dir / 'boxplot_markers_cv.png', format='png', dpi=500)

p5 = ggplot(data.loc[(data['Markers'] == 'PlGF(%f)') & (data['Phase'] == 'PRO-129')], aes(x='Phase-Class', y='Conc', color='Class'))\
    + geom_boxplot(outlier_shape='') \
    + labs(x='Class', y='Percentage', title='Clinical Verification (PRO-129): PlGF(%f)') \
    + geom_point(alpha=0.3, position=position_jitter(height=0, width=0.1))\
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=6), axis_text_y=element_text(size=6)) \
    + plotnine.scales.scale_y_log10() \
    + scale_color_manual(values=class_cols) \
    + scale_x_discrete(labels=['', ''])

p5.save(output_dir / 'boxplot_plgf_log_cv.png', format='png', dpi=500)


# density plots
phases = {'Optimization': 'opt', 'PRO-129': 'ver'}
classes = {'Non-PE' : 'nonpe', 'PE': 'pe'}

for s in phases.keys():
     p6 = ggplot(data.loc[(data['Phase'] == s) & data['Markers'].isin(markers2)], aes(x='Conc', color='Class', fill='Class'))\
         + geom_density(alpha=0.1) \
         + labs(x='Concentration (pg/mL)', y='Density', title=f'{s}: Protein Marker Concentration Distributions') \
         + facet_wrap('Markers', scales='free') \
         + theme_bw()\
         + theme(axis_text_x=element_text(angle=90, hjust=1, size=6), axis_text_y=element_text(size=6), subplots_adjust={'wspace': 0.25, 'hspace': 0.6}) \
         + plotnine.scales.scale_x_log10() \
         + scale_color_manual(values=class_cols) \
         + scale_fill_manual(values=class_cols) \

     p6.save(output_dir / f'density1_{phases[s]}_class.png', format='png', dpi=500)

phase_cols = [{'Optimization': '#3182bd', 'PRO-129': '#31bd3f'}, {'Optimization': '#e6e20d', 'PRO-129': '#e6550d'}]
i = 0
for s in classes.keys():
     p7 = ggplot(data.loc[(data['Class'] == s) & data['Markers'].isin(markers2)], aes(x='Conc', color='Phase', fill='Phase'))\
         + geom_density(alpha=0.1) \
         + labs(x='Concentration (pg/mL)', y='Density', title=f'{s}: Protein Marker Concentration Distributions') \
         + facet_wrap('Markers', scales='free') \
         + theme_bw()\
         + theme(axis_text_x=element_text(angle=90, hjust=1, size=6), axis_text_y=element_text(size=6), subplots_adjust={'wspace': 0.25, 'hspace': 0.6}) \
         + plotnine.scales.scale_x_log10() \
         + scale_color_manual(values=phase_cols[i]) \
         + scale_fill_manual(values=phase_cols[i]) \
   
     p7.save(output_dir / f'density1_{classes[s]}_phase.png', format='png', dpi=500)
     i += 1


p8 = ggplot(data.loc[data['Markers'] == 'PlGF(%f)'], aes(x='Conc', color='Class', fill='Class'))\
    + geom_density(alpha=0.1) \
    + labs(x='Percentage', y='Density', title='PlGF(%f)') \
    + facet_wrap('Phase') \
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=4), axis_text_y=element_text(size=4), subplots_adjust={'wspace': 0.25, 'hspace': 0.6}) \
    + plotnine.scales.scale_x_log10() \
    + scale_color_manual(values=class_cols) \
    + scale_fill_manual(values=class_cols) \

p8.save(output_dir / 'density2_plgf_percent.png', format='png', dpi=500)

p8 = ggplot(data.loc[data['Markers'] == 'PlGF(%f)'], aes(x='Conc', color='Phase', fill='Phase'))\
    + geom_density(alpha=0.1) \
    + labs(x='Percentage', y='Density', title='PlGF(%f)') \
    + facet_wrap('Class', scales='free') \
    + theme_bw()\
    + theme(axis_text_x=element_text(angle=90, hjust=1, size=4), axis_text_y=element_text(size=4), subplots_adjust={'wspace': 0.25, 'hspace': 0.6}) \
    + plotnine.scales.scale_x_log10() \
    + scale_color_manual(values=phase_cols[0]) \
    + scale_fill_manual(values=phase_cols[0]) \

p8.save(output_dir / 'density2_plgf_percent_phase.png', format='png', dpi=500)


# KS - Test

pe_stat = []
pe_p = []
nonpe_stat = []
nonpe_p = []

opt_stat = []
opt_p = []
ver_stat = []
ver_p = []

ks_result_class = pd.DataFrame(columns = ['Marker', 'PE-statistic', 'PE-pvalue', 'Non-PE-statistic', 'Non-PE-pvalue',])
ks_result_phase = pd.DataFrame(columns = ['Marker', 'Opt-statistic', 'Opt-pvalue', 'Clin-statistic', 'Clin-pvalue',])

for m in markers:
    print(m)
    opt_nonpe = data.loc[(data['Phase-Class'] == 'Opt-NonPE') & (data['Markers'] == m)]
    ver_nonpe = data.loc[(data['Phase-Class'] == 'PRO-129-NonPE')  & (data['Markers'] == m)]

    opt_pe    = data.loc[(data['Phase-Class'] == 'Opt-PE')  & (data['Markers'] == m)]
    ver_pe    = data.loc[(data['Phase-Class'] == 'PRO-129-PE')  & (data['Markers'] == m)]

    pe_stat.append(ks_2samp(opt_pe['Conc'], ver_pe['Conc']).statistic)
    pe_p.append(ks_2samp(opt_pe['Conc'], ver_pe['Conc']).pvalue)
    nonpe_stat.append(ks_2samp(opt_nonpe['Conc'], ver_nonpe['Conc']).statistic)
    nonpe_p.append(ks_2samp(opt_nonpe['Conc'], ver_nonpe['Conc']).pvalue)

    opt_stat.append(ks_2samp(opt_pe['Conc'], opt_nonpe['Conc']).statistic)
    opt_p.append(ks_2samp(opt_pe['Conc'], opt_nonpe['Conc']).pvalue)

    ver_stat.append(ks_2samp(ver_pe['Conc'], ver_nonpe['Conc']).statistic)
    ver_p.append(ks_2samp(ver_pe['Conc'], ver_nonpe['Conc']).pvalue)


ks_result_class['Marker'] = markers
ks_result_class['PE-statistic'] = pe_stat
ks_result_class['PE-pvalue'] = pe_p
ks_result_class['Non-PE-statistic'] = nonpe_stat
ks_result_class['Non-PE-pvalue'] = nonpe_p

print(ks_result_class)
ks_result_class.to_csv('ks_results_class.tsv', sep='\t', index=None)

ks_result_phase['Marker'] = markers
ks_result_phase['Opt-statistic'] = opt_stat
ks_result_phase['Opt-pvalue'] = opt_p
ks_result_phase['Clin-statistic'] = ver_stat
ks_result_phase['Clin-pvalue'] = ver_p

print(ks_result_phase)
ks_result_phase.to_csv('ks_results_phase.tsv', sep='\t', index=None)

    



