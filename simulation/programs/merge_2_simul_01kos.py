#1/bin/usr/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('simul1', help='Directory with simulations for first comparate')
parser.add_argument('simul1_set', help='Which set of simulations for first comparate to consider')
parser.add_argument('simul2', help='Directory with simulations for second comparate')
parser.add_argument('simul2_set', help='Which set of simulations for second comparate to consider')
parser.add_argument('outdir', help='Output directory')
args = parser.parse_args()

proj_dir = '/mnt/c/Users/ksherbina/mcintyre_lab/BASE_2020/simulation'
df_c1_files = pd.DataFrame([args.simul1+'/'+i for i in os.listdir(os.path.expanduser(proj_dir+'/g3_sim_output/'+args.simul1)) if 'set_'+str(args.simul1_set) in i], columns = ['c1'])
df_c1_param = df_c1_files['c1'].str.split('_qsim', 1, expand=True)
df_c1_files['c1_alpha'] = df_c1_param[0].str.split('alpha', expand=True)[1].str.replace('_', 'alpha1_')
df_c1_files['param'] =  df_c1_param[1]

df_c2_files = pd.DataFrame([args.simul2+'/'+i for i in os.listdir(os.path.expanduser(proj_dir+'/g3_sim_output/'+args.simul2)) if 'set_'+str(args.simul2_set) in i], columns = ['c2'])
df_c2_param = df_c2_files['c2'].str.split('_qsim', 1, expand=True)
df_c2_files['c2_alpha'] = df_c2_param[0].str.split('alpha', expand=True)[1].str.replace('_', 'alpha2_')
df_c2_files['param'] =  df_c2_param[1]

df_both_comparates = pd.merge(df_c1_files, df_c2_files, on='param', how='inner')
print(df_both_comparates)

for index, row in df_both_comparates.iterrows():
    df_c1_simul = pd.read_csv(proj_dir+'/g3_sim_output/'+row['c1'])
    df_c2_simul = pd.read_csv(proj_dir+'/g3_sim_output/'+row['c2'])
    df_c2_simul.drop(['comparison', 'FEATURE_ID'], axis=1, inplace=True)
    df_c2_simul.columns = df_c2_simul.columns.str.replace('c1', 'c2')
    filename = (row['c1_alpha'] + '_' + row['c2_alpha'] + '_qsim'+row['param']).replace('qsim_', 'qsim-')
    print(filename)
    pd.concat([df_c1_simul, df_c2_simul], axis=1, sort=False).to_csv(proj_dir+'/g3_sim_output/'+args.outdir+'/'+filename, index=False)
