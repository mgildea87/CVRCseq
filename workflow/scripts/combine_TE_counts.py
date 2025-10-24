import pandas as pd, yaml

def main():
  print('concatenating TEcount files')
  #get config file sample_file
  with open('config/config.yaml', 'r') as file:
  	config = yaml.safe_load(file)
  workflow = config['workflow']
  sample_file = config['sample_file']
  sample_table = pd.read_table(sample_file)
  combine_TE_counts(sample_table, workflow)
  
def combine_TE_counts(sample_table, workflow):
  df = pd.read_csv('%s/results/TEcount/%s.cntTable' % (workflow, sample_table['Sample_name'][0]), sep='\t')
  sample_1 = sample_table['Sample_name'][0]
  df.rename(columns={df.columns[0]: 'gene'}, inplace=True)
  df.rename(columns={df.columns[1]: sample_1}, inplace=True)
  for sample in sample_table['Sample_name']:
    column_data = []
    with open('%s/results/TEcount/%s.cntTable' % (workflow, sample)) as out:
      for line in out:
        column_data.append(line.split('\t')[1].strip())
      column_data.pop(0)
      out.close()
      df[sample] = column_data
  df.to_csv('%s/results/TEcount/count_table_all.csv' % (workflow), index=False)
        
if __name__ == "__main__":
    main()
