import sys
import pandas as pd
import matplotlib.pyplot as plt

column_names=['chr', 'readname', 'unit', 'meth_pct']

df = pd.read_csv(sys.argv[1], sep='\t', header=None, usecols=[0, 1, 2, 3], names=column_names, index_col=None)
print(df.head())
plt.hist(df['meth_pct'], bins=50, color='blue', edgecolor='black')
plt.xlabel('Methylation %')
plt.ylabel('Frequency')
plt.title('Distribution of unit methylation %')
plt.savefig(sys.argv[2])
plt.show()
