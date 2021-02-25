# Genomic position of human residues of interest

> NB. there is no Python module or similar here. It is a file storage repo.

Human post translational modifications etc. in genomic position.

## Getting amino acids of interest

Using the protein module for Michelanglo, (https://github.com/matteoferla/MichelaNGLo-protein-module):

```python
from michelanglo_protein import ProteinCore
import os

ProteinCore.settings.verbose = True #False
ProteinCore.settings.startup(data_folder=os.environ['PROTEIN_DATA'])
```
This module creates and handles pickled version of the Uniprot database enriched with other stuff.
If you are not Matteo and want the data see https://github.com/matteoferla/MichelaNGLo-human-protein-data.
```python
p = ProteinCore(uniprot ='Q99624', taxid=9606).load()
print({'recommended_name': p.recommended_name, 
                      'uniprot': p.uniprot,
                      'gene_name': p.gene_name, 
                      'length': len(p)})
```

    {'recommended_name': 'Sodium-coupled neutral amino acid transporter 3', 'uniprot': 'Q99624', 'gene_name': 'SLC38A3', 'length': 504}

Parsing for something (redacted):
```python
import pandas as pd
import os

offset = 7 # how many residues before and after to search.
keepers = ['gene_name', 'recommended_name', 'uniprot', 'ENSG', 'ENST', 'ENSP']
data = pd.DataFrame([], columns=keepers+['residue_idx', 'aa_to_interface'])
for file in os.listdir(f'{ProteinCore.settings.pickle_folder}/taxid9606'):
    if file[-2:] != '.p':
        continue
    p = ProteinCore(uniprot =file[:-2], taxid=9606).load()
    for feat in p.features['👾👾']:
        data = data.append({**{key: getattr(p, key) for key in keepers},
                            **{'residue_idx': feat['x'],
                               'residue_name': p.sequence[feat['x'] - 1]}
                           }, ignore_index=True)
# adding pLI
# Variant ID	Gene	Verdict
pLI = pd.read_csv('gnomad.v2.1.1.lof_metrics.by_gene.txt', sep='\t')
ENST2pLI = dict(zip(pLI.transcript.values, pLI.pLI.values))
def get_pLI(ENST): return ENST2pLI[ENST] if ENST in ENST2pLI else float('nan')
data['pLI'] = data.ENST.apply(get_pLI)
data = data.sort_values('pLI', ascending=False)
data.to_csv('👾👾.csv')
```

To reverse translate using BackLocate (http://lindenb.github.io/jvarkit/BackLocate.html),
first I made the tab-separate file:

```python
data['mutation'] = data.apply(lambda row: f'{row.residue_name}{row.residue_idx}X'), 1)
mini = data.loc[data.gene_name.apply(len) != 0][['gene_name','mutation', 'ENST']]
mini.to_csv('👾👾.tsv', sep='\t', index=False, header=False)
```
One cannot search by ENST, but one can always filtere out the non-matches.
```bash
cat 👾👾.tsv | java -jar 👾👾👾/backlocate.jar -R 👾👾👾.fna -g 👾👾👾.gtf.gz > backlocated.txt
```

Back in the hygge of Python

```python
backed = pd.read_csv('backlocated.txt', sep='\t')
backed = backed.loc[backed['transcript.AA'] == '👾']\
               .sort_values(['chromosome','index0.in.genomic'])\
               .drop_duplicates(['chromosome','index0.in.genomic'])
backed['genomic'] = backed['index0.in.genomic'].astype(int)
backed['petide'] = backed['petide.pos.1'].astype(int)
backed = backed[['transcript.name', 'petide', 'chromosome','genomic']]
backed.to_csv('👾👾.genomic.csv')
```

Filtering by pLI beforehand is called for, in GEL at least. [pLI.json](pLI.json) are the genes with pLI+pRec ≥ 1.
