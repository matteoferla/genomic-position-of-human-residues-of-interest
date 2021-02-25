## PTMs

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
    if 'PSP_modified_residues' not in p.features:
        continue
    for feat in p.features['PSP_modified_residues']:
        data = data.append({**{key: getattr(p, key) for key in keepers},
                            **{'residue_idx': feat['residue_index'],
                               'residue_name': feat['from_residue'],
                               'ptm': feat['ptm']}
                           }, ignore_index=True)
# adding pLI
# Variant ID	Gene	Verdict
pLI = pd.read_csv('gnomad.v2.1.1.lof_metrics.by_gene.txt', sep='\t')
ENST2pLI = dict(zip(pLI.transcript.values, pLI.pLI.values))
def get_pLI(ENST): return ENST2pLI[ENST] if ENST in ENST2pLI else float('nan')
data['pLI'] = data.ENST.apply(get_pLI)
data = data.sort_values('pLI', ascending=False)
data.to_csv('ptms.csv')
```
