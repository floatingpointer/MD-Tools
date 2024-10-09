"""
After using CHARMM to prep a PDB, this removes Hydrogens and renames
residues as needed (currently CYM -> CYS) to allow H++ to create a pqr
"""

import numpy as np
import pandas as pd
import sys

pdb_cols = [
    ['RType', (0,6)],
    ['Serial', (6,11)],
    ['AName', (12,16)],
    ['Altname', (16,17)],
    ['Resname', (17,21)],
    ['Chain', (21,22)],
    ['Resnum', (22,26)],
    ['Rescode', (26,28)],
    ['X', (30,38)],
    ['Y', (38,46)],
    ['Z', (46,54)],
    ['Occupancy', (54,60)],
    ['Tempfactor', (60,66)],
    ['Segid', (72,76)],
    ['Element', (76,78)],
    ['Charge', (78,80)]
]

pdb_colnames = np.transpose(np.array(pdb_cols, dtype=object))[0].astype(str).tolist()
pdb_colspecs = np.transpose(np.array(pdb_cols, dtype=object))[1].tolist()

pdb = pd.read_fwf(sys.argv[1], skiprows=3, names=pdb_colnames, colspecs=pdb_colspecs)

# delete hydrogens, fix CYM
for i in range(len(pdb)):
  if str(pdb['Resname'][i])=='CYM':
      pdb.loc[i, 'Resname'] = 'CYS'
  if str(pdb['AName'][i])[0]=='H':
    pdb.drop(i, inplace=True)


# Delete TER, END, and all ZN and Water
pdb.drop(pdb[pdb['RType']=='TER'].index, inplace=True)
pdb.drop(pdb[pdb['RType']=='END'].index, inplace=True)
pdb.drop(pdb[pdb['Resname']=='ZN2'].index, inplace=True)
pdb.drop(pdb[pdb['Resname']=='TIP'].index, inplace=True)

# fix Serial numbering
pdb['Serial'] = list(range(1,1+len(pdb)))

for i, row in pdb.iterrows():
  r = list(row)

  print(
      f'{r[0]:<6}' + # RType
      f'{int(r[1]):>5}  ' + # Serial
      f'{r[2]:<3}' + # AName
      f'{"" if np.isnan(r[3]) else r[3]:1}' + # Altname
      f'{r[4]:<4}' + # Resname 
      f'{r[5]:1}' + # Chain
      f'{int(r[6]):>4d}' + # Resnum
      f'{"" if np.isnan(r[7]) else r[7]:2}  ' + # Rescode
      f'{r[8]:>8.3f}' + # X
      f'{r[9]:>8.3f}' + # Y
      f'{r[10]:>8.3f}' + # Z
      f'{r[11]:>6.2f}' + # Occupancy
      f'{r[12]:>6.2f}      ' + # Tempfactor
      f'{r[13]:<4}' + # Segid
      f'{"" if np.isnan(r[14]) else r[14]}' + # Element
      f'{"" if np.isnan(r[15]) else r[15]}' # Charge
  )
  if r[2][0] not in 'CONS':
    print('ATOM found not C, O, N, or S:', r[2])

print(f'{"TER":<6}{int(pdb.iloc[len(pdb)-1]["Serial"]):>5}      {pdb.iloc[len(pdb)-1]["Resname"]:>3}\nEND')
