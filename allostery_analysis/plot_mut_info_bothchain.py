import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from util.plot_structure_util import plot_vmd_cylinder_from_inds

system = 'dmi_bothchain'

labels_file = f'protease_residue_labels_{system}.npy'
out = f'{system}.npy'
out_name = f'plots/{system}'
o_type = 'pdf'

labels = np.load(open('output/' + labels_file, "rb"))

out00 = np.load(open('output/' + out, "rb"))

out000 = np.copy(out00)

#out000 -= out000.diagonal() * np.eye(*out000.shape)

out000[out000 < 0.05] = np.nan

minr = labels.min()
maxr = labels.max()

plt.imshow(np.log(out000), origin='lower', interpolation='none', cmap='viridis_r', extent=[minr, maxr, minr, maxr])
cb = plt.colorbar()
plt.xlabel('Residue Number', fontsize=16)
plt.ylabel('Residue Number', fontsize=16)
cb.set_label('Log(MI)', fontsize=16)
plt.tight_layout()
plt.savefig(out_name + '_matrix.' + o_type)
plt.savefig('../DESRES_protease_manuscript/raw_figures/dmi_bothchain_matrix.pdf')
plt.clf()

pdb = '../../DESRES_protease_chainid.pdb'
out000 -= out000.diagonal() * np.eye(*out000.shape)

result = np.transpose(np.where(out000 > 0.33))
print('number of cylinders drawn is: ', result.shape[0])

plot_vmd_cylinder_from_inds(pdb, result, f'plots/{system}_over33', residue=True, color='red', width=5)

result1 = np.transpose(np.where((0.33 > out000) & (out000 > 0.13)))
print('number of cylinders drawn is: ', result1.shape[0])

plot_vmd_cylinder_from_inds(pdb, result1, f'plots/{system}_over13', residue=True, width=5)
