import numpy as np
import matplotlib.pyplot as plt

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
plt.clf()
