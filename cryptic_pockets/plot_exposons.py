import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
sys.path.append("..")
from util.util import backup

mi_name = 'sasa_mi.npy'
exp_name = 'exposons.npy'
out_name = 'plots/exposons_try1'
o_type = 'pdf'

mi = np.load(open(mi_name,"rb"))
exp = np.load(open(exp_name,"rb"))

#out000 -= out000.diagonal() * np.eye(*out000.shape)

#lab_arr = obj1.labels + 1
#lab_list = lab_arr.tolist()

mi[mi < 0.005] = np.nan

#minr = lab_arr.min()
#maxr = lab_arr.max()

plt.imshow(mi, origin='lower', interpolation='none', cmap='viridis_r') #extent=[minr, maxr, minr, maxr])
cb = plt.colorbar()
plt.xlabel('Residue Number', fontsize=16)
plt.ylabel('Residue Number', fontsize=16)
cb.set_label('Log(MI)', fontsize=16)
plt.tight_layout()
plt.savefig(out_name + '_matrix.' + o_type)
plt.clf()

plt.plot(exp)
plt.xlabel('Residue Number', fontsize=16)
plt.ylabel('exposon', fontsize=16)
plt.tight_layout()
plt.savefig(out_name + '_exposon.' + o_type)
plt.clf()


def print_vmd_string(inds):
    string = f'(residue {inds[0]}'
    for i in range(1, len(inds)):
        string = string + f' or residue {inds[i]}'
    return string + ') and (sidechain or type CA)'

unique, un_inds, un_counts = np.unique(exp, return_index=True, return_counts=True)
order = np.argsort(un_counts)
name = out_name + '_VMD_selection.txt'
backup(name)
f = open(name,'w')
f.write('#selection strings to highlight exposons in VMD\n')
f.write('\n')

for i in range(2,12):
    inds = np.where(exp == unique[order][-i])[0]
    new_inds = print_vmd_string(inds)
    f.write(new_inds + '\n')
    f.write('\n')
