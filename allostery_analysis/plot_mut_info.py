import numpy as np
import pickle
import matplotlib.pyplot as plt

obj = 'dmi_back_try2_both.pkl'
out = 'dmi_back_output_try2_both.pkl'
out_name = 'dmi_back_both_log'
o_type = 'pdf'

with open(obj, 'rb') as handle:
        obj1 = pickle.load(handle)

with open(out, 'rb') as handle:
        out00 = pickle.load(handle)


out000 = np.copy(out00)

#out000 -= out000.diagonal() * np.eye(*out000.shape)

lab_arr = obj1.labels + 1
lab_list = lab_arr.tolist()

out000[out000 < 0.02] = np.nan

#out000[out000 > 1.5] = 1.5

minr = lab_arr.min()
maxr = lab_arr.max()

plt.imshow(np.log(out000), origin='lower', interpolation='none', extent=[minr, maxr, minr, maxr], cmap='viridis_r')
cb = plt.colorbar()
plt.xlabel('Residue Number', fontsize=16)
plt.ylabel('Residue Number', fontsize=16)
cb.label('Log(MI)', fontsize=16)
plt.savefig(out_name + '_matrix.' + o_type)
plt.clf()
