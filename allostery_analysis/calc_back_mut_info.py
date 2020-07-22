import mdentropy.metrics
import mdtraj as md
import numpy as np
import pickle

#traj_num = [ '0000', '0020', '0040', '0060', '0080']
traj_num = [f'{i:04d}' for i in range(100)]
traj_path = '../../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]
pdb = '../../DESRES_protease_chainid.pdb'

traj = md.load(traj_list, stride=5, top=pdb)

dmi = mdentropy.metrics.DihedralMutualInformation(types=['phi', 'psi'], n_bins=3, method='knn', normed=True)
out = dmi.partial_transform(traj)

print('entropy calulation is complete for backbone dihedral')
print('output is ', out)

p = pickle.HIGHEST_PROTOCOL
with open('dmi_back_normed.pkl','wb') as handle:
        pickle.dump(dmi, handle,protocol=p)

with open('dmi_back_output_normed.pkl','wb') as handle:
        pickle.dump(out, handle,protocol=p)

x = np.copy(out)
np.savetxt('back_dmi_normed.dat', x, fmt='%12.4f')
