import mdtraj as md
import numpy as np
import pickle
import os
import enspara.geometry.pockets as poc

#traj_num = [ '0000', '0020', '0040', '0060', '0080']
traj_num = [f'{i:04d}' for i in range(100)]
traj_path = '../../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]
pdb = '../../DESRES_protease_chainid.pdb'

traj = md.load(traj_list, stride=5, top=pdb)

pockets = poc.get_pockets(traj)

save_num = 10

for i in range(save_num):
    pockets[i].save(f'pockets/pockets{i}.pdb')
