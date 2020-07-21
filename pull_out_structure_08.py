# pulls out structures to get a better interpretation of what the top tica dimentions are
import numpy as np
import mdtraj as md
import pyemma.coordinates as coor
import pyemma.plots as pyemma_plots
import pyemma
from util.util import check_iter_of_sequences, KDTree
import matplotlib.pyplot as plt

def heatmap_label(tica_data_cat=None, cl=None, cl_list=None, ic2=1):
    fig, ax = plt.subplots()
    fig, ax, cb = pyemma_plots.plot_density(tica_data_cat[:,0], tica_data_cat[:,ic2], ax=ax, logscale=True)
    ax.scatter(tica_data_cat[0,0], tica_data_cat[0,ic2], marker='x', s=200, c='k')
    if cl_list == None:
        cl_list = list(range(cl.n_clusters))
    ax.scatter(cl.cluster_centers_[cl_list, 0],cl.cluster_centers_[cl_list, ic2], c='k')
    for i in cl_list:
        ax.annotate(f'{i}', (cl.cluster_centers_[i,0], cl.cluster_centers_[i,ic2]), fontsize=16, weight='bold', textcoords="offset points", xytext=(0,10), ha='center')
    ax.set_xlabel(f'IC 1', fontsize=16)
    ax.set_ylabel(f'IC {ic2 + 1}', fontsize=16)
    cb['cbar'].set_label('Density', fontsize=16)
    fig.tight_layout()
    return fig, ax

def colorTIME_label(tica_data_cat=None, cl=None, cl_list=None, ic2=1, stride=2):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=[12.8, 4.8])
    time = list(range(100000))
    time_us = [i/1000 for i in time]
    mid = int(len(tica_data_cat)/2)
    print(mid)
    cbb = ax1.scatter(tica_data_cat[:mid:stride,0], tica_data_cat[:mid:stride,ic2], c=time_us[::stride])
    ax1.scatter(tica_data_cat[0,0], tica_data_cat[0,ic2], marker='x', s=200, c='k')
    if cl_list == None:
        cl_list = list(range(cl.n_clusters))
    ax1.scatter(cl.cluster_centers_[cl_list, 0],cl.cluster_centers_[cl_list, ic2], c='k')
    for i in cl_list:
        ax1.annotate(f'{i}', (cl.cluster_centers_[i,0], cl.cluster_centers_[i,ic2]), fontsize=16, weight='bold', textcoords="offset points", xytext=(0,10), ha='center')
    ax1.set_xlabel(f'IC 1', fontsize=16)
    ax1.set_ylabel(f'IC {ic2 + 1}', fontsize=16)
    ax1.set_title('Chain A', fontsize=18)
    cb = ax2.scatter(tica_data_cat[mid::stride,0], tica_data_cat[mid::stride,ic2], c=time_us[::stride])
    ax2.scatter(tica_data_cat[0,0], tica_data_cat[0,ic2], marker='x', s=200, c='k')
    ax2.scatter(cl.cluster_centers_[cl_list, 0],cl.cluster_centers_[cl_list, ic2], c='k')
    for i in cl_list:
        ax2.annotate(f'{i}', (cl.cluster_centers_[i,0], cl.cluster_centers_[i,ic2]), fontsize=16, weight='bold', textcoords="offset points", xytext=(0,10), ha='center')
    ax2.set_xlabel(f'IC 1', fontsize=16)
    ax2.set_title('Chain B', fontsize=18)
    cb2 = fig.colorbar(cbb, ax=ax1)
    cb2.set_label(f'Time (\u03BCs)', fontsize=16)
    cb1 = fig.colorbar(cb, ax=ax2)
    cb1.set_label(f'Time (\u03BCs)', fontsize=16)
    fig.tight_layout()
    return fig, (ax1, ax2)

sys = 'fdis'
n_clusters = 50
tica_data = coor.load('tica_data_05/fdis_tica_data.h5')
tica_data_cat = np.concatenate(tica_data)
cl = pyemma.load(f'{sys}_{n_clusters}_mini_cluster_object.h5')

#label_list = [12, 41, 46] 
label_list = [19, 48]
pull_structure = True

fig1, ax1 = heatmap_label(tica_data_cat=tica_data_cat, ic2=2, cl=cl, cl_list=label_list)
fig1.savefig(f'structure_08/{sys}_heatmap_label_IC3_two.pdf')

fig2, ax2 = colorTIME_label(tica_data_cat=tica_data_cat, ic2=2, cl=cl, cl_list=label_list)
fig2.savefig(f'structure_08/{sys}_color_label_ic3_two.pdf')

if pull_structure == True:
    traj_num = [f'{i:04d}' for i in range(100)]
    traj_path = '../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
    traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]
    pdb = '../DESRES_protease_chainid.pdb'
    tree = KDTree(tica_data)
    dist, index_list = tree.query(cl.cluster_centers_[label_list], k=1)
    chain_check = []
    for i,j in zip(index_list, label_list):
        if i[0] >= 100:
            chain_check.append('b')
            chain_b = i[0] - 100
            new_pdb = md.load_frame(traj_list[chain_b], index=i[1], top=pdb) # if using a stride multiply index by stride
            new_pdb_chain = new_pdb.atom_slice(new_pdb.topology.select('chainid == 1'))
            new_pdb_chain.save_pdb(f'structure_08/structure_ic3_chainB_index_{j}.pdb')
        else:
            chain_check.append('a')
            new_pdb = md.load_frame(traj_list[i[0]], index=i[1], top=pdb) # if using a stride multiply index by stride
            new_pdb_chain = new_pdb.atom_slice(new_pdb.topology.select('chainid == 0'))
            new_pdb_chain.save_pdb(f'structure_08/structure_ic3_chainA_index_{j}.pdb')
    if len(set(chain_check)) == 1:
        print(f'All structures are from chain {chain_check[0]}')
    else:
        print('WARNING: structures are from different chains')
    print('indexes are ', index_list)
    print('distances from center are ', dist)
