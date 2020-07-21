#! /usr/bin/env/ python
# plotting utility

import numpy as np
import mdtraj as md


def plot_vmd_cylinder_from_inds(structure_file, inds, outname, residue=False, matrix=False):
    '''writes a tcl file which can draw cylinders in vmd
    structure file should be pdb
    inds should be 0 indexed atom indicies
    or if residue = True residue sequence (generally 1-based, but depends on topology)
    if matrix=true inds can be a matrix

    instructuns:
    1. open structure_file in vmd
    2. open tk consol: extentions > tk consol
    3. default color is dark blue for other colors type: set color 3
        number is vmd color code
    4. source outname.tcl'''

    t = md.load(structure_file)
    top = t.topology
    
    first_frame = t.xyz

    if outname.endswith('.tcl'):
        f = open(outname,'w')
    else:
        f = open(outname + '.tcl','w')
    f.write('set center {0 0 0}\n')
    f.write('draw color blue\n')
    if matrix == True:
        loop = inds[0]
    else:
        loop = inds

    for i in range(len(loop)):
        if matrix == True:
            j = inds[0][i]
            k = inds[1][i]
        else:
            j = inds[i,0]
            k = inds[i,1]
        if residue == True:
            l = top.select('resSeq ' + str(j) + ' and name CA')
            m = top.select('resSeq ' + str(k) + ' and name CA')
            n = first_frame[0,l,:] * 10
            o = first_frame[0,m,:] * 10
        else:
            n = first_frame[0,j,:] * 10
            o = first_frame[0,k,:] * 10
        start = '{'
        end = '}'
        n1 = '\n'
        f.write(f'graphics top line {start} {n[0]:.4f} {n[1]:.4f} {n[2]:.4f} {end} {start} {o[0]:.4f} {o[1]:.4f} {o[2]:.4f} {end} width 3 style solid {n1}' )

    f.close()

def plot_pymol_cylinder_from_inds(structure_file, inds, outname, residue=False, matrix=False):
    '''writes a python file which can draw cylinders in pymol
    structure file should be pdb
    inds should be 0 indexed atom indicies
    or if residue = True residue sequence (generally 1-based, but depends on topology)
    if matrix=true inds can be a matrix

    instructuns:
    1. open pymol
    2. go to command line
    3. run outname.py'''
    

    t = md.load(structure_file)
    top = t.topology
    
    first_frame = t.xyz
    if outname.endswith('.py'):
        f = open(outname,'w')
    else:
        f = open(outname + '.py','w')
    f.write('from pymol import cmd\n')
    f.write('from pymol.cgo import *\n')
    f.write("cmd.load('" + str(structure_file) + "', 'prot')\n")
    f.write("cmd.show('cartoon')\n")
    f.write('obj=[]\n')

    if matrix == True:
        loop = inds[0]
    else:
        loop = inds

    for i in range(len(loop)):
        if matrix == True:
            j = inds[0][i]
            k = inds[1][i]
        else:
            j = inds[i,0]
            k = inds[i,1]
        if residue == True:
            l = top.select('resSeq ' + str(j) + ' and name CA')
            m = top.select('resSeq ' + str(k) + ' and name CA')
            n = first_frame[0,l,:] * 10
            o = first_frame[0,m,:] * 10
        else:
            n = first_frame[0,j,:] * 10
            o = first_frame[0,k,:] * 10
        n1 = '\n'
        f.write(f'obj.extend([CYLINDER, {n[0]:.4f}, {n[1]:.4f}, {n[2]:.4f}, {o[0]:.4f}, {o[1]:.4f}, {o[2]:.4f}, 0.15, 0.3917, 0.3917, 0.9961, 0.3917, 0.3917, 0.9961, ]){n1}' )

    f.write("cmd.load_cgo(obj, 'contacts')")
    f.close()
