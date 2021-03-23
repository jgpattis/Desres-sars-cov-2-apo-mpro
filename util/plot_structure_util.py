#! /usr/bin/env/ python
# plotting utility

import numpy as np
import mdtraj as md


def plot_vmd_cylinder_from_inds(structure_file, inds, outname, residue=False, color='blue', width=3):
    '''writes a tcl file which can draw cylinders in vmd
    structure file should be pdb
    inds should be 0 indexed atom indicies
    or if residue = True 0 indexed residue indicie
    color should be string of VMD color

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
    
    start = '{'
    end = '}'
    bk = '\n'
    
    f.write('set center {0 0 0}\n')
    f.write(f'draw color {color} {bk}')

    for i in range(len(inds)):
        j = inds[i,0]
        k = inds[i,1]
        if residue == True:
            l = top.select('resid ' + str(j) + ' and name CA')
            m = top.select('resid ' + str(k) + ' and name CA')
            n1 = first_frame[0,l,:] * 10
            o1 = first_frame[0,m,:] * 10
            n = n1[0]
            o = o1[0]
        else:
            n = first_frame[0,j,:] * 10
            o = first_frame[0,k,:] * 10
        f.write(f'graphics top line {start} {n[0]:.4f} {n[1]:.4f} {n[2]:.4f} {end} {start} {o[0]:.4f} {o[1]:.4f} {o[2]:.4f} {end} width {width} style solid {bk}' )

    f.close()

def plot_pymol_cylinder_from_inds(structure_file, inds, outname, residue=False):
    '''writes a python file which can draw cylinders in pymol
    structure file should be pdb
    inds should be 0 indexed atom indicies
    or if residue = True 0 indexed residue indicie

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


    for i in range(len(inds)):
        j = inds[i,0]
        k = inds[i,1]
        if residue == True:
            l = top.select('resid ' + str(j) + ' and name CA')
            m = top.select('resid ' + str(k) + ' and name CA')
            n1 = first_frame[0,l,:] * 10
            o1 = first_frame[0,m,:] * 10
            n = n1[0]
            o = o1[0]
        else:
            n = first_frame[0,j,:] * 10
            o = first_frame[0,k,:] * 10
        bk = '\n'
        f.write(f'obj.extend([CYLINDER, {n[0]:.4f}, {n[1]:.4f}, {n[2]:.4f}, {o[0]:.4f}, {o[1]:.4f}, {o[2]:.4f}, 0.15, 0.3917, 0.3917, 0.9961, 0.3917, 0.3917, 0.9961, ]){bk}' )

    f.write("cmd.load_cgo(obj, 'contacts')")
    f.close()
