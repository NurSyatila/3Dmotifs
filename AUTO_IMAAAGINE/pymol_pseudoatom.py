#!/usr/bin/env python
import glob
import os
import sys
import optparse
from optparse import OptionParser
from pymol import cmd, CmdException
import __main__
__main__.pymol_argv =['pymol','-qc']
from time import sleep
import pymol
pymol.finish_launching()
#usage:pymol -qrc superposed_image.py  -- cluster_* 
#pymol -qrc /Users/nursyatila/assam/ASS_EXE/imaaagine_pymol_pseudoatoms.py -- 2nd {}
strand = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
trans = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','EOF':''}

sidechainatoms = { 
			"A":["ca","cb"],
			"C":["ca","sg"],
			"D":["cb","od1+od2"],
			"E":["cg","oe1+oe2"],
			"F":["cg","cz"],
			"G":["ca","n+c"],
			"H":["cg","ce1+ne2"],
            "I":["cb","cg2+cd1"],
            "K":["cg","nz"],
            "L":["cb","cd1+cd2"],
            "M":["ca","sd"],
            "N":["cb","od1+nd2"],
            "P":["ca","cd+cg"],
            "Q":["cg","oe1+ne2"],
            "R":["cd","nh1+nh2"],
            "S":["ca","og"],
            "T":["ca","og1+cg2"],
            "V":["ca","cg1+cg2"],
            "W":["cd1","cz2+cz3"],
            "Y":["cg","cz"]}

def remred():
	pdb_list = glob.glob("*.pdb")
	for pdb in pdb_list:
		cmd.load(pdb)
	cmd.run("/usr/local/bin/all_against_all.py")
	cmd.do("align_all_to_all(selection='all',cutoff=0,cycles=0,full_matrix=1,method='align')")

def generate_pseudo_pdbs(): #three representative atom
	pdb_list = glob.glob("*.pdb")
	for pdbx in pdb_list: 
		cmd.load(pdbx)
		pdb_lines = [line[21:22]+trans[line[17:20]]+line[22:26].replace(" ","") for line in open(pdbx,"r").readlines() if line[:4] == "ATOM" and line[13:15]=="CA"]	
		for x in pdb_lines:
			cmd.pseudoatom(pdbx.replace(".pdb",""),selection="{}/{}/{}".format(x[0],x[2:],sidechainatoms[x[1]][0]),name="PS1",resn=strand[x[1]],resi=x[2:],chain=x[0],hetatm=0,color="gray")
			cmd.pseudoatom(pdbx.replace(".pdb",""),selection="{}/{}/{}".format(x[0],x[2:],sidechainatoms[x[1]][1]),name="PS2",resn=strand[x[1]],resi=x[2:],chain=x[0],hetatm=0,color="gray")
			cmd.pseudoatom(pdbx.replace(".pdb",""),selection="{}/{}/{}".format(x[0],x[2:],sidechainatoms[x[1]][0]+"+"+sidechainatoms[x[1]][1]),name="PS3",resn=strand[x[1]],resi=x[2:],chain=x[0],hetatm=0,color="gray")
		cmd.save(filename="{}.pdb".format(pdbx.replace(".pdb","")),selection=pdbx.replace(".pdb",""),format="pdb")
		cmd.delete(pdbx.replace(".pdb",""))

def images(B):
	#cmd.do("cd {}".format(B))
	pdb_list = glob.glob("*.pdb")
	for pdb in pdb_list:
		cmd.load(pdb)
	cmd.run("/usr/local/bin/all_against_all.py")
	cmd.select("all_ps","name PS1+PS2+PS3")
	cmd.do("align_all_to_all(selection='all_ps',cutoff=0,cycles=0,full_matrix=1,method='align')")
	cmd.remove("het")
	cmd.remove("all_ps")
	cmd.show("sticks","all")
	cmd.hide("sticks","all_ps")
	cmd.hide("sphere","all_ps")
	cmd.set("sphere_scale",0.30,"all_ps")
	cmd.orient("all")
	cmd.multisave("{}.pdb".format(B))
	cmd.set("ray_opaque_background", 0)
	cmd.ray(2400)
	cmd.set("antialias", 1)
	cmd.png("{}.png".format(B),1200,1200,300)
	cmd.quit()


if sys.argv[1] == "remred": 
	remred()

if sys.argv[1] == "pseudo": 
	generate_pseudo_pdbs()

if sys.argv[1] == "superposition":
	B=sys.argv[2]
	images(B)


"""

if len(sys.argv[1:]) == 1:
	generate_pseudo_pdbs()
if len(sys.argv[1:]) == 2:
	images(sys.argv[2])
"""
