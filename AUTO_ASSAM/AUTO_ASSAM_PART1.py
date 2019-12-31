#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A script used for analysis of sub-structures in PDB using ASSAM program.
"""
import itertools, sys, os, subprocess, shutil, glob, numpy as np, re, collections, operator, datetime, optparse, csv
#import Bio
trans = {'ALA':'A','CYS':'C','CYH':'C','CSS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','UNK':'X','MSE':'M'}

atom_list = {
			"ALA":('N  ','C  ','O  ','CA ','CB '), "CYS":('N  ','C  ','O  ','CA ','CB ','SG '),"ASP":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','OD2'),"GLU":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','OE2'),"PHE":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ '),"GLY":('N  ','C  ','O  ','CA '),"HIS":('N  ','C  ','O  ','CA ','CB ','CG ','ND1','CD2','CE1','NE2'),
            "ILE":('N  ','C  ','O  ','CA ','CB ','CG1','CG2','CD1'),"LYS":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','CE ','NZ '),"LEU":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2'),"MET":('N  ','C  ','O  ','CA ','CB ','CG ','SD ','CE '),"ASN":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','ND2'),"PRO":('N  ','C  ','O  ','CA ','CB ','CG ','CD '),"GLN":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','NE2'),"ARG":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','NE ','CZ ','NH1','NH2'),
            "SER":('N  ','C  ','O  ','CA ','CB ','OG '),"THR":('N  ','C  ','O  ','CA ','CB ','OG1','CG2'),"VAL":('N  ','C  ','O  ','CA ','CB ','CG1','CG2'),"TRP":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),"TYR":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ ','OH ')}
#=========PART 1: GET DRUG BINDING SITES===========
#1 - Download PDB from Advanced Search (TabularResult.csv
#2 - Get binding sites from protein-drug complexes
def set_pdb_complexes(complex_dir):
	drug_file = csv.reader(open("drugTable.csv","r"),delimiter=",")
	drug_list = [i for i in drug_file][1:-1]
	drug_list2 = sorted(list(set([n for k in [i[4].split(",") for i in drug_list] for n in k]))) #pdbligandid
	remlist=["GOL","EOH"]
	drug_list3 = [" "*(3-len(i))+i for i in drug_list2 if i not in remlist]
	output = open("pair_drug.txt","w")
	for pdb in [f for f in glob.glob("{}/*.pdb".format(complex_dir))]:
		print "....pdb {}".format(pdb)
		acrs = sorted(list(set([line[17:22] for line in open(pdb,"r").readlines() if line[:6] == "HETATM" and line[17:20] in drug_list3])))
		if acrs != []:
			for acr in acrs:
				acrx = acr[:3]
				pdbx = pdb.replace(complex_dir+"/","").replace(".pdb","").upper()+"_"+acr[3:].replace(" ","")
				print>>output, pdb.replace(complex_dir+"/","").replace(".pdb","").upper()+"_"+acr[3:].replace(" ","")+"\t"+acrx

def unique_items(ls):
	seen = set()
	result = []
	for item in ls:
		if len(item) == 2:
			it = item[1]
			if it not in seen:
				result.append(item[0])
				seen.add(it)
		else:
			result.append(item[0])
	return result
	
#create a dictionary for each drug-protein pairs
def set_representative():
	pair_drug = [i.replace("\n","").split("\t") for i in open("pair_drug.txt","r").readlines()]
	pair_dict = collections.defaultdict(list)
	for p in pair_drug: pair_dict[" "*(3-len(p[1]))+p[1]].append(p[0][:4].lower())
	pair_dict = {k:sorted(list(set(v))) for k,v in pair_dict.iteritems()}
	return pair_dict

side_chain_centroid = {'ALA':['CB '], 'CYS':['CB ','SG '], 'ASP':['CB ','CG ','OD1','OD2'],
	 'GLU':['CB ','CG ','CD ','OE1','OE2'], 'PHE':['CB ','CG ','CD1 ','CD2 ','CE1 ','CE2 ','CZ '], 'GLY':['CA '],
	 'HIS':['CB ','CG ','ND1','CD2','CE1','NE2'], 'ILE':['CB ','CG1','CG2','CD1'], 'LYS':['CB ','CG ','CD ','CE ','NZ '],
	 'LEU':['CB ','CG ','CD1','CD2'], 'MET':['CB ','CG ','SD ','CE '], 'ASN':['CB ','CG ','OD1','ND2'],
	 'PRO':['CB ','CG ','CD '], 'GLN':['CB ','CG ','CD ','OE1','NE2'], 'ARG':['CB ','CG ','CD ','NE ','CZ ','NH1','NH2'],
	 'SER':['CB ','OG '], 'THR':['CB ','OG1','CG2'], 'VAL':['CB ','CG1','CG2'],
	 'TRP':['CA ','CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'], 'TYR':['CB ','CG ','CD1','CD2','CE1 ','CE2 ','CZ ','OH ']}

#calculation of distance between coordinates of x,y,z
def get_distance(a,b):
    dist = ((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)**0.5
    return dist

#get binding sites for all pdb-drug complexes
def get_binding_sites_pdb(complex_dir):
	pair_drug_dict_two = set_representative()
	for drug in sorted(pair_drug_dict_two.keys()):
		list_pdbs = pair_drug_dict_two[drug]
		for pdb in list_pdbs:
			print "getting interfaces from ", drug, pdb
			test_pdb = complex_dir+"/"+pdb+".pdb"
			sc_coord = []
			het_coord = []
			pdb_lines = [i.replace("\n","") for i in open(test_pdb,"r").readlines() if i[:4] == "ATOM" ]+[i.replace("\n","") for i in open(test_pdb,"r").readlines() if i[:6] == "HETATM"]
			for line in pdb_lines:
				if line[:6] == "HETATM" and line[17:20] == drug:
					het_coord.append([line[13:26], float(line[30:38]), float(line[38:46]),float(line[46:54])])
				for sc in side_chain_centroid.keys():
					if line[:4] == "ATOM" and line[17:20] == sc and line[13:16] in side_chain_centroid[sc]:
						sc_coord.append([line[13:26], float(line[30:38]), float(line[38:46]),float(line[46:54])])
			E = collections.defaultdict(list)
			for h in het_coord: E[h[0][4:]].append(h)
			hetxx = [list(i) for i in E.items()]
			bs_residues = []
			for het in hetxx:
				combinations_aa_het = [i for i in itertools.product(sc_coord,het[1])]
				ligand_binding_residues = []
				for l in combinations_aa_het:
					if get_distance(l[0],l[1])<4.0: ########
						ligand_binding_residues.append(l[0][0][4:])
				ligand_binding_residues2 = sorted(list(collections.OrderedDict.fromkeys(list(set(ligand_binding_residues)))),key=lambda x:int(x[5:]))
				if len(ligand_binding_residues2)>2:
					bs_residues.append([het[0]]+ligand_binding_residues2)
			for per_set in bs_residues:
				sc_coord2 = []
				Ksc2 = []
				scfound = set()
				output_dir = "renew_bs/"
				if not os.path.exists(output_dir): os.mkdir(output_dir)
				output = open(output_dir+pdb+"_"+per_set[0][4:5]+"_"+per_set[0].replace(" ","")+"_"+str(len(per_set[1:]))+".pdb","w")#'-'.join([k[5:].strip(" ") for k in per_set[1:]])+".pdb","w")
				atom_list = {
		    		"ALA":('N  ','C  ','O  ','CA ','CB '), "CYS":('N  ','C  ','O  ','CA ','CB ','SG '),"ASP":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','OD2'),"GLU":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','OE2'),"PHE":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ '),"GLY":('N  ','C  ','O  ','CA '),"HIS":('N  ','C  ','O  ','CA ','CB ','CG ','ND1','CD2','CE1','NE2'),
            		"ILE":('N  ','C  ','O  ','CA ','CB ','CG1','CG2','CD1'),"LYS":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','CE ','NZ '),"LEU":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2'),"MET":('N  ','C  ','O  ','CA ','CB ','CG ','SD ','CE '),"ASN":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','ND2'),"PRO":('N  ','C  ','O  ','CA ','CB ','CG ','CD '),"GLN":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','NE2'),"ARG":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','NE ','CZ ','NH1','NH2'),
            		"SER":('N  ','C  ','O  ','CA ','CB ','OG '),"THR":('N  ','C  ','O  ','CA ','CB ','OG1','CG2'),"VAL":('N  ','C  ','O  ','CA ','CB ','CG1','CG2'),"TRP":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),"TYR":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ ','OH ')}
				atoms = [atom_list[u[:3]] for u in per_set[1:]]
				V = [i for k in atoms for i in k]
				with open(test_pdb,"r") as pdb_input2:
					pdb_lines = pdb_input2.readlines()
					t = [line.replace("\n","") for line in pdb_lines if line[:4] =="ATOM" and any(item in line[17:26] for item in per_set[1:]) and any(thing in line[13:16] for thing in V)]
					for i in t: sc_coord2.append(i)
					s = [line.replace("\n","") for line in pdb_lines if line[:6] == "HETATM" and line[17:26] == per_set[0]]
					for r in s: sc_coord2.append(r)
				Ksc1 = [i[:16]+i[16:17].replace(i[16:17]," ")+i[17:] for i in sc_coord2]
				for item in Ksc1:
					n = item[13:26]
					if n not in scfound:
						Ksc2.append(item)
						scfound.add(n)
				for i in Ksc2:
					print>>output, i.replace("\n","")

#group binding sites by residue number into different directory
def sorted_by_residue_number():
	pdbs = [f for f in glob.glob("renew_bs/*.pdb")]
	os.mkdir("renew_bs_sorted/")
	for pat in pdbs:
		new_dir = "renew_bs_sorted/{}/".format(str(int(pat.replace(".pdb","").split("_")[-1])))
		if not os.path.exists(new_dir): os.mkdir(new_dir)
		shutil.copy(pat, pat.replace("renew_bs/",new_dir))

#split binding site if residues more than 12 -> first by chain, and then cluster by distance between residues
def split_clusters_by_chains():
	for num in glob.glob("renew_bs_sorted/*"):
		if int(num.split("/")[-1])>12:
			for pdb_name in glob.glob(num+"/*.pdb"):
				print pdb_name
				chain_dict = collections.defaultdict(list)
				pdb_lines = [line.replace("\n","") for line in open(pdb_name,"r") if line[:4]=="ATOM"]
				pdb_hets = [line.replace("\n","") for line in open(pdb_name,"r") if line[:6]=="HETATM"]
				for line in pdb_lines:
					chain_dict[line[20:22]].append(line)
				chain_dict = {k:v+pdb_hets for k,v in chain_dict.items()}
				for chain in chain_dict.keys():
					len_res = sorted(list(set([line[17:26] for line in chain_dict[chain] if line[:4]=="ATOM"])))
					output = open(num+"/"+pdb_name.split("/")[-1].replace(".pdb","_"+chain.replace(" ","")+"_"+str(len(len_res))+".pdb"),"w")
					for x in chain_dict[chain]: print>>output, x
				os.remove(pdb_name)

# move to new directory
def relocate_dict_clusters():
	for num in glob.glob("renew_bs_sorted/*"):
		if int(num.split("/")[-1])>12:
			for pdb_name in glob.glob(num+"/*.pdb"):
				print pdb_name
				pdb_lines = sorted(list(set([line[17:26] for line in open(pdb_name,"r") if line[:4]=="ATOM"])))
				if len(pdb_lines)>2:
					dir_res_num="renew_bs_sorted/"+str(len(pdb_lines))+"/"
					if not os.path.exists(dir_res_num): os.mkdir(dir_res_num) 
					new_pdb="renew_bs_sorted/"+str(len(pdb_lines))+"/"+pdb_name.split("/")[-1]
					shutil.move(pdb_name,new_pdb)
				if len(pdb_lines)<3:
					os.remove(pdb_name)

renumber_bss=[[r]+[12]+[r-12] if r<25 else [r]+[12]+[12]+[r-24] for r in range(13,37)]
renumber_bs=collections.defaultdict(list)
for rn in renumber_bss: renumber_bs[rn[0]].append(rn[1:])

#renumber binding sites/ rename PDB file containing binding sites
def renumber_bs_clusters():
	for num in glob.glob("renew_bs_sorted/*"):
		if int(num.split("/")[-1])>12:
			for pdb_name in glob.glob(num+"/*.pdb"):
				print pdb_name
				numres=int(num.split("/")[-1])
				sc_coord = []
				het_coord = []
				pdb_lines = [line.replace("\n","") for line in open(pdb_name,"r") if line[:4]=="ATOM"]
				het_lines = [line.replace("\n","") for line in open(pdb_name,"r") if line[:6]=="HETATM"]
				for line in pdb_lines:
					for sc in side_chain_centroid.keys():
						if line[17:20] == sc and line[13:16] in side_chain_centroid[sc]:
							sc_coord.append([line[13:26], float(line[30:38]), float(line[38:46]),float(line[46:54])])
				for line in het_lines:
					het_coord.append([line[13:26], float(line[30:38]), float(line[38:46]),float(line[46:54])])
				E = collections.defaultdict(list)
				combinations_aa = [i for i in itertools.product(sc_coord,het_coord)]
				aha = [[l[0][0]]+[l[1][0]]+[ "{0:.1f}".format(get_distance(l[0],l[1]))] for l in combinations_aa]
				all_dict = collections.defaultdict(list)
				for ah in aha: all_dict[ah[0][4:]].append(ah[2])
				all_dict = {k:sorted(v,key=lambda x:float(x))[0] for k,v in all_dict.iteritems()}
				if len(renumber_bs[numres][0])==2:
					all_keys = [n[0] for n in sorted([[n]+[all_dict[n]] for n in all_dict.keys()],key=lambda x:float(x[1]),reverse=True)[:renumber_bs[numres][0][0]]]
					output1 = open(pdb_name.replace(".pdb","_{}.pdb".format(renumber_bs[numres][0][0])),"w")
					output2 = open(pdb_name.replace(".pdb","_{}.pdb".format(renumber_bs[numres][0][1])),"w")
					remove_lines = [line.replace("\n","") for line in open(pdb_name,"r") for n in all_keys if line[:4]=="ATOM" and line[17:26]==n]+het_lines
					new_lines = [line for line in pdb_lines if line not in remove_lines]+het_lines
					for y in remove_lines: print>>output1, y
					for x in new_lines: print>>output2, x
					os.remove(pdb_name)
				else:
					all_keys1 = [n[0] for n in sorted([[n]+[all_dict[n]] for n in all_dict.keys()],key=lambda x:float(x[1]),reverse=True)[:renumber_bs[numres][0][0]]]
					all_keys2 = [n[0] for n in sorted([[n]+[all_dict[n]] for n in all_dict.keys()],key=lambda x:float(x[1]),reverse=True)[renumber_bs[numres][0][0]:renumber_bs[numres][0][1]+12]]
					all_keys3 = [n[0] for n in sorted([[n]+[all_dict[n]] for n in all_dict.keys()],key=lambda x:float(x[1]),reverse=True)[renumber_bs[numres][0][1]+12:renumber_bs[numres][0][2]+24]]
					output1 = open(pdb_name.replace(".pdb","_{}.pdb".format(renumber_bs[numres][0][0])),"w")
					output2 = open(pdb_name.replace(".pdb","_12_12.pdb"),"w")
					output3 = open(pdb_name.replace(".pdb","_12_12_{}.pdb".format(renumber_bs[numres][0][2])),"w")
					remove_lines1 = [line.replace("\n","") for line in open(pdb_name,"r") for n in all_keys1 if line[:4]=="ATOM" and line[17:26]==n]+het_lines
					remove_lines2 = [line.replace("\n","") for line in open(pdb_name,"r") for n in all_keys2 if line[:4]=="ATOM" and line[17:26]==n]+het_lines
					remove_lines3 = [line.replace("\n","") for line in open(pdb_name,"r") for n in all_keys3 if line[:4]=="ATOM" and line[17:26]==n]+het_lines
					for y in remove_lines1: print>>output1, y
					for x in remove_lines2: print>>output2, x
					for z in remove_lines3: print>>output3, z
					os.remove(pdb_name)

#merge all binding sites					
def combine_all_bs():
	os.mkdir("BINDING_SITES/")
	for bs in glob.glob("renew_bs_sorted/*/*.pdb"):
		new_bs = "BINDING_SITES/"+bs.split("/")[-1]
		shutil.copy(bs, new_bs)
		
#sort by residue number of ASSAM searches
def sorted_by_residue_number_finalized():
	pdbs = [f for f in glob.glob("renew_bs/*.pdb")]
	for pat in pdbs:
		new_dir = "renew_bs_sorted/{}/".format(str(int(pat.replace(".pdb","").split("_")[-1])))
		if not os.path.exists(new_dir): os.mkdir(new_dir)
		shutil.copy(pat, pat.replace("renew_bs/",new_dir))

"""
complex_dir="pdbs"
set_pdb_complexes(complex_dir)
get_binding_sites_pdb(complex_dir)
sorted_by_residue_number()
split_clusters_by_chains()
relocate_dict_clusters()
renumber_bs_clusters()
relocate_dict_clusters()
combine_all_bs()
"""
