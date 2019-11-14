#!/usr/bin/env python
# -*- coding: utf-8 -*-
#A program to parse user.LP
import itertools, sys, os, subprocess, shutil, glob, numpy as np, re, collections, operator, datetime, optparse, csv
#import Bio
trans = {'ALA':'A','CYS':'C','CYH':'C','CSS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','UNK':'X','MSE':'M'}

atom_list = {
			"ALA":('N  ','C  ','O  ','CA ','CB '), "CYS":('N  ','C  ','O  ','CA ','CB ','SG '),"ASP":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','OD2'),"GLU":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','OE2'),"PHE":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ '),"GLY":('N  ','C  ','O  ','CA '),"HIS":('N  ','C  ','O  ','CA ','CB ','CG ','ND1','CD2','CE1','NE2'),
            "ILE":('N  ','C  ','O  ','CA ','CB ','CG1','CG2','CD1'),"LYS":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','CE ','NZ '),"LEU":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2'),"MET":('N  ','C  ','O  ','CA ','CB ','CG ','SD ','CE '),"ASN":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','ND2'),"PRO":('N  ','C  ','O  ','CA ','CB ','CG ','CD '),"GLN":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','NE2'),"ARG":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','NE ','CZ ','NH1','NH2'),
            "SER":('N  ','C  ','O  ','CA ','CB ','OG '),"THR":('N  ','C  ','O  ','CA ','CB ','OG1','CG2'),"VAL":('N  ','C  ','O  ','CA ','CB ','CG1','CG2'),"TRP":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),"TYR":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ ','OH ')}

#=========PART 3a: SAVE INFORMATION===========
#save information of annotated drug binding sites for Drug ReposER application
#BINDING_INTERFACES_NEW, BINDING_INTERFACES_CLUSTERS, BINDING_INTERFACES_ASSAM, BINDING_INTERFACES_EXACT
def save_output_binding_interfaces():
	csv_reader=csv.reader(open("pdbdescription.csv","r"),delimiter=",")
	rows=[row for row in csv_reader][1:-2]
	rows2=[row[:4]+[";".join(sorted(list(set([z+"("+y+")" for z,y in zip(list(set([n.split()[0] for n in row[4].split(", ")])),row[5].split(", "))])))) if row[4] != "" else "None"]+row[6:] for row in rows]
	drugbank_mappings = [[i.replace("\n","").split("\t")[0]]+[i.replace("\n","").split("\t")[1]+"("+i.replace("\n","").split("\t")[2]+")"] for i in open("/Users/nursyatila/NSAG_PART2/db/pdb_drugbank.txt","r").readlines()]
	#dictionary of drugbank annotation
	dict_drugbank_mappings = collections.defaultdict(list)
	#dictionary of drug indication
	dict_indication_mappings = collections.defaultdict(list)
	for d in drugbank_mappings: 
		dict_drugbank_mappings[d[0]].append(d[1]+"("+d[2]+")")
		dict_indication_mappings[d[0]].append(d[3])
	#dictionary of source organism
	dict_source_mappings = collections.defaultdict(list)
	for s in rows2: dict_source_mappings[s[0].lower()+s[1]].append(s[3])
	#dictionary of pdb-macromolecule mappings
	dict_macromolecule_mappings = collections.defaultdict(list)
	for m in rows2: dict_macromolecule_mappings[m[0].lower()+m[1]].append(m[2])
	#pfam annotation
	pfam_annotation = collections.defaultdict(list)
	for p in rows2: pfam_annotation[p[0].lower()+p[1]].append(p[4])

	csv_writer = csv.writer(open("BINDING_INTERFACES.csv","w",),delimiter=",")
	for test_file in glob.glob("BINDING_SITES/*.pdb"):
		pdbs = test_file.split("/")[-1]
		pdb = pdbs[:4].lower() #pdbid
		dreposed_id = pdbs.replace(".pdb","").upper() #drreposed id
		binding_residues = ";".join(sorted(list(set([line[17:26] for line in open(test_file,"r") if line[:4] == "ATOM"])),key=lambda x:int(x[5:]))) #binding residues in 'ATOM' record
		hetatm_residues = list(set([line[17:26] for line in open(test_file,"r").readlines() if line[:6] == "HETATM"]))[0] #drug molecule in 'HETATM' record
		pdb_ligand_id = hetatm_residues[:3].replace(" ","") #pdb ligand id
		pdbchains0 = sorted(list(set([pdb+line[21:22] for line in open(test_file,"r").readlines()])))
		drugbank_id = ";".join(sorted(list(set(dict_drugbank_mappings[pdb_ligand_id.replace(" ","")])))) #drugbank id
		if drugbank_id == "": drugbank_id = "-"
		indications = ";".join(sorted(list(set(dict_indication_mappings[pdb_ligand_id]))))
		organism = ";".join(sorted(list(set([k for n in [dict_source_mappings[pdbchain] for pdbchain in pdbchains0] for k in n])))).replace(",",";")
		if organism == "": organism = "-"
		macromolecule = ";".join(sorted(list(set([k for n in [dict_macromolecule_mappings[pdbchain] for pdbchain in pdbchains0] for k in n])))).upper()
		if macromolecule == "": macromolecule = "-"
		pfam_id = ";".join(sorted(list(set([k for n in [pfam_annotation[pdbchain] for pdbchain in pdbchains0] for k in n]))))
		compiled_details = [dreposed_id]+[pdb]+[pdb_ligand_id]+[drugbank_id]+[indications]+[organism]+[macromolecule]+[pfam_id]+[binding_residues]+[hetatm_residues] #arrangement in csv file
		#print compiled_details
		csv_writer.writerow(compiled_details)
	csv_writer2 = csv.writer(open("BINDING_INTERFACES_CLUSTERS.csv","w",),delimiter=",")
	for test_file in glob.glob("BINDING_SITES_CLUSTERS/*.pdb"):
		pdbs = test_file.split("/")[-1]
		pdb = pdbs[:4].lower() #pdbid
		dreposed_id = pdbs.replace(".pdb","").upper() #drreposed id
		binding_residues = ";".join(sorted(list(set([line[17:26] for line in open(test_file,"r") if line[:4] == "ATOM"])),key=lambda x:int(x[5:]))) #binding residues in 'ATOM' record
		hetatm_residues = list(set([line[17:26] for line in open(test_file,"r").readlines() if line[:6] == "HETATM"]))[0] #drug molecule in 'HETATM' record
		pdb_ligand_id = hetatm_residues[:3].replace(" ","") #pdb ligand id
		pdbchains0 = sorted(list(set([pdb+line[21:22] for line in open(test_file,"r").readlines()])))
		drugbank_id = ";".join(sorted(list(set(dict_drugbank_mappings[pdb_ligand_id.replace(" ","")])))) #drugbank id
		if drugbank_id == "": drugbank_id = "-"
		indications = ";".join(sorted(list(set(dict_indication_mappings[pdb_ligand_id]))))
		organism = ";".join(sorted(list(set([k for n in [dict_source_mappings[pdbchain] for pdbchain in pdbchains0] for k in n])))).replace(",",";")
		if organism == "": organism = "-"
		macromolecule = ";".join(sorted(list(set([k for n in [dict_macromolecule_mappings[pdbchain] for pdbchain in pdbchains0] for k in n])))).upper()
		if macromolecule == "": macromolecule = "-"
		pfam_id = ";".join(sorted(list(set([k for n in [pfam_annotation[pdbchain] for pdbchain in pdbchains0] for k in n]))))
		compiled_details = [dreposed_id]+[pdb]+[pdb_ligand_id]+[drugbank_id]+[indications]+[organism]+[macromolecule]+[pfam_id]+[binding_residues]+[hetatm_residues] #arrangement in csv file
		#print compiled_details
		csv_writer.writerow(compiled_details)

def change_het_res(hetatm_residues):
	if "from" in hetatm_residues:
		if "site" not in hetatm_residues:
			rr = hetatm_residues
			het_resx = [" ".join([r.split(" from ")[1].replace("s","")[-9:]]+["("+r.split(" from ")[0].replace(" ","")+")" if "-" in r.split(" from ")[0] else "( "+r.split(" from ")[0].replace(" ","")+")" ]) if float(re.findall("[0-9]+[.][0-9]+",r.split(" from ")[0])[0])<5.0 else "None" for r in rr.split(";")]
			if all(item == "None" for item in het_resx) == True: het_res = "None"
			else: het_res = ";".join(het_resx)
		else: het_res = "None"
	else: het_res = "None"
	return het_res

#BINDING_INTERFACES_ASSAM
def save_output_binding_interfaces_assam():
	pdb_dict=collections.defaultdict(list)
	pdb_new_dict=collections.defaultdict(list)
	pdbs_all=[["_".join(i.split("/")[-1].split("_")[:3]).upper()]+[i.split("/")[-1].replace(".pdb","")] for i in glob.glob("renew_bs_finalized/*.pdb")]
	for pdb in pdbs_all: pdb_dict[pdb[0]].append(pdb[1])
	for k in pdb_dict.keys():
		for index, pdb in enumerate(pdb_dict[k]):
			pdb_new_dict[pdb].append(k+"_"+str(index))

	csv_writer=csv.writer(open("BINDING_INTERFACES_ASSAM.csv","w"),delimiter=",")
	all_cs=[cs for cs in glob.glob("/Users/nursyatila/NSAG_PART2/ASS_EXE/drreposer_output/*/*/*_sum_2.csv")]	for cs in all_cs:
		csv_readerx=[pdb_new_dict[row[0]]+row[1:] for row in csv.reader(open(cs,"r"),delimiter=",")]
		#rowsx=[pdb_new_dict[row[0]]+row[1:6]+[";".join([change_het_res(hetatm_residues) for hetatm_residues in row[6].split(";")])]+row[7:] for row in csv_readerx]
		for r in csv_readerx:
			csv_writer.writerow(r)
			
#BINDING_INTERFACES_EXACT
def binding_interfaces_exact():
	csv_readerx=csv.reader(open("BINDING_INTERFACES_CLUSTERS.csv","r"),delimiter=",")
	rowsx=[[row[0]]+[row[8]] for row in csv_readerx]
	dreposer_id_dict=collections.defaultdict(list)
	for r in rowsx: dreposer_id_dict[r[0]].append(len(r[1].split(";")))
	csv_readerxx=csv.reader(open("BINDING_INTERFACES_ASSAM.csv","r"),delimiter=",")
	rowsxx=[row for row in csv_readerxx if int(dreposer_id_dict[row[0]][0])==len(row[5].split(";"))]
	csv_writerxx=csv.writer(open("BINDING_INTERFACES_ASSAM_EXACT.csv","w"),delimiter=",")
	for r in rowsxx:
		csv_writerxx.writerow([str(dreposer_id_dict[r[0]][0])]+r)

#=========PART 3b: SAVE PDB-FORMATTED AND GRAPH PATTERNS FOR SPRITE SEARCHES===========
#rename pdb file containing binding sites
def rename_pdb_dbs():
	pdbs_to_edit = [i for i in glob.glob("BINDING_SITES_CLUSTERS/*.pdb")]
	pdbs_to_edit2 = [[i]+[i.split("/")[-1].split("_")[0].lower()] for i in pdbs_to_edit]
	d_dct = collections.defaultdict(list)
	for n in pdbs_to_edit2: d_dct[n[1]].append(n[0])
	keys = sorted([key for key in d_dct.keys()])
	for key in keys:
		for ren, pdb in enumerate(d_dct[key]):
			output = open("PDB_DBS/"+key+"_0"+str(ren)+".pdb","w")
			pdb_lines = ["REMARK SOURCE "+key[:4]]+open(pdb,"r").readlines()
			for li in pdb_lines: print>>output, li.replace("\n","")
			print pdb, key+"_0"+str(ren)


#add information in pat files
def pik_for_sprite():
	proc = subprocess.Popen('source dreposed_sprite_pat.sh',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	proc.communicate()
	csv_reader=csv.reader(open("BINDING_INTERFACES_CLUSTERS.csv","r"),delimiter=",")
	rows=[row for row in csv_reader]
	dictprotein=collections.defaultdict(list)
	for m in rows: dictprotein[m[0]].append(m[6])
	pdbs_to_edit = [i for i in glob.glob("BINDING_SITES_CLUSTERS/*.pdb")]
	pdbs_to_edit2 = [[i]+[i.split("/")[-1].split("_")[0].lower()] for i in pdbs_to_edit]
	d_dct = collections.defaultdict(list)
	for n in pdbs_to_edit2: d_dct[n[1]].append(n[0])
	keys = sorted([key for key in d_dct.keys()])
	for key in keys:
		for ren, pdb in enumerate(d_dct[key]):
			map_bsclus = pdb.split("/")[-1].replace(".pdb","")
			pikfile = "PIK_DBS/"+key+"_0"+str(ren)+".pat"
			to_add="COMPND source {}   ".format(key[:4])+dictprotein[map_bsclus][0]+" ("+map_bsclus+")"
			print pikfile, to_add
			output = open(pikfile,"w")
			pdb_lines = [to_add]+open(pikfile,"r").readlines()[1:]
			for li in pdb_lines: print>>output, li.replace("\n","")
