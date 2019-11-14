#!/usr/bin/env python
# -*- coding: utf-8 -*-
#A program to conduct automatic search and clustering of amino acid patterns
import itertools, sys, os, subprocess, shutil, glob, re, collections, operator, datetime, optparse, csv, random

import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import fcluster
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from scipy.cluster import hierarchy as hc
from scipy.spatial import distance as dc
from scipy.spatial import distance_matrix
ba_dict=collections.defaultdict(list)
baalls=[i.replace("\n","").split() for i in open("biological_assemblies.txt","r").readlines()]
for b in baalls: ba_dict[b[0]].append(b[1])

trans = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','EOF':''}
strand = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
opt_dict = {'all':'ACDEFGHIKLMNPQRSTVWY','charged-amide':'DEHKRQN','charged-aromatic-amide':'FWYDEHKRQN','charged':'DEHKR','basic':'HKR','acidic':'DE','amide':'QN','aromatic':'FWY','hydrophobic':'LVIAPMFWY','hydroxyl':'CST'}		

#----------------------FIRST OPTIONS: SEARCH FOR HYPOTHETICAL SYMMETRICAL PATTERNS-------------------------
#options for search: 1)R: residue number, 2)N: list of residues for search e.g. DEHRK or -P:exact pattern for search 

trans = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','EOF':''}
strand = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
opt_dict = {'all':'ACDEFGHIKLMNPQRSTVWY','charged-amide':'DEHKRQN','charged-aromatic-amide':'FWYDEHKRQN','charged':'DEHKR','basic':'HKR','acidic':'DE','amide':'QN','aromatic':'FWY','hydrophobic':'LVIAPMFWY','hydroxyl':'CST'}		

def input_files(d,r,combinations): # 02-CREATE INPUT FILES FOR IMAAAGINE SEARCHES
	for comb in combinations:
		#d = 3.5 #3.5+1.5 distance tolerance -> 2-5 -> 3.5 = 1.5-4+1.5 = 5.5
		x = "\n%.2f" % d #required distance
		y = "\n " #leave blank for wildcard distance, or Y for close contact
		if r == 3: distance = "{}".format(x*3)
		if r == 4: distance = "{}{}{}".format(x*2,y*2,x*2)
		if r == 5: distance = "{}{}{}{}{}{}{}".format(x,y*2,x*2,y*2,x,y,x)
		if r == 6: distance = "{}{}{}{}{}{}{}{}{}".format(x,y*3,x*2,y*3,x,y*2,x,y,x)
		if r == 7: distance = "{}{}{}{}{}{}{}{}{}{}{}".format(x,y*4,x*2,y*4,x,y*3,x,y*2,x,y,x)
		if r == 8: distance = "{}{}{}{}{}{}{}{}{}{}{}{}{}".format(x,y*5,x*2,y*5,x,y*4,x,y*3,x,y*2,x,y,x)
	inputs = [''.join(''.join(str(r))+'\n'+''.join(comb)+''.join(distance)+'\n'+">>EOF") for comb in combinations]
	return inputs
	
def imaaagine(r,n_dict,inputs, output_dir): # 03-RUN IMAAAGINE SEQUENTIALLY
	for input in inputs:
		proc = subprocess.Popen('source ass_imaaagine.sh \n source hyp_assam.sh',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		proc.stdin.write(input)
		while True:
			out = proc.stdout.read()
			if out == '' and proc.poll() !=None: break
			if out != '':
				sys.stdout.write(out)
				sys.stdout.flush()
				src = ''.join([x.strip('\n>> .1234567890') for x in input])
				dict = ''.join(trans[(''.join(src))[x:x+3]] for x in range(0, len(''.join(src)), 3))
				infile1 = 'user.LP'
				infile2 = 'user.SUMAj'
				output_directory = "{}{}/{}".format(output_dir,n_dict,r)
				if not os.path.exists(output_directory): os.makedirs(output_directory)
				outfile1 = '{}/{}.LP'.format(output_directory,dict)
				outfile2 = '{}/{}.SUMAj'.format(output_directory,dict)
				shutil.move(infile1, outfile1)
				shutil.move(infile2, outfile2) # 04-SAVE THE OUTPUT FILES CONTAIN PDBID AND MATCHED PATTERN

def remove_no_hits(r,n_dict,output_dir): #06-REMOVE PATTERNS W/O HIT
	os.chdir("{}{}/{}".format(output_dir,n_dict,r))
	no_hits_directory = "no_hits"
	if not os.path.exists(no_hits_directory): os.mkdir(no_hits_directory)
	for filename in glob.glob("*.LP"):
		input = open(filename,"r")
		all = re.findall(r"^(   [PATTERN].*\s     [0-9].*\s.*\s.*\s.*\s.*\s.*\s.*\s.*\szap)",input.read(), re.MULTILINE)
		patterns = collections.OrderedDict.fromkeys([i.replace(i[:54],i[34:38]).replace(i[-4:],"") for i in all]).keys()
		sorted_hits = [[i.split("\n")[0]]+sorted([w[20:21]+w[21:26]+w[10:14]+w[56:60] for w in i.split("\n")[1:]]) for i in patterns]
		non_redundant_hits = sorted(list(set([tuple(i) for i in sorted_hits])))
		if len(non_redundant_hits) == 0:
			os.rename(filename, "no_hits/{}".format(filename))
			os.rename(filename.replace(".LP",".SUMAj"), "no_hits/{}".format(filename.replace(".LP",".SUMAj")))
	os.chdir("/Users/nursyatila/NSAG_PART1/ASS_EXE/")

#----------------------SECOND OPTIONS: PARSE OUTPUT FILES-------------------------

#01-parse Lp file and get patterns from hits
def parse_output(filename):
	input = open(filename,"r")
	all = re.findall(r"^(   [PATTERN].*\s     [0-9].*\s.*\s.*\s.*\s.*\s.*\s.*\s.*\szap)",input.read(), re.MULTILINE)
	patterns = [s.split("\n") for s in collections.OrderedDict.fromkeys([i.replace(i[:54],i[34:38]).replace(i[-4:],"") for i in all]).keys()]
	sorted_hits=[[pat[0][:4]]+sorted([al[10:13]+al[19:25] for al in pat[1:]],key=lambda x:x[4:])+sorted(list(set([al[62:71] for al in pat[1:] if float(re.findall(r"[-+]?[0-9]+[.][0-9]+", al)[0])<4.0 ]))) for pat in patterns if len("".join(sorted(list(set([al[19:21].strip() for al in pat[1:]])))))>1]
	non_redundant_hits = sorted(list(set([tuple(i) for i in sorted_hits])))
	indexed_hits = [i.split("\t") for i in sorted(list(set(["\t".join(x) for x in [[item[0][:4]] + list(i for i in item[1:]) for item in non_redundant_hits]])))]
	return indexed_hits

#02-generate PDB-formatted patterns and generate PDB-formatted patterns with pseudoatoms - will later be used for superposition
def get_coordinates(indexed_hits, filename): # 02- GENERATE TEMPORARY PDB FILE FOR EACH MOTIFS
	for m in indexed_hits:
		pdb = m[0]
		#pdb_filename = m[0]+"_"+'-'.join(i[4:5]+trans[i[:3]][0]+str(int(i[5:10])).replace("-","") for i in m[1:len(filename.replace(".LP",""))+1])+"_"+'-'.join(i[:3].replace(" ","") for i in m[len(filename.replace(".LP",""))+1:])
		pdb_filename = m[0]+"_"+'-'.join(i[4:5]+trans[i[:3]][0]+str(int(i[5:10])).replace("-","") for i in m[1:len(filename.replace(".LP",""))+1])
		ATOM_HETATM = K = V = []
		found = set()
		atom_list = {
			"ALA":('N  ','C  ','O  ','CA ','CB '), "CYS":('N  ','C  ','O  ','CA ','CB ','SG '),"ASP":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','OD2'),"GLU":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','OE2'),"PHE":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ '),"GLY":('N  ','C  ','O  ','CA '),"HIS":('N  ','C  ','O  ','CA ','CB ','CG ','ND1','CD2','CE1','NE2'),
            "ILE":('N  ','C  ','O  ','CA ','CB ','CG1','CG2','CD1'),"LYS":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','CE ','NZ '),"LEU":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2'),"MET":('N  ','C  ','O  ','CA ','CB ','CG ','SD ','CE '),"ASN":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','ND2'),"PRO":('N  ','C  ','O  ','CA ','CB ','CG ','CD '),"GLN":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','NE2'),"ARG":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','NE ','CZ ','NH1','NH2'),
            "SER":('N  ','C  ','O  ','CA ','CB ','OG '),"THR":('N  ','C  ','O  ','CA ','CB ','OG1','CG2'),"VAL":('N  ','C  ','O  ','CA ','CB ','CG1','CG2'),"TRP":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),"TYR":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ ','OH ')}
		atoms = [atom_list[u[:3]] for u in m[1:len(filename.replace(".LP",""))+1]]
		V = [i for k in atoms for i in k]
		#with open("/Users/nursyatila/NSAG_PART1/ba_all/{}.pdb".format(pdb)) as pdb_file:
		with open("ba_all/{}.pdb".format(pdb)) as pdb_file:
			pdb_lines = pdb_file.readlines()
			t = [line[:22]+line[22:26].replace("-"," ")+line[26:27].replace(line[26:27]," ")+line[27:] for line in pdb_lines if line[:4] == "ATOM" and any(item in line[17:26] for item in m[1:len(filename.replace(".LP",""))+1]) and any(thing in line[13:16] for thing in V)]
			for i in t: ATOM_HETATM.append(i)
			het = [line for line in pdb_lines if line[:6] == "HETATM" and line[17:26] in m[len(filename.replace(".LP",""))+1:]]
			for i in het: ATOM_HETATM.append(i)
		for item in ATOM_HETATM:
			n = '{} {}'.format(item[13:16],item[17:26])
			if n not in found:
				K.append(item)
				found.add(n)
		K2 = list(collections.OrderedDict.fromkeys([i[:16]+i[16:17].replace(i[16:17]," ")+i[17:] for i in K]))
		if len(K2)-len(het) == len(V): 
			temp_output = open("{}.pdb".format(pdb_filename),"w")
			for line in K2: 
				temp_output.write(line)
	return V
#03-remove redundant patterns, rmsd=0.0
def remred_one():
	proc = subprocess.Popen("pymol -qrc pymol_pseudoatom.py -- remred > remred.txt", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	proc.communicate()

def remred_two():
	mat="remred.txt"
	allpdbs=sorted([pdb for pdb in glob.glob("*.pdb")])
	xxx="       ".join(open(mat,"r").read().split("       ")[1:])
	rmsdval=[r for r in re.findall(r"\d+\.\d+",xxx)]
	names=[s for s in itertools.combinations(allpdbs,2)]
	macs=[[n[0][0]]+[n[0][1]]+[n[1]] for n in list(zip(names, rmsdval))]
	macs2=[[n[0][0]]+[n[0][1]]+[n[1]] for n in list(zip(names, rmsdval))]+[[n[0][1]]+[n[0][0]]+[n[1]] for n in list(zip(names, rmsdval))]+[[c]+[c]+["0.0"] for c in allpdbs]
	data=pd.DataFrame({"a":[a[0] for a in macs2],"b":[a[1] for a in macs2],"rmsd":[float(a[2]) for a in macs2]})
	data_piv = data.pivot("a", "b", "rmsd")
	piv_arr = data_piv.as_matrix()
	dist_mat = piv_arr + np.transpose(piv_arr)
	Z=hcl.linkage(squareform(dist_mat),method="average")
	max_d=0.5
	clusters=fcluster(Z,max_d,criterion='distance')
	merges_cl=[list(i) for i in zip(clusters,allpdbs)]
	cluster_dict=collections.defaultdict(list)
	for c in merges_cl:
		cluster_dict[c[0]].append(c[1])
	clusters = [i[1] for i in list(cluster_dict.items())]
	for s in clusters:
		if len(s)>1:
			if len(list(set([k[:4] for k in s])))==1:
				for p in s[1:]:
					os.remove(p)

#03-include pseudoatoms in patterns
def generate_pseudoatoms():
	proc = subprocess.Popen("pymol -qrc pymol_pseudoatom.py -- pseudo", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	proc.communicate()

#04-superposition using pseudoatoms and cluster similar patterns of RMSD value lesser than 1.5A
def superposition():
	proc = subprocess.Popen("R < superposition_and_clustering.R --no-save", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	proc.communicate()

def clustering(matrix_file):
	pymol_clusters = []
	input = open(matrix_file,"r")
	f_cluster = [[i.replace("\n","").split()[0]]+[i.replace("\n","").split()[1]] for i in input.readlines()]
	dictionary = collections.defaultdict(list)
	for i in f_cluster: dictionary[i[1]].append(i[0])
	clusters = sorted([[i[0]]+i[1] for i in list(dictionary.items())],key=lambda x:len(x[1:]), reverse=True)
	rmsd_clusters=[i for i in clusters if len(i) > 2 if len(list(set([x[:4] for x in i[1:]])))>1]
	single_clusters = list(set([i[1] for i in clusters if len(i) == 2]+[k for n in [i[1:] for i in clusters if len(list(set([x[:4] for x in i[1:]])))==1] for k in n]))
	for i in single_clusters:
		single_directory = "singletons"
		if not os.path.exists(single_directory): os.mkdir(single_directory)
		infile = '{}'.format(i)
		outfile = '{}/{}'.format(single_directory,i)
		os.rename(infile,outfile)
	for j,h in enumerate(rmsd_clusters):
		dir_name = "cluster_"+str(j)
		pymol_clusters.append(dir_name)
		cluster_members = h[1:]
		output_directory = '{}'.format(dir_name)
		if not os.path.exists(output_directory): os.makedirs(output_directory)
		for m in cluster_members:
			infile = '{}'.format(m)
			outfile = '{}/{}'.format(output_directory,m)
			os.rename(infile,outfile)
	return rmsd_clusters, single_clusters, pymol_clusters

#06- group similar patterns together
def pymol(pymol_clusters): #06- SUPERPOSE PATTERNS IN THE CLUSTER i.e PATTERNS WITH SIMILAR RESIDUE CONFORMATIONS
	if pymol_clusters != []:
		#output_pymol = open("rmsd_clusters.txt","a")
		for i in pymol_clusters:
			proc = subprocess.Popen("pymol -qrc pymol_pseudoatom.py -- superposition {}".format(i), shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
			proc
			output_pymol = open("{}/pymol_matrix.txt".format(i),"w")
			while True:
				out = proc.stdout.read()
				if out == '' and proc.poll() !=None: break
				if out != '': output_pymol.write(out)

#07-generate list of details for hits
def list_details(filename,item,group): #07- EXTRACT PFAM ANNOTATION FOR EACH PDB CHAIN
#macromolecule name mappings
	items = item.split("_")
	#pdbchains0 = ",".join(sorted(list(set([items[0][:4].lower()+i[0] for i in items[1].split("-")]))))
	pdbchainsx = sorted(list(set([items[0][:4].lower()+i[0] for i in items[1].split("-")])))
	pdbchains0 = ",".join(sorted(list(set([ba_dict[c][0] if c in ba_dict.keys() else c for c in pdbchainsx]))))
	csv_reader=csv.reader(open("pdb_description.csv","r"),delimiter=",")
	rows=[row for row in csv_reader]
	dict_macromolecule_mappings = collections.defaultdict(list)
	dict_source_mappings = collections.defaultdict(list)
	dict_classification_mappings = collections.defaultdict(list)
	dict_pfam_mappings = collections.defaultdict(list)
	#dict_pubmed_mappings = collections.defaultdict(list)
	for r in rows: 
		dict_macromolecule_mappings[r[0].lower()+r[1]].append(r[3].upper())
		dict_classification_mappings[r[0].lower()+r[1]].append(r[2].upper())
		dict_source_mappings[r[0].lower()+r[1]].append(r[4])
		#dict_pubmed_mappings[r[0].lower()+r[1]].append(r[5])
		dict_pfam_mappings[r[0].lower()+r[1]].append(r[6])
	macromolecule=";".join(sorted(list(set([k for n in [dict_macromolecule_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n]))))
	source=";".join(sorted(list(set([k for n in [dict_source_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n]))))
	classification=";".join(sorted(list(set([k for n in [dict_classification_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n]))))
	pfamx=";".join(sorted(list(set(";".join(sorted(list(set([k for n in [dict_pfam_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n])))).replace(" ","").split(";")))))
	pfam=",".join(sorted(list(set(pfamx.split(",")))))
	#pfam=";".join(sorted(list(set([k for n in [dict_pfam_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n]))))
	#pubmed=";".join(sorted(list(set([k for n in [dict_pubmed_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n]))))
	if classification == "": classification = "-"
	if macromolecule == "": macromolecule = "-"
	if source == "": source = "-"
	if pfam == "": pfam = "-"	
	if group == "none":
		matched_residues = ";".join([line[17:26] for line in open(item,"r").readlines() if line[:4] == "ATOM" and line[13:15] == "CA"])
		hetatm_residues = ";".join(sorted(list(set([line[17:26] for line in open(item,"r").readlines() if line[:6] == "HETATM"]))))	
	else:
		matched_residues = ";".join([line[17:26] for line in open(group+"/"+item,"r").readlines() if line[:4] == "ATOM" and line[13:15] == "CA"])
		hetatm_residues = ";".join(sorted(list(set([line[17:26] for line in open(group+"/"+item,"r").readlines() if line[:6] == "HETATM"]))))		
	details = [item[:4]]+[classification]+[macromolecule]+[source]+[pfam]+[matched_residues]+[hetatm_residues]
	return details

#08-combine annotation into output file - annotation for hits with singletons or only 1 hit
def no_clusters(filename,indexed_hits):
	N = []
	for m in indexed_hits: 
		item = m[0]+"_"+'-'.join(i[0]+str(int(i[1:5])) for i in m[1:])+'_'+'-'.join(i[-3:].strip(" ") for i in m[1:]).strip(" ")
		N.append(list_details(filename,item,"none"))
	return N

#09-compile annotation from all hits - from clusters and singletons
def compile_singles_clusters(filename,single_clusters,rmsd_clusters):
	SS_comp = []
	for item in sorted(single_clusters):
		single = [filename.replace(".LP","")]+["single"]+list_details(filename,item,"singletons")
		SS_comp.append(single)
	for n,cl in enumerate(rmsd_clusters):
		k = cl[1:]
		for item in k:
			cluster = [filename.replace(".LP","")]+["cluster_"+str(n)]+list_details(filename,item,"cluster_"+str(n))
			SS_comp.append(cluster)
	return SS_comp
#11-compile functions altogether

def aascad_output(filename):
	start_time = datetime.datetime.now()
	temp_directory = "{}".format((filename).replace('.LP','')) 
	if not os.path.exists(temp_directory): os.makedirs(temp_directory) 	#create directory for each LP file
	shutil.copy(filename,"{}/{}".format(temp_directory,filename))  #copy Lp and SUMAj files
	shutil.copy(filename.replace(".LP",".SUMAj"),"{}/{}".format(temp_directory,filename.replace(".LP",".SUMAj")))
	os.chdir(temp_directory) #enter directory
	indexed_hits = parse_output(filename)  #get list of hits
	summary_file = open("{}_summary.csv".format(filename.replace(".LP","")),"a")  #create output file
	summary_writer = csv.writer(summary_file)
	if len(indexed_hits) != 0:
		print "analysing {}...".format(filename)
		print "generate PDB-formatted patterns..."
		get_coordinates(indexed_hits, filename)  #generate PDB files
		if len([fil for fil in glob.glob("*.pdb")]) > 1:
			print "remove redundant patterns, rmsd=0.0..."
			remred_one()
			remred_two()
			print "superposition and clustering..."
			generate_pseudoatoms()
			superposition()
			print "clustering..."
			#cluster hits with similar patterns - rmsd cut off of 1.0A
			matrix_file = "clusters_matrix.txt"
			rmsd_clusters, single_clusters, pymol_clusters = clustering(matrix_file)
			#merge similar patterns into one group/folder
			pymol(pymol_clusters)
			#delete pseudoatoms
			for pdbx in glob.glob("cluster_*/*.pdb"):
				with open(pdbx,"r") as fi: lines=fi.readlines()
				with open(pdbx,"w") as fi:
					linenew=[li for li in lines if li[12:14] != "PS"]
					for x in linenew: fi.write(x)			
			for pdbx in glob.glob("singletons/*.pdb"):
				with open(pdbx,"r") as fi: lines=fi.readlines()
				with open(pdbx,"w") as fi:
					linenew=[li for li in lines if li[12:14] != "PS"]
					for x in linenew: fi.write(x)
			print "printing singletons and clusters..."
			#compile output for each cluster member and singleton
			All_comp = compile_singles_clusters(filename,single_clusters,rmsd_clusters)
			for i in All_comp: summary_writer.writerow(i)
			print "printing chains and unique chains..."
			for n,cl in enumerate(rmsd_clusters):
				k = cl[1:]
				rmsd_chain_set = []
				for item in k:
					pdb_chain2 = list(set([item.split("_")[0][:4]+k[0] for k in item.split("_")[1].split("-")]))
					for ch in pdb_chain2: 
						if ch in ba_dict.keys():
							rmsd_chain_set.append(ba_dict[ch][0])
						else:
							rmsd_chain_set.append(ch)
				if len(sorted(list(set(rmsd_chain_set))))>1:
					output_per_cluster = open("cluster_{}.txt".format(n),"w")
					for ch2 in sorted(list(set(rmsd_chain_set))):
						print>>output_per_cluster, ch2
			#split superposed pattern in each pdb
			for pdb_ex in glob.glob("cluster_*/cluster*.pdb"):
				lines=open(pdb_ex,"r").read().split("HEADER")[1:]
				for new_pdb_line in lines:
					new_pdbname=new_pdb_line.split("\n")[0].replace(" ","")+"_x.pdb"
					new_pdb_line2="HEADER   "+new_pdb_line
					output = open("/".join(pdb_ex.split("/")[:-1])+"/"+new_pdbname,"w")
					print>>output, new_pdb_line2
		#if there is only one hit the use then compile the output using the following script
		if len([fil for fil in glob.glob("*.pdb")]) == 1:
			ind1 = [filename.replace(".LP","")]+["single"]+list_details(filename,fil,"none")
			summary_writer.writerow(ind1)
			single_directory = "singletons"
			if not os.path.exists(single_directory): os.mkdir(single_directory)
			for pdbx in glob.glob("*.pdb"):
				with open(pdbx,"r") as fi: lines=fi.readlines()
				with open(pdbx,"w") as fi:
					linenew=[li for li in lines if li[12:14] != "PS"]
					for x in linenew: fi.write(x)				
				shutil.move(pdbx,"singletons/"+pdbx)
	end_time = datetime.datetime.now()
	print('Duration: {}'.format(end_time - start_time))

#11-compile functions altogether
def aascad_output(filename):
	start_time = datetime.datetime.now()
	temp_directory = "{}".format((filename).replace('.LP','')) 
	if not os.path.exists(temp_directory): os.makedirs(temp_directory) 	#create directory for each LP file
	shutil.copy(filename,"{}/{}".format(temp_directory,filename))  #copy Lp and SUMAj files
	shutil.copy(filename.replace(".LP",".SUMAj"),"{}/{}".format(temp_directory,filename.replace(".LP",".SUMAj")))
	os.chdir(temp_directory) #enter directory
	indexed_hits = parse_output(filename)  #get list of hits
	summary_file = open("{}_summary.csv".format(filename.replace(".LP","")),"a")  #create output file
	summary_writer = csv.writer(summary_file)
	if len(indexed_hits) != 0:
		print "analysing {}...".format(filename)
		print "generate PDB-formatted patterns..."
		get_coordinates(indexed_hits, filename)  #generate PDB files
		if len([fil for fil in glob.glob("*.pdb")]) > 1:
			print "superposition and clustering..."
			generate_pseudoatoms()
			superposition()
			print "clustering..."
			#cluster hits with similar patterns - rmsd cut off of 1.0A
			matrix_file = "clusters_matrix.txt"
			rmsd_clusters, single_clusters, pymol_clusters = clustering(matrix_file)
			#merge similar patterns into one group/folder
			pymol(pymol_clusters)
			#delete pseudoatoms
			for pdbx in glob.glob("cluster_*/*.pdb"):
				with open(pdbx,"r") as fi: lines=fi.readlines()
				with open(pdbx,"w") as fi:
					linenew=[li for li in lines if li[12:14] != "PS"]
					for x in linenew: fi.write(x)			
			for pdbx in glob.glob("singletons/*.pdb"):
				with open(pdbx,"r") as fi: lines=fi.readlines()
				with open(pdbx,"w") as fi:
					linenew=[li for li in lines if li[12:14] != "PS"]
					for x in linenew: fi.write(x)
			print "printing singletons and clusters..."
			#compile output for each cluster member and singleton
			All_comp = compile_singles_clusters(filename,single_clusters,rmsd_clusters)
			for i in All_comp: summary_writer.writerow(i)
			print "printing chains and unique chains..."
			for n,cl in enumerate(rmsd_clusters):
				k = cl[1:]
				rmsd_chain_set = []
				for item in k:
					pdb_chain2 = list(set([item.split("_")[0][:4]+k[0] for k in item.split("_")[1].split("-")]))
					for ch in pdb_chain2: 
						if ch in ba_dict.keys():
							rmsd_chain_set.append(ba_dict[ch][0])
						else:
							rmsd_chain_set.append(ch)
				if len(sorted(list(set(rmsd_chain_set))))>1:
					output_per_cluster = open("cluster_{}.txt".format(n),"w")
					for ch2 in sorted(list(set(rmsd_chain_set))):
						print>>output_per_cluster, ch2
			#split superposed pattern in each pdb
			for pdb_ex in glob.glob("cluster_*/cluster*.pdb"):
				lines=open(pdb_ex,"r").read().split("HEADER")[1:]
				for new_pdb_line in lines:
					new_pdbname=new_pdb_line.split("\n")[0].replace(" ","")+"_x.pdb"
					new_pdb_line2="HEADER   "+new_pdb_line
					output = open("/".join(pdb_ex.split("/")[:-1])+"/"+new_pdbname,"w")
					print>>output, new_pdb_line2
		#if there is only one hit the use then compile the output using the following script
		if len([fil for fil in glob.glob("*.pdb")]) == 1:
			ind1 = [filename.replace(".LP","")]+["single"]+list_details(filename,fil,"none")
			summary_writer.writerow(ind1)
			os.mkdir("singletons")
			for pdbx in glob.glob("*.pdb"):
				with open(pdbx,"r") as fi: lines=fi.readlines()
				with open(pdbx,"w") as fi:
					linenew=[li for li in lines if li[12:14] != "PS"]
					for x in linenew: fi.write(x)				
				shutil.move(pdbx,"singletons/"+pdbx)
	end_time = datetime.datetime.now()
	print('Duration: {}'.format(end_time - start_time))

#----------------------SECOND OPTIONS: FOLD COMPARISON (LINUX)-------------------------
#01-use dalilite to compare folds
def dalilite_readbrk(clustertext):
	proc2 = subprocess.Popen('DaliLite -AllAll {}'.format(clustertext),shell=True,stdout=subprocess.PIPE) #perform pairwise fold comparison
	out=proc2.communicate()
	if "dali.default" in [fi for fi in glob.glob("*")]: os.remove("dali.default")
	if "dali.lock" in [fi for fi in glob.glob("*")]: os.remove("dali.lock")

#02-parse dccp files to generate cluster and matrix of fold
def parse_dccp(fdir,clustertext):
	os.chdir(fdir)
	print "running DaliLite for fold comparison.."
	dalilite_readbrk(clustertext)
	print "compiled Z-scores.."
	chains = [chain.replace("\n","") for chain in open("{}".format(clustertext),"r").readlines()]
	combinations = ["\t".join(sorted(i)) for i in itertools.combinations(chains,2)]
	if [files for files in glob.glob("*.dccp")] !=[]:
		N = []
		output = open(clustertext.replace(".txt","_fold.txt"),"w") #get list of pairwise fold similarity value in Z-score
		for dccp_file in glob.glob("*.dccp"):
			input = open(dccp_file,"r")
			dccp_1 = [line[69:80].split(" ")+[line[1:5]]+[int(line[5:9])]+[float(line[9:18])]+[float(line[18:22])]+[int(line[22:26])]+[float(line[26:34])]+[int(line[43:46])]+[int(line[46:53])] for line in input.readlines() if line.startswith(" DCCP")]
			dccp_2 = list(set(['\t'.join(i[:2]+[str(i[7])]+[str(i[8])]) for i in dccp_1]))
			for i in dccp_2: N.append(i)
		dc = sorted(sorted(sorted(list(set(N)),key=lambda x: x.split("\t")[3]),key=lambda x: x.split("\t")[2]),key=lambda x: sorted(x.split("\t")[:2]))
		for i in dc: N.append(i)	
		N2 = [n.split("\t") for n in sorted(list(set(N)))]
		dc_dict = collections.defaultdict(list)
		for n in N2: dc_dict["\t".join(n[:2])].append("\t".join(n[2:]))
		dc_dict2 = {k:sorted(sorted(v,key=lambda x:float(x.split("\t")[0]),reverse=True),key=lambda x:float(x.split("\t")[1]),reverse=True)[0] for k,v in dc_dict.iteritems()}
		dc_dict3 = sorted(["\t".join(list(i)) for i in dc_dict2.items()])
		for k in dc_dict3: print>>output, k
	for dssp_file in glob.glob("*.dccp"): os.remove(dssp_file)
	#os.chdir("/Users/nursyatil/NSAG_PART1/ASS_EXE/")
	os.chdir("../../")

def fold_comparison(fdir):
	print "Start fold comparison..."
	fdirall=sorted([i for i in glob.glob("{}*/cluster_*.txt".format(fdir)) if "clusters_pseudoatoms" not in i],reverse=True)
	for fdira in fdirall:
		print fdira
		fdir="/".join(fdira.split("/")[:-1])+"/"
		clustertext=fdira.split("/")[-1]
		parse_dccp(fdir,clustertext)

#----------------------THIRD OPTIONS: SEQUENCE COMPARISON-------------------------
def generate_sequence_identity(matrix_dir):
	for dir in [s for s in glob.glob(matrix_dir+"4_x/*") if ".LP" not in s if ".SUMAj" not in s if "no_hits" not in s]:
		os.chdir(dir)
		lst_cl = [s for s in glob.glob("cluster_*.txt") if "fold" not in s]
		if lst_cl !=[]:
			print dir
			proc = subprocess.Popen("R < /Users/nursyatila/NSAG_PART1/ASS_EXE/sequence_comparison.R --no-save", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
			proc.communicate()
		os.chdir(matrix_dir)
 
def save_matrix(matrixdir):
	outputnf=open("single_fold.txt","a")
	outputns=open("single_seq.txt","a")
	for fdira in glob.glob("{}/*/*/*_fold.txt".format(matrixdir)):
		print fdira
		#fdira = "/Users/nursyatila/NSAG_PART1/ASS_EXE/symmetrical_patterns/single/6/AAAAAA/cluster_0_fold.txt"
		allchains=[i.replace("\n","") for i in open(fdira.replace("_fold.txt",".txt"),"r").readlines()]
		#fold comparison
		if len([s for s in open(fdira,"r").readlines()]) ==1 and "No data" in open(fdira,"r").read():  #NO FOLD
			if len([s for s in open(fdira.replace("_fold","_seqid"),"r").readlines()]) ==1 and "No data" in open(fdira.replace("_fold","_seqid"),"r").read(): #NO FOLD, NO SEQUENCE
				outputmatrix=open(fdira.replace("_fold.txt","_mat.txt"),"w")
				print>>outputmatrix, "No data obtained.\n"
				print>>outputmatrix, "==========Fold Similarity (Z-score)==========\n"
				print>>outputmatrix, "No data obtained.\n"
				print>>outputmatrix, "\nSequences: Non-homologous proteins (Seq.Identity < 30%)"
				print>>outputns, fdira
				print>>outputmatrix, "Fold: Non-homologous proteins (Z-score < 2.0)"
				print>>outputnf, fdira
			else:	#NO FOLD, YES SEQUENCES
				seq_data = [r for k in [[float(n.replace("\n","")) for n in s.split("\t")[1:]] for s in open(fdira.replace("_fold","_seqid"),"r").readlines()[1:]] for r in k]
				seqid=np.array([[float(n.replace("\n","")) for n in s.split("\t")[1:]] for s in open(fdira.replace("_fold","_seqid"),"r").readlines()[1:]])
				dfseq = pd.DataFrame(seqid, index=allchains, columns=allchains)
				outputmatrix=open(fdira.replace("_fold.txt","_mat.txt"),"w")
				print>>outputmatrix, "==========Sequence Identity (%)==========\n"
				print>>outputmatrix, dfseq
				print>>outputmatrix, "\n\n"
				print>>outputmatrix, "==========Fold Similarity (Z-score)==========\n"
				print>>outputmatrix, "No data obtained.\n"
				seqscheck=[k for k in seq_data if float(k)<0.3]
				if seqscheck !=[]: 
					print>>outputmatrix, "\nSequences: Non-homologous proteins (Seq.Identity < 30%)"
					print>>outputns, fdira 
				print>>outputmatrix, "Fold: Non-homologous proteins (Z-score < 2.0)"
				print>>outputnf, fdira			
		else: #YES FOLD
			fold_dict=collections.defaultdict(list)
			lst_fold=[i.replace("\n","").split("\t") for i in open(fdira,"r").readlines()]
			for l in lst_fold: fold_dict[l[0]+"-"+l[1]].append(l[2])
			lst_chain_comb=["-".join(k) for k in itertools.product(allchains,repeat=2)]
			fold_comb=[k.split("-")+fold_dict[k] if k in fold_dict.keys() else k.split("-")+["0.0"] for k in lst_chain_comb]
			datafold=np.array([i[2] for i in fold_comb])
			nfold=len(allchains)
			shapefold=(nfold,nfold)
			dfold = datafold.reshape(shapefold)
			dffold = pd.DataFrame(dfold, index=allchains, columns=allchains)
			#seq comparison
			if len([s for s in open(fdira.replace("_fold","_seqid"),"r").readlines()]) ==1 and "No data" in open(fdira.replace("_fold","_seqid"),"r").read():
				outputmatrix=open(fdira.replace("_fold.txt","_mat.txt"),"w")
				print>>outputmatrix, "==========Sequence Identity (%)==========\n"
				print>>outputmatrix, "No data obtained.\n"
				print>>outputmatrix, "\n\n"
				print>>outputmatrix, "==========Fold Similarity (Z-score)==========\n"
				print>>outputmatrix, dffold
				print>>outputmatrix, "\nSequences: Non-homologous proteins (Seq.Identity < 30%)\n"
				print>>outputns, fdira
				foldcheck=[k for k in [i for i in fold_comb if len(sorted(list(set(i[:-1]))))==2] if float(k[2])<2.0]
				if foldcheck !=[]: 
					print>>outputmatrix, "\nFold: Non-homologous proteins (Z-score < 2.0)\n"
					print>>outputnf, fdira				
			else:
				seq_data = [r for k in [[float(n.replace("\n","")) for n in s.split("\t")[1:]] for s in open(fdira.replace("_fold","_seqid"),"r").readlines()[1:]] for r in k]
				seqid=np.array([[float(n.replace("\n","")) for n in s.split("\t")[1:]] for s in open(fdira.replace("_fold","_seqid"),"r").readlines()[1:]])
				dfseq = pd.DataFrame(seqid, index=allchains, columns=allchains)
				outputmatrix=open(fdira.replace("_fold.txt","_mat.txt"),"w")
				print>>outputmatrix, "==========Sequence Identity (%)==========\n"
				print>>outputmatrix, dfseq
				print>>outputmatrix, "\n\n"
				print>>outputmatrix, "==========Fold Similarity (Z-score)==========\n"
				print>>outputmatrix, dffold
				seqscheck=[k for k in seq_data if float(k)<0.3]
				foldcheck=[k for k in [i for i in fold_comb if len(sorted(list(set(i[:-1]))))==2] if float(k[2])<2.0]
				if seqscheck !=[]: 
					print>>outputmatrix, "\nSequences: Non-homologous proteins (Seq.Identity < 30%)\n"
					print>>outputns, fdira
				if foldcheck !=[]: 
					print>>outputmatrix, "\nFold: Non-homologous proteins (Z-score < 2.0)\n"
					print>>outputnf, fdira
