#!/usr/bin/env python
# -*- coding: utf-8 -*-
#A program to conduct automatic search and clustering of amino acid patterns
import itertools, sys, os, subprocess, shutil, glob, numpy as np, re, collections, operator, datetime, optparse, csv, pandas as pd
from scipy.cluster import hierarchy as hc
from scipy.spatial import distance as dc

trans = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','EOF':''}
strand = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
opt_dict = {'all':'ACDEFGHIKLMNPQRSTVWY','charged-amide':'DEHKRQN','charged-aromatic-amide':'FWYDEHKRQN','charged':'DEHKR','basic':'HKR','acidic':'DE','amide':'QN','aromatic':'FWY','hydrophobic':'LVIAPMFWY','hydroxyl':'CST'}		

#----------------------FIRST OPTIONS: SEARCH FOR HYPOTHETICAL PATTERNS-------------------------
#options for search: 1)R: residue number, 2)N: list of residues for search e.g. DEHRK or -P:exact pattern for search 
#01-get patterns
def input_files(r,type,n_dict): # 02-CREATE INPUT FILES FOR IMAAAGINE SEARCHES
	d = 3.5
	#d = 3.0 #3.5+1.5 distance tolerance -> 2-5 -> 3.5 = 1.5-4+1.5 = 5.5
	x = "\n%.2f" % d #required distance
	y = "\n " #leave blank for wildcard distance, or Y for close contact
	if r == 3: distance = "{}".format(x*3)
	if r == 4: distance = "{}{}{}".format(x*2,y*2,x*2)
	if r == 5: distance = "{}{}{}{}{}{}{}".format(x,y*2,x*2,y*2,x,y,x)
	if r == 6: distance = "{}{}{}{}{}{}{}{}{}".format(x,y*3,x*2,y*3,x,y*2,x,y,x)
	if r == 7: distance = "{}{}{}{}{}{}{}{}{}{}{}".format(x,y*4,x*2,y*4,x,y*3,x,y*2,x,y,x)
	if r == 8: distance = "{}{}{}{}{}{}{}{}{}{}{}{}{}".format(x,y*5,x*2,y*5,x,y*4,x,y*3,x,y*2,x,y,x)
	if type=="list":
		combinations=["\n".join([strand[x] for x in i]) for i in itertools.combinations_with_replacement(n_dict,r)]
	if type=="exact":
		combinations = ["\n".join(strand[x] for x in n_dict)]
	for comb in combinations: inputs = [''.join(''.join(str(r))+'\n'+''.join(comb)+''.join(distance)+'\n'+">>EOF") for comb in combinations]
	return inputs
#02-run imaaagine
def imaaagine(r,n_dict,inputs,output_dir): # 03-RUN IMAAAGINE SEQUENTIALLY
	for input in inputs:
		proc = subprocess.Popen('source ass_imaaagine2.sh \n source hyp_assam.sh',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
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
				if not os.path.exists(output_dir): os.makedirs(output_dir)
				outfile1 = '{}/{}.LP'.format(output_dir,dict)
				outfile2 = '{}/{}.SUMAj'.format(output_dir,dict)
				shutil.move(infile1, outfile1)
				shutil.move(infile2, outfile2) # 04-SAVE THE OUTPUT FILES CONTAIN PDBID AND MATCHED PATTERN
#03-remove hits
def remove_no_hits(r,n_dict,output_dir): #06-REMOVE PATTERNS W/O HIT
	os.chdir(output_dir)
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

def auto_imaaagine(r,type,n_dict):
	print "Generating hypthetical patterns!"
	inputs = input_files(input_files(r,type,n_dict))
	if type=="list": output_directory = "{}/{}".format(n_dict,r)
	if type=="exact": output_directory = "{}".format(n_dict)
	print "Search for all patterns of amino acid in PDB match the query protein!"
	imaaagine(r,n_dict,inputs,output_dir)
	remove_no_hits(r,n_dict,output_dir)
	print "Done search!"

#----------------------SECOND OPTIONS: PARSE OUTPUT FILES-------------------------

#01-parse Lp file and get patterns from hits
def parse_output(filename):
	input = open(filename,"r")
	all = re.findall(r"^(   [PATTERN].*\s     [0-9].*\s.*\s.*\s.*\s.*\s.*\s.*\s.*\szap)",input.read(), re.MULTILINE)
	patterns = collections.OrderedDict.fromkeys([i.replace(i[:54],i[34:38]).replace(i[-4:],"") for i in all]).keys()
	sorted_hits = [[i.split("\n")[0]]+sorted([w[10:13]+w[19:25] for w in i.split("\n")[1:]],key=lambda x:x[4:])+sorted(list(set([w[62:71] for w in i.split("\n")[1:]])),key=lambda x:x[4:]) for i in patterns]
	non_redundant_hits = sorted(list(set([tuple(i) for i in sorted_hits]))) 
	indexed_hits = [i.split("\t") for i in sorted(list(set(["\t".join(x) for x in [[item[0][:4]] + list(i for i in item[1:]) for item in non_redundant_hits]])))]
	indexed_hits2=[i for i in indexed_hits if len(set([k[4:5] for k in i[1:len(filename.replace(".LP",""))+1]]))>1]
	return indexed_hits2

#02-generate PDB-formatted patterns and generate PDB-formatted patterns with pseudoatoms - will later be used for superposition
def get_coordinates(indexed_hits, filename): # 02- GENERATE TEMPORARY PDB FILE FOR EACH MOTIFS
	for m in indexed_hits:
		pdb = m[0][:4]
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
		with open("/Users/nursyatila/NSAG_PART1/ba_all/{}.pdb".format(pdb)) as pdb_file:
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

#03-include pseudoatoms in patterns
def generate_pseudoatoms():
	proc = subprocess.Popen("pymol -qrc /Users/nursyatila/NSAG_PART1/ASS_EXE/imaaagine_pymol_pseudoatoms.py -- 1st", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	proc.communicate()

#04-superposition using pseudoatoms and cluster similar patterns of RMSD value lesser than 1.5A
def superposition_and_clustering():
	proc = subprocess.Popen("R < /Users/nursyatila/NSAG_PART1/ASS_EXE/imaaagine_superposition_pseudoatoms.R --no-save", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	proc.communicate()
		
#05-cluster similar patterns into groups
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
			proc = subprocess.Popen("pymol -qrc /Users/nursyatila/NSAG_PART1/ASS_EXE/imaaagine_pymol_pseudoatoms.py -- 2nd {}".format(i), shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
			proc
			output_pymol = open("{}/rmsd_clusters.txt".format(i),"w")
			while True:
				out = proc.stdout.read()
				if out == '' and proc.poll() !=None: break
				if out != '': output_pymol.write(out)

#07-generate list of details for hits
def list_details(filename,item,group): #07- EXTRACT PFAM ANNOTATION FOR EACH PDB CHAIN
#macromolecule name mappings
	items = item.split("_")
	pdbchains0 = ",".join(sorted(list(set([items[0][:4].lower()+i[0] for i in items[1].split("-")]))))
	csv_reader=csv.reader(open("/Users/nursyatila/NSAG_PART1/ASS_EXE/pdbdescription_ba_all.csv","r"),delimiter=",")
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
	pfam=";".join(sorted(list(set([k for n in [dict_pfam_mappings[pdbchain] for pdbchain in pdbchains0.split(",")] for k in n]))))
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

#10-get unique pdb chains for later fold or sequence comparison
def pdb_chain_unique(filename):
	XY = []
	for item in glob.glob("*.pdb"):
		kk = list(set([item.split("_")[0][:4]+k[0] for k in item.split("_")[1].split("-")]))
		for i in kk: XY.append(i)
	output = open("{}_chain_unique.txt".format(filename.replace(".LP","")),"w")
	for k in sorted(list(set(XY))): print>>output, k.replace("\n","")

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
			#generate_representative_patterns_R()  #get representative patterns
			#if len([fil for fil in glob.glob("*.pdb")]) > 1:
			print "include pseudoatoms..."
			generate_pseudoatoms()  #generate pseudoatoms to be used for superposition
			print "superposition and clustering..."
			superposition_and_clustering()
			print "clustering..."
			#cluster hits with similar patterns - rmsd cut off of 1.5A/0.8A
			matrix_file = "clusters_pseudoatoms.txt"
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
					for ch in pdb_chain2: rmsd_chain_set.append(ch)
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


#13-command line arguments 
if __name__ == '__main__': # 13- PARSE COMMAND LINE ARGUMENTS
	p1="\n USAGE: \n\n python auto_imaaagine.py -q search -r residuenumber -n residuelist e.g. python auto_imaaagine.py -q search -r 4 -n DERHK \n                                  OR \n "
	p2="python auto_imaaagine.py -q search -r residuenumber -p exactpattern e.g. python auto_imaaagine.py -q search -r 3 -p SHD \n                                  OR \n "
	p3="python auto_imaaagine.py -q analyse -o output_dir -f filename e.g. python auto_imaaagine.py -q analyse -o /Users/nursyatila/NSAG_PART1/ASS_EXE/imaaagine_output/ -f DEDE.LP \n                                  OR \n "
	p4="python auto_imaaagine.py -q analyse -o output_dir -f directory e.g. python auto_imaaagine.py -q analyse -o /Users/nursyatila/NSAG_PART1/ASS_EXE/imaaagine_output/ -f 4\n"
	usage=p1+p2+p3+p4
	#usage = "\n USAGE: \n\n python auto_imaaagine.py -q search -r residuenumber -n residuelist \n OR \n python auto_imaaagine.py -q search -r residuenumber -p exactpattern \n OR \n python auto_imaaagine.py -q analyse -o output_dir -f filename \n OR \n python auto_imaaagine.py -q analyse -o output_dir -f directory \n"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option("-q", "--query", dest="query", type=str, help="query")
	parser.add_option("-r", "--resnum", dest="resnum", type=str, help="residue number")
	parser.add_option("-n", "--residuelist", dest="residuelist", type=str, help="list of residues for searching")
	parser.add_option("-p", "--exactpattern", dest="exactpattern", type=str, help="exact pattern for searching")
	parser.add_option("-o", "--output_dir", dest="output_dir", type=str, help="output_directory")
	parser.add_option("-f", "--filename", dest="filename", type=str, help="user.LP")
	parser.add_option("-d", "--directory", dest="directory", type=str, help="directory")
	(options, args) = parser.parse_args()

	if (options.query):
		if options.query=="search":
			if (options.query and options.resnum and options.residuelist):
				r=int(options.resnum)
				n_dict=options.residuelist
				type=="list"
				auto_imaaagine(r,type,n_dict)
			if (options.query and options.resnum and options.exactpattern):
				r=int(options.resnum)
				n_dict=options.exactpattern
				type=="exact"
				auto_imaaagine(r,type,n_dict)
		if options.query=="analyse":
			if (options.output_dir and options.filename):
				os.chdir("{}{}".format(options.output_dir,len((options.filename).replace(".LP",""))))
				aascad_output(options.filename)
				print "done!"
			if (options.output_dir and options.directory):
				os.chdir("{}{}".format(options.output_dir,options.directory))
				for filename in glob.glob("*.LP"):
					aascad_output(filename)
					os.chdir("../")
				print "done!"			
			if not (options.output_dir and options.filename) or (options.output_dir and options.directory):
				print "Not enough input arguments supplied"
				print usage
				quit()
	else:
		print "Not enough input arguments supplied"
		print usage
		quit()
	