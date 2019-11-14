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


#=========PART 2a: SEQUENTIALLY RUN ASSAM SEARCHES===========
def res_run_assam(output_dir,res):
	input_files = [i for i in glob.glob("BINDING_SITES/{}/*.pdb".format(res))]
	for input_file in input_files:
		download_dir = input_file.split("/")[-1].replace(".pdb","")
		download_dir_path = "{}/{}/{}".format(output_dir,res,download_dir)
		#download_dir_path = "drreposer_11_04_19_output/{}/{}".format(res,download_dir)
		if not os.path.exists(download_dir_path): os.makedirs(download_dir_path)
		os.chdir(download_dir_path)
		shutil.move(input_file,"hebe.pdb") #move
		shutil.copy("/Users/nursyatila/NSAG_PART2/ASS_EXE/assam-run.sh","assam-run.sh")
		print "searching for similar patterns as ",download_dir,"..."
		proc = subprocess.Popen('perl /Users/nursyatila/NSAG_PART2/ASS_EXE/res_no.pl {} {}'.format(res,download_dir),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		proc.communicate()
		proc = subprocess.Popen('source assam-run.sh > apa.txt 2>&1',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		proc.communicate()
		list_to_remove = ["AUTO.LPA","assam-run.sh","AUTO.PAT","AUTO.SUMA","AUTO.SUMS","AUTO.XUMS","AUTOASP.INP","fort.22","auto.vek","PATTERN.ID","PREVIOUS.LPA","RESULT_AAMP.txt"]
		for l in list_to_remove: os.remove(l)
		print "done search!"
		os.chdir("/Users/nursyatila/NSAG_PART2/ASS_EXE/")

#=========PART 2b: PARSE OUTPUT FILES===========
#01-remove patterns with no hits
def remove_no_hits():
	if not os.path.exists("no_hits"): os.makedirs("no_hits") #consisting all PDB-formatted patterns
	print "remove filenames with no hits.."
	for filename in glob.glob("*/AUTO.LPS"):
		filedir = filename.split("/")[0]
		patterns01 = re.findall(r"KERB.*coo.*",open(filename,"r").read(),re.S)
		if patterns01 == []:
			os.makedirs("no_hits/{}".format(filedir))
			for fi in glob.glob("{}/*".format(filedir)): shutil.move(fi,"no_hits/{}/{}".format(filedir,fi.split("/")[-1]))
			shutil.rmtree(filedir)
		if patterns01 != []:
			patterns02 = [line for line in patterns01[0].split("\n") if line.startswith("KERBm")]
			patterns03 = [line for line in patterns02 if filedir[:4].lower() not in line]
			if patterns03 == []:
				os.makedirs("no_hits/{}".format(filedir))
				for fi in glob.glob("{}/*".format(filedir)): shutil.move(fi,"no_hits/{}/{}".format(filedir,fi.split("/")[-1]))
				shutil.rmtree(filedir)
	
def pfam_annotation(pdbchains0):
	#pdbchains = "4grvA" or "4grvA,4grvB"
	pdbchains_all_annot = []
	pfam_annot = [[i.split("\t")[0]+i.split("\t")[1]]+[i.split("\t")[4].split(".")[0]]+[i.split("\t")[5]] for i in open("/Users/nursyatila/NSAG_PART2/ASS_EXE/pfampdball.txt","r").readlines()[1:]]
	for pdbchain in pdbchains0.split(","):
		#pfam_annot_matched = ", ".join(sorted(list(set([i[1]+"("+i[2]+")" for i in pfam_annot if i[0]==pdbchain.upper()]))))
		pfam_annot_matched = ";".join(sorted(list(set([i[1]+"("+i[2]+")" for i in pfam_annot if i[0]==pdbchain.upper()]))))
		if pfam_annot_matched == "": pfam_annot_matched = "no annotation"
		pdbchains_all_annot.append(pfam_annot_matched)
	return";".join(pdbchains_all_annot)

#06a-get sequences
def get_sequence(pdbchain):
	sequences = [i for i in re.findall(r">[0-9][a-z0-9][a-z0-9][a-z0-9][a-z0-9A-Z]\n.*",open("/Users/nursyatila/NSAG_PART2/db/pdb_sequence.txt","r").read())]
	dict_sequence_mappings = collections.defaultdict(list)
	for s in sequences: dict_sequence_mappings[s.split("\n")[0][1:]].append(s)
	return "".join(dict_sequence_mappings[pdbchain])

#06b-run sequence similarity comparison using Bio3D package in R
def sequence_alignment():
	proc = subprocess.Popen('R {} --no-save < /Users/nursyatila/NSAG_PART2/ASS_EXE/sequence_comparison.R',shell=True,stdout=subprocess.PIPE)
	out=proc.communicate()

#07-sorted by RMSD values
def save_output(output_dir,filedir):
	#drugbank mappings
	drugbank_mappings = [[i.replace("\n","").split("\t")[0]]+[i.replace("\n","").split("\t")[1]+"("+i.replace("\n","").split("\t")[2]+")"] for i in open("/Users/nursyatila/NSAG_PART2/db/pdb_drugbank.txt","r").readlines()]
	dict_drugbank_mappings = collections.defaultdict(list)
	for d in drugbank_mappings: dict_drugbank_mappings[d[0]].append(d[1])
	#source mappings
	source_mappings = [[i.replace("\n","").split("\t")[0]]+[i.replace("\n","").split("\t")[1]] for i in open("/Users/nursyatila/NSAG_PART2/db/pdb_source.txt","r").readlines()]
	dict_source_mappings = collections.defaultdict(list)
	for s in source_mappings: dict_source_mappings[s[0]].append(s[1])
	#macromolecule mappings
	macromolecule_mappings = [[i.replace("\n","").split("\t")[0]]+[i.replace("\n","").split("\t")[1]] for i in open("/Users/nursyatila/NSAG_PART2/db/pdb_macromolecule.txt","r").readlines()]
	dict_macromolecule_mappings = collections.defaultdict(list)
	for m in macromolecule_mappings: dict_macromolecule_mappings[m[0]].append(m[1])
	os.chdir(filedir)
	#generate all directories
	if not os.path.exists("domains"): os.makedirs("domains") #preparing for Dali analyses
	patterns01 = re.findall(r"KERB.*coo.*",open("AUTO.LPS","r").read(),re.S)
	if patterns01 != []:
		pdb_chains_compilation = [] #to generate single-chain PDB file for Dali analyses
		output = open(filedir+"_sum.csv","a") #compilation of first output without sequence and fold similarity values
		pdb_chains_list_output1 = open("domains/all_chains.txt","w") #list of all chains for Dali analyses
		pdb_chains_list_output2 = open("all_chains.txt","w") #list of all chains for sequence analyses
		csv_writer = csv.writer(output,delimiter=",")
		output_compilations_all = []
		patterns0x = [i.split("\n") for i in patterns01[0].split("\n \nKERBm")]
		patterns02 = [line for line in patterns0x if filedir[:4].lower() not in line]
		print "analyse hits for {}".format(filedir)
		for bil,pat in enumerate(patterns02): #perform analysis per hit
			patt = [line.split("matches")[1].split(">")[0].replace("  <","").replace("    ","") for line in pat if line.startswith("         ")] 
			pdb_id = re.findall(r"[0-9][a-z0-9][a-z0-9][a-z0-9]", pat[0].split("^")[0])[0].lower() 
			rmsd = pat[1].split("=")[1].replace(" ","")
			pdbchains0 = sorted(list(set([pdb_id+n[0] for n in patt])))
			for chain in pdbchains0: pdb_chains_compilation.append(chain)
			organism = ";".join([";".join(dict_source_mappings[pdbchain]) for pdbchain in pdbchains0]).replace(",",";")
			if organism == "": organism = "-"
			macromolecule = ";".join([";".join(dict_macromolecule_mappings[pdbchain]) for pdbchain in pdbchains0]).upper()
			if macromolecule == "": macromolecule = "-"
			pfam_id = pfam_annotation(",".join(pdbchains0))
			matched_residues = ";".join([n[-3:]+" "+n[0]+" "*(4-len(n[1:5].strip(" ")))+n[1:5].strip(" ") for n in patt])
			pattx = sorted(list(set([line.split("matches")[1].split(">")[1].split("from")[1][-12:-3] for line in pat if line.startswith("         ") if line.split("matches")[1].split(">")[1] != " " ]))) #get HETATM codes
			if pattx == []: hetatm_residues = "None"
			else: hetatm_residues = ";".join([p.split("(")[-1].replace(")","")[:-2] for p in pat if "matches" in p])
			if "from" in hetatm_residues:
				if "site" not in hetatm_residues:
					rr = hetatm_residues
					het_resx = [" ".join([r.split(" from ")[1].replace("s","")[-9:]]+["("+r.split(" from ")[0].replace(" ","")+")" if "-" in r.split(" from ")[0] else "( "+r.split(" from ")[0].replace(" ","")+")" ]) if float(re.findall("[0-9]+[.][0-9]+",r.split(" from ")[0])[0])<5.0 else "None" for r in rr.split(";")]
					if all(item == "None" for item in het_resx) == True: het_res = "None"
					else: het_res = ";".join(het_resx)
				else: het_res = "None"
			else: het_res = "None"
			r = pat[2].split()[1:]
			matrix =  ";".join([r[0]]+[r[1]]+[r[2]]+[r[9]]+[r[3]]+[r[4]]+[r[5]]+[r[10]]+[r[6]]+[r[7]]+[r[8]]+[r[11]]+[str(0.0)]+[str(0.0)]+[str(0.0)]+[str(1.0)])
			csv_writer.writerow([filedir]+[pdb_id]+[organism]+[macromolecule]+[pfam_id]+[matched_residues]+[het_res]+[rmsd]+[matrix])
		print "generating sequences for sequence comparison.."
		output2 = open("hits.fasta","a") #get sequences from hits
		for s in sorted(list(set(pdb_chains_compilation))): 
			print>>pdb_chains_list_output1, s
			print>>pdb_chains_list_output2, s		
			print>>output2, get_sequence(s)
		query_chains = sorted(list(set([filedir[:4].lower()+line[21:22] for line in open("hebe.pdb","r").readlines() if line[:4] == "ATOM"])))
		output1 = open("query.fasta","w") #get sequences from query pattern (a/w)
		output1a = open("domains/query_chains.txt","w")
		for i in query_chains:
			print>>output1, get_sequence(i)
			print>>output1a, i	
	print "done saving output!"

#07-get unique items
def unique_items(L):
    found = set()
    for item in L:
        if '\t'.join(item[0:2]) not in found:
            yield item
            found.add('\t'.join(item[0:2]))

#08--convert pdb to DAT for DaliLite
def dalilite_readbrk(filedir,query_chains):
	os.chdir("domains")
	for x in query_chains:
		proc2 = subprocess.Popen('DaliLite -list {} all_chains.txt'.format(x),shell=True,stdout=subprocess.PIPE) #perform pairwise fold comparison
		out=proc2.communicate()
		if "dali.default" in [fi for fi in glob.glob("*")]: os.remove("dali.default")
		if "dali.lock" in [fi for fi in glob.glob("*")]: os.remove("dali.lock")

#09-parse dccp files to generate cluster and matrix of fold
def parse_dccp(filedir,query_chains):
	chains = [chain.replace("\n","") for chain in open("all_chains.txt","r").readlines()]
	combinations = [s for n in [[filename_chain+"\t"+chain for chain in chains] for filename_chain in query_chains] for s in n]
	if [files for files in glob.glob("*.dccp")] !=[]:
		N = []
		output = open("fold_comparison.txt","w") #get list of pairwise fold similarity value in Z-score
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
	for dssp_file in glob.glob("*.dssp"): os.remove(dssp_file)
	for pdb_file in glob.glob("*.pdb"): os.remove(pdb_file)

#10-combine all command for fold comparison
def compile_all(output_dir):
	os.chdir(output_dir)
	for filedirs in glob.glob("*/AUTO.LPS"):
		filedir = filedirs.split("/")[0]
		os.chdir(filedir)
		query_chains = sorted(list(set([filedir[:4]+line[21:22] for line in open("hebe.pdb","r").readlines() if line[:4] == "ATOM"])))
		print filedir
		print "fold comparison using DaliLite"
		dalilite_readbrk(filedir,query_chains)
		print "get pairwise Z-score and sequence identity"
		parse_dccp(filedir,query_chains)
		print "done!"
		os.chdir(output_dir)

#11-add sequence and fold information in csv files
def insert_fold_and_sequence(filedir,output_dir):
	os.chdir(filedir)
	if "domains/fold_comparison.txt" in [fi for fi in glob.glob("domains/*")]:
		query_chains = sorted(list(set([filedir[:4].lower()+line[21:22] for line in open("hebe.pdb","r").readlines() if line[:4] == "ATOM"]))) #get query chains ["1dctA","1dctB"]
		list_chain = sorted(list(set([x+"-"+i.replace("\n","") for i in open("domains/all_chains.txt","r").readlines() for x in query_chains]))) #get hit chains
		sum1 = [row for row in csv.reader(open(filedir+"_sum1.csv","r"),delimiter=",")] #get original csv rows
		sum = [row+[sorted(list(set([x+"-"+row[1]+k[4:5] for k in row[5].split("; ") for x in query_chains])))] for row in sum1] #edit rows containing hits #1dctA-
		#fold = [[i.replace("\n","").split("\t")[0]+" "+i.replace("\n","").split("\t")[1]]+[i.replace("\n","").split("\t")[2]] for i in open("domains/fold_comparison.txt","r").readlines()] #get pairwise z-score
		fold = [[i.replace("\n","").split("\t")[0][:4].lower()+i.replace("\n","").split("\t")[0][4:5]+" "+i.replace("\n","").split("\t")[1]]+[i.replace("\n","").split("\t")[2]] for i in open("domains/fold_comparison.txt","r").readlines()] #get pairwise z-score
		seq = [[j.replace("\n","").split(" ")[0]+" "+j.replace("\n","").split(" ")[1]]+[str("{0:.2f}".format(float(j.replace("\n","").split(" ")[2])))] for j in open("sequence_identity.txt","r").readlines()] #get pairwise seq id
		fold_dict = collections.defaultdict(list) #get list of fold similarity value
		seq_dict = collections.defaultdict(list) #get list of sequence identity values
		for f in fold: fold_dict[f[0].replace(" ","-")].append(f[1])
		for s in seq: seq_dict[s[0].replace(" ","-")].append(s[1])
		sums2 = [row[:8]+[";".join([k+":"+"".join(fold_dict[k]) if k in fold_dict.keys() else k+":undetectable" for k in row[9]])]+[";".join([k+":"+"".join(seq_dict[k]) for k in row[9]])]+[row[8]] for row in sum]
		csv_writer = csv.writer(open(filedir+"_sum2.csv","w"),delimiter=",")
		for k in sums2: csv_writer.writerow(k) #compile second output, now containing fold and sequence similarity values
		print filedir, "done compiling fold and sequences!"
	os.chdir(output_dir)

#12-command line arguments
if __name__ == '__main__': # 13- PARSE COMMAND LINE ARGUMENTS
	usage = "USAGE: python AUTO_ASSAM_PART2.py -o output_dir -s (search, hits, fold or merge)"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option("-s", "--selection_o", dest="selection_o", type=str, help="hits or fold or merge")
	parser.add_option("-o", "--output_dir", dest="output_dir", type=str, help="output_directory")
	(options, args) = parser.parse_args()

	if (options.output_dir and options.selection_o == "search"):
	for r in range(3,13):
		print "running assam for res {}".format(res)
		res_run_assam(options.output_dir,res)
		
	if (options.output_dir and options.selection_o == "hits"):
		os.chdir(options.output_dir)
		remove_no_hits()
		for filedirs in glob.glob("*/AUTO.LPS"):
			filedir = filedirs.split("/")[0]
			save_output(options.output_dir, filedir)
			if "hits.fasta" in glob.glob("*.fasta"):
				sequence_alignment()
			os.chdir(options.output_dir)
		os.chdir(options.output_dir)
		print "done sort hits!!!"
	
	if (options.output_dir and options.selection_o == "fold"):
		compile_all(options.output_dir)
		print "done comparing fold and sequence!!!"
		
	if (options.output_dir and options.selection_o == "merge"):
		os.chdir(options.output_dir)
		for filedirs in glob.glob("*/AUTO.LPS"):
			filedir = filedirs.split("/")[0]
			insert_fold_and_sequence(filedir,options.output_dir)
		os.chdir(options.output_dir)
		print "done inserting sequence and fold similarity!!!"
	
	if not (options.output_dir and options.selection_o):
		print "Not enough input arguments supplied"
		print usage
		quit()