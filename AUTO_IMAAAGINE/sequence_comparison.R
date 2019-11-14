#!/usr/bin/env Rscript
#use: calculate rmsd and cluster similar conformations for all pdbs in pdb_filess
#use: "R < /Users/nursyatila/assam/ASS_EXE/imaaagine_superposition_pseudoatoms.R --no-save"

options(warn = -1)
library("bio3d")
library("stringr")
myfiles <- list.files(pattern="^cluster_[0-9]+.txt$")
for(i in 1:length(myfiles)) {
	pdbids <- readLines(myfiles[i])
	pdbidsnew <- paste0(substring(pdbids, 1, 4), "_", substring(pdbids, 5, 6))
	seqs <- get.seq(pdbidsnew)
	aln <- seqaln(seqs)
	matrix <- seqidentity(aln)
	ide.mat <- seqidentity(aln)
	rownames(ide.mat) <- paste0(substring(rownames(ide.mat), 5, 8), "_", substring(rownames(ide.mat), 10,11))
	colnames(ide.mat) <- paste0(substring(colnames(ide.mat), 5, 8), "_", substring(colnames(ide.mat), 10,11))
	write.table(ide.mat, file=str_replace(myfiles[i],".txt","_seqid.txt"), quote=FALSE, col.names=TRUE, sep="\t")
}

#library("bio3d")
#library("stringr")
#myfile="cluster_3.txt"
#pdbids <- readLines(myfile)
#pdbidsnew <- paste0(substring(pdbids, 1, 4), "_", substring(pdbids, 5, 6))
#seqs <- get.seq(pdbidsnew)
#aln <- seqaln(seqs)
#matrix <- seqidentity(aln)
#ide.mat <- seqidentity(aln)
#rownames(ide.mat) <- paste0(substring(rownames(ide.mat), 5, 8), "_", substring(rownames(ide.mat), 10,11))
#colnames(ide.mat) <- paste0(substring(colnames(ide.mat), 5, 8), "_", substring(colnames(ide.mat), 10,11))
#write.table(ide.mat, file=str_replace(myfile,".txt","_seqid.txt"), quote=FALSE, col.names=TRUE, sep="\t")
