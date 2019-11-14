#!/usr/bin/env Rscript
#use: calculate rmsd and cluster similar conformations for all pdbs in pdb_filess
#use: "R < /Users/nursyatila/assam/ASS_EXE/imaaagine_superposition_pseudoatoms.R --no-save"

options(warn = -1)
library("bio3d")
myfiles <- list.files(".", ".pdb$", full.names=TRUE)
xyz <- NULL
for(i in 1:length(myfiles)) {
	pdbs <- read.pdb(myfiles[i])
	all.inds <- atom.select(pdbs, elety=c("PS1","PS2","PS3"))
    xyz <- rbind(xyz, pdbs$xyz[,all.inds$xyz])
}
rd <- rmsd(xyz,fit=TRUE)
hc <- hclust(as.dist(rd), method="average")
groups <- cutree(hc, h=1.0)
dd <- data.frame(groups)
rownames(rd) <- basename(myfiles)
rownames(dd) <- basename(myfiles)
write.table(rd, file="rmsd_matrix.txt", quote=FALSE, col.names=FALSE, sep="\t")
write.table(dd, file="clusters_matrix.txt", quote=FALSE, col.names=FALSE, sep="\t")