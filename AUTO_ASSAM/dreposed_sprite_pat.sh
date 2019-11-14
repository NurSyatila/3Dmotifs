#!/bin/sh

#  assam_new.sh
#  
#
#  Created by Nur Syatila Ab Ghani on 2/19/16.
#
source ass_assam.sh

for file in PDB_DBS/*.pdb; do
bname=$(basename "$file")
extension="${bname##*.}"
filenamewoext="${bname%.*}"
newfilename="${filenamewoext}.pat"

$ASS_EXE/aamc_pat  <<EOF
$file
${filenamewoext}
PIK_DBS/$newfilename
EOF
rm "${filenamewoext}.vek"
done

#PETE_HIT.NOTE: aspreg (creates a .LPA) > aspsup (L/R superposition, creates a .LPS) > LPS2ras (examine with rasmol, source R.ras)