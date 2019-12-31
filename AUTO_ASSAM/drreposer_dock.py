#!/usr/bin/env python
import sys
import os
import chimera
from chimera import runCommand
from AddCharge import initiateAddCharges
initiateAddCharges(method='am1-bcc', nogui=True)
#use: chimera --nogui --script drreposer_run_vina.py OR open with UCSF Chimera
#Have to install Chimera first

#prepare ligand molecule for docking - ligand derived from complexes
#e.g. 1sn0 - protein-drug complex
#e.g. 602.B - drug molecule
runCommand('open 1sn0; sel :602.B; sel invert; del sel; del solvent; del ions; wait')
runCommand('swapaa same sel preserve true ignoreOtherModels true; addh; addcharge s; wait')
#save files to current directory
runCommand('write format mol2 0 predictedbs.ligand.mol2; wait')

#prepare receptor structure for docking - give a pdb id and only select the first model in the case of nmr structure
#e.g 2k0z - protein to be search for its potential of binding to ligand
runCommand('open 2k0z; del #1.2-100; sel protein; sel invert sel; del sel; del solvent; del ions; wait')
runCommand('swapaa same sel preserve true ignoreOtherModels true; addh; addcharge s; wait')
#save files to current directory
runCommand('write format mol2 1 predictedbs.receptor.mol2; wait')

#display predicted residues similar to known drug binding site - colored in orange
runCommand('open 2k0z; del #2.2-100; sel #2:6.A,74.A,70.A,73.A; color orange sel; display sel; ~ribbon sel; represent b+s sel; sel invert sel; del sel; wait')

#run AutoDock Vina and retrieve potential protein-drug conformations
cwd = os.getcwd()
runCommand('open predictedbs.ligand.mol2; open predictedbs.receptor.mol2; wait')
runCommand('vina docking receptor #4 ligand #3 output {}/predictedbs'.format(cwd))



