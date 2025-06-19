# Tutorial-go-bonds

Hi welcome :)

Something like an intro to go bonds and 2beg...
...
...
...



## -- Setup for martinize2 --
Before we start we need to make sure that the new martinize2 script can be called.
So we need to install martinize2 from the github page:
pip install git+https://github.com/niek-fraaije/vermouth-martinize.git@go-implementation
If you want you can check the installation with: which martinize2

Other dependencies:
mdtraj
MDAnalysis (version 2.8.0)
Gromacs version: 2024.4

NOTE: It is advised to work in a virtual environment (pyhton=3.12)


## -- Tutorial --
### Step 1)
First we need to set up our experiment.
So we need to copy our template folder, (e.g. cp -r template/ ./2beg_stacking)

### step 2)
Head over to the rot_trans_fitting directory (cd rot_trans_fitting)
You will find two files here, 2beg.pdb and Stacking.py
2beg.pdb is the first NMR model of the 2BEG entry in the PDB database
Stacking.py is the script that we will use to create a new fiber, based on how the fibrils in 2beg.pdb are stacked ontop of each other.
The script is set to make a fiber composed of 26 fibrils/monomers, the structure of the fiber is saved as multiple_stacked_chains.pdb.

Run the script (python3 Stacking.py)
The Stacking.py script might give some warnings about deprecation of some packages, those can be ignored

NOTE: This script uses the the second and third chain to calculate how the third chain needs to translate and rotate to fit 'best' onto the forth chain ('best' because we now use the rmsd of the alpha carbons as metric for the fitting, but there is no real evidence this is the best metric). With the translation and rotation of one chain onto another, you can stack as many as you want.

### Step 3)
The multiple_stacked_chains.pdb needs to be moved/copied to the init directory (the init directory is located inside the copied template; 2beg_stacking)
Move or copy the file (e.g. mv rot_trans_fitting/multiple_stacked_chains.pdb 2beg_stacking/init/)

### Step 4)
Make sure you work in the copied template directory (2beg_stacking)
Navigate to the itps folder (cd itps) here you will generate all the additional itp files with martinize2.

Run the following command (might take a few minutes):
martinize2 -f ../init/multiple_stacked_chains.pdb -o ../init/topol.top -x ../init/martinized.pdb -ss CEEEEEEEEESCCSEEEEEEEEEEEC -p backbone -ff martini3001 -go -go-ff martini_v3.0.0_go.itp -go-write-file ../init/cont.txt -go-eps-inter 5.0 -go-eps-intra 15.0 -go-up-intra 1.2

The multiple_stacked_chains.pdb is taken as input for martinize2, which generates the course grain structure of the input. The course grain structure (Martini3) is saved as martinized.pdb.
The -go flag indicates that go_bonds should be included, -go-ff is the forcefield file where the non-bonded parameters between course grain beads are defined. 
The -go-write-file is optional, it will write the contact map which is used to create the go-bonds. -go-eps-inter is the epsilon value for all inter molecular go-bonds.
-go-eps-intra is the epsilon value for all intra molecular go-bonds. -go-up-intra is the cutoff, until what distance go-bonds should be included (default is 1.1).
-p backbone makes sure the position restraints are on the backbone beads. And -ss sets the secondary structure of the backbone beads (coil, beta-sheet and bends; C, E and S respectifely)

If everything went well, there should now be additional itp files in you itps folder. 
Namely, go_nbparams.itp, go_atomtypes.itp and molecule_0.itp
The molecule_0.itp describes one of the fibrils (beads, bonds, angles, etc.). There are also many virtual sites added to all the backbone beads, these are the beads that form the go-bonds.
And these virtual sites are not defined in the forcefield, thus are defined in the go_atomtypes.itp. The specific go-bonds are defined in go_nbparams.itp.

### Step 5)
martinize2 does not automatically generate correct topology files. The amount of molecules are not correct and not all itps are referenced or referenced correctly. Thus we need to patch the topol.top.
You should navigate to the copied template directory (2beg_stacking). Here you can run the python script patch_topol.py: python3 patch_topol.py molecule_0 26
molecule_0 is the name of one fibril (defined in molecule_0.itp) and 26 is the amount of fibrils in our system.

### Step 6)
Before we can run the simulation we should convert the output from martinize2 to a gro-file, head over to the init folder (cd init). 
Here you will find martinize.pdb which we can convert to martinize.gro with: gmx editconf -f init/martinized.pdb -o init/martinized.gro -box 20 20 20

The box we simulate is 20 by 20 by 20 nm

### Step 7)
Now you are set to run script.sh in your folder (2beg_stacking)
so navigate to the folder, then run the script (bash script.sh)

NOTE: script.sh will run all gromacs commands for you, from minimization to equilibration and production run. And assumes in the production run that there are 10 freen GPU nodes. The production run will generate 1 microsecond
