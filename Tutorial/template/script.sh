#!/bin/bash

set -e

# Set working directory
work_dir=$(pwd)

# Minimize in vacuum
gmx grompp -p init/topol.top -f mdps/minimization.mdp -c init/martinized.gro -o em_vac/minimization-vac.tpr -r init/martinized.gro -maxwarn 1

# Navigate to the em_vac directory, to run minimization
cd em_vac || { echo "Directory em_vac not found!"; exit 1; }
gmx mdrun -deffnm minimization-vac -v
cd $work_dir

# Solvate the system
gmx solvate -cp em_vac/minimization-vac.gro -cs init/water.gro -radius 0.21 -o em/solvated.gro -p init/topol.top

# Add ions
gmx grompp -f mdps/minimization.mdp -p init/topol.top -c em/solvated.gro -r em/solvated.gro -o em/ions.tpr
echo W | gmx genion -s em/ions.tpr -p init/topol.top -neutral -conc 0.15 -o em/ions.gro

# Minimize with water and ions
gmx grompp -p init/topol.top -c em/ions.gro -f mdps/minimization.mdp -o em/minimization.tpr -r em/ions.gro -maxwarn 1

# Navigate to the em directory, to run minimization
cd em || { echo "Directory em not found!"; exit 1; }
gmx mdrun -deffnm minimization -v
cd $work_dir

gmx grompp -p init/topol.top -c em/minimization.gro -f mdps/equilibration.mdp -o eq/equilibration.tpr -r em/solvated.gro -maxwarn 3

# Navigate to the eq directory, to run equilibration
cd eq || { echo "Directory eq not found!"; exit 1; }
gmx mdrun -deffnm equilibration -v
cd $work_dir

gmx grompp -p init/topol.top -c eq/equilibration.gro -f mdps/dynamic.mdp -o md/dynamic.tpr

# Navigate to the eq directory, to run dynamics
cd md || { echo "Directory md not found!"; exit 1; }
gmx convert-tpr -s dynamic.tpr -until 1000000 -o dynamic1us.tpr
gmx mdrun -deffnm dynamic -s dynamic1us.tpr -cpi dynamic.cpt -nt 10 -pin on -v

# Generate xtc and gro files with only Protein
echo 1 | gmx trjconv -s dynamic1us.tpr -f dynamic.gro -pbc mol -o no_water.gro
echo 1 | gmx trjconv -s dynamic1us.tpr -f dynamic.xtc -pbc mol -o no_water.xtc

cd $work_dir







