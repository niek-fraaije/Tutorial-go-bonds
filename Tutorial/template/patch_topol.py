import sys

if len(sys.argv) < 3:
	print('no mol_name and/or amount given to patch_topol.py')
	sys.exit(1)
	
mol_name, amount = sys.argv[1], sys.argv[2]

topol = []
with open('init/topol.top', 'r') as file:
	for line in file:
		if line.startswith('#include "martini.itp"'):
			topol.append('#include "../itps/martini_v3.0.0_go.itp"\n')
			topol.append('#include "../itps/martini_v3.0.0_ions_v1.itp"\n')
			topol.append('#include "../itps/martini_v3.0.0_solvents_v1.itp"\n')
			topol.append('#include "../itps/martini_v3.0.0_phospholipids_v1.itp"\n')
		elif line.startswith('#include'):
			parts = line.split()
			include_file = parts[1].strip('"')
			new_path = f'"../itps/{include_file}"'
			line = f'#include {new_path}\n'
			topol.append(line)

		elif line.startswith(mol_name):
			line = line.split()
			line[1] = amount
			line = ' '.join(line)
			topol.append(line)
		else:
			topol.append(line)
	topol.append('\n')
		
		
with open('init/topol.top', 'w') as file:
	file.writelines(topol)
	
