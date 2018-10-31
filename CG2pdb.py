import os

parts_dir = "./inputs/parts/"	# Path to directory containing parts
out_dir = "./output_atomistic"	# name of directory to create and hold atomistic files
w_dir = "./output/pdbs/"	        # name of directory where CG files are

prefix = 'model.'		# File prefix, e.g. "model.0" = "model."
suffix = ''		# File suffix excluding ".pdb"

n_first = 0 	# First file
n_last = 99    # Last file

parts = os.listdir(parts_dir)

for part in parts:
	print "loaded "+part
	cmd.load(parts_dir+part,part)

for template in range(n_first,n_last+1):
	cmd.load(w_dir+prefix+str(template)+".pdb",str(template))
	
	cmd.align("TS_1WBJ_A.pdb","chain A")
	cmd.align("TS_1WBJ_B.pdb","chain B")
	cmd.align("TS_1WBJ_C.pdb","chain C")
	cmd.align("TS_1WBJ_D.pdb","chain D")
	
	cmd.delete(str(template))
	
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	cmd.save("./"+out_dir+"/"+str(template).zfill(2)+"_atomistic.pdb")
	
