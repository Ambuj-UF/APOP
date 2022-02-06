import numpy
import os
import shutil
from src import *

def apop_cli(filename):
    
    output_folder = "APOP_" + filename.split(".")[0]
    try:
        os.system("mkdir %s" %output_folder)
    except:
        os.system("rm -r %s" %output_folder)
        os.system("mkdir %s" %output_folder)
    
    Model = Allostery(filename, chain = "All", cutoff=10.0, active_site=[])
    store_pockets = Model.get_pockets()
    with open(output_folder + "/apop_output.txt", "w") as fp:
        for key, val in store_pockets.items():
            fp.write("%s\n" %key)
            fp.write("Pocket name: %s\n" %val[0])
            fp.write("APOP score: %s\n" %val[1])
            fp.write("Residues: %s\n" %", ".join(val[2]))
            fp.write("---------------------------------------------------\n\n")
            
            src = filename.split(".")[0] + "_out/pockets/" + val[0]
            shutil.copy2(src, output_folder)
            
    shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(filename.split(".")[0] + "_out")
    shutil.rmtree(output_folder)
            
#apop_cli(filename="2V7A.pdb")
