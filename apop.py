import os
import sys
import shutil
from src import *
import argparse
import urllib

def apop_cli(filename, pdbid, chain, cutoff, output_filepath=None):

    try:
        mol = molecule.load_structure(filename)
    except:
        try:
            molecule.download_structure(pdbid, ftype="pdb")
        except urllib.error.HTTPError:
            sys.exit("Failed to download structure for pdb id: {}".format(pdbid))

    filename = pdbid + ".pdb"

    if chain != "All":
        mol = molecule.load_structure(filename)
        res = mol[0][chain].get_residues()
        with open(pdbid+".pdb", "w") as fp:
            for residue in res:
                for i in residue.get_atoms():
                    fp.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(i.get_id(),i.get_name(),i.get_parent().get_name(),i.get_parent().get_parent().get_id(),i.get_parent().get_id(),around(i.get_location()[0],decimals=3),around(i.get_location()[1],decimals=3),around(i.get_location()[2],decimals=3),i.get_occupancy(),i.get_bfactor(),'',i.get_element(),''))

    output_folder = "APOP_" + filename.split(".")[0]
    try:
        os.system("mkdir %s" %output_folder)
    except:
        os.system("rm -r %s" %output_folder)
        os.system("mkdir %s" %output_folder)

    shutil.copy2(filename, output_folder+"/"+filename)

    Model = Allostery(filename, pdbid=pdbid, chain = chain, cutoff=cutoff, active_site=[])
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

    if output_filepath is not None:
        zip_filepath = output_filepath.split('.')[0]
    else:
        zip_filepath = output_folder

    shutil.make_archive(zip_filepath, 'zip', output_folder)
    shutil.rmtree(filename.split(".")[0] + "_out")
    shutil.rmtree(output_folder)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,
                        help='input filename')
    parser.add_argument('--terminallog', type=str,
                        help='redirect stdout here')
    parser.add_argument('--chain', type=str, default='All',
                        help='chain')
    parser.add_argument('--pdbid', type=str,
                        help='pdb id')
    parser.add_argument('--cutoff', type=float, default=10.0,
                        help='cutoff')
    parser.add_argument('--output', type=str, help='Path and filename write output to')
    args = parser.parse_args()

    if args.terminallog:
        sys.stdout = open(args.terminallog, 'w')

    if "/" in args.input or "\\" in args.input:
        shutil.copy2(args.input, os.getcwd())
        input_file = os.path.basename(args.input)
    else:
        input_file = args.input

    apop_cli(input_file, args.pdbid, args.chain, args.cutoff, args.output)

    if args.terminallog:
        sys.stdout.close()

if __name__ == "__main__":
    main()
#apop_cli(filename="2V7A.pdb")
