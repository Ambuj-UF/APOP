
from packman import molecule
from .gnm import GNM

from functools import reduce
import numpy
from numpy import around
import os
import glob
import sys
import shutil
import urllib



class Allostery():
    def __init__(self, filename, pdbid=None, chain="All", cutoff=10.0, active_site=None):
        self.cutoff = cutoff
        self.filename = filename
        self.chain = chain
        self.active_site = active_site
        self.pdbid = pdbid


    def pull_pockets(self):
        try:
            os.system("fpocket -f" + self.filename)
        except:
            sys.exit("fpocket installation required to run the program\n")

        files = glob.glob(self.filename[:-4] + "_out/pockets/*pdb")
        self.pocket_data = dict()
        for pocket_file in files:
            mol = molecule.load_structure(pocket_file)
            residues = [str(i.get_id())+str(i.get_parent().get_id()) for i in mol[0].get_residues()]
            self.pocket_data[pocket_file] = list(set(residues))

        return True



    def prep_mol(self):
        #print(self.chain)
        self.mol    = molecule.load_structure(self.filename)
        #self.calpha = [i.get_calpha() for i in self.mol[0][self.chain].get_residues()]
        self.allchains = list()
        self.tip = list()
        if self.chain == "All":
            self.calpha = list()
            for i in self.mol[0].get_residues():
                if i.get_calpha() != None:
                    self.calpha.append(i.get_calpha())
                    self.tip.append(i.get_tip())
                    self.allchains.append(i.get_calpha().get_parent().get_parent().get_id())

            self.allchains = list(set(self.allchains))
            #print(self.allchains)
            #self.calpha = [i.get_tip() for i in self.mol[0].get_residues()]
        else:
            self.calpha = list()
            for chainID in self.chain:
                for i in self.mol[0][chainID].get_residues():
                    self.calpha.append(i.get_calpha())
                    self.tip.append(i.get_tip())

        self.ids = [str(x.get_parent().get_id())+str(x.get_parent().get_parent().get_id()) for x in self.calpha]

        self.coords = list()
        for x in self.calpha:
            try:
                self.coords.append(x.get_location())
            except AttributeError:
                continue

        print("res count", len(self.calpha))
        self.coords = numpy.array(self.coords)
        self.coord_tips = [x.get_location() for x in self.tip]
        #print(len(self.coords), len(self.calpha))
        #self.coords = [x.get_location() for x in self.calpha]
        #print("Number of residues =", self.ids)
        return True

    def get_center(self, points):
        x_comp = list()
        y_comp = list()
        z_comp = list()
        for obj in points:
            x_comp.append(obj.get_location()[0])
            y_comp.append(obj.get_location()[1])
            z_comp.append(obj.get_location()[2])

        return numpy.array([numpy.mean(x_comp), numpy.mean(y_comp), numpy.mean(z_comp)])


    def chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def download(self):
        pdb_id = self.filename.split("/")[-1][:-4]
        chain = self.chain

        try:
            molecule.download_structure(pdb_id, ftype="pdb")
        except urllib.error.HTTPError:
            sys.exit("Failed to download structure for pdb id:", pdb_id)

        mol = molecule.load_structure(pdb_id+".pdb")
        element = mol[0].get_atoms()
        with open(pdb_id+"_pre.pdb", "w") as fp:
            for i in element:
                fp.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(i.get_id(),i.get_name(),i.get_parent().get_name(),i.get_parent().get_parent().get_id(),i.get_parent().get_id(),around(i.get_location()[0],decimals=3),around(i.get_location()[1],decimals=3),around(i.get_location()[2],decimals=3),i.get_occupancy(),i.get_bfactor(),'',i.get_element(),''))

        fdata = open(pdb_id+"_pre.pdb")
        with open(self.filename, "w") as fp:
            for lines in fdata:
                if "HETATM" in lines:
                    continue
                else:
                    fp.write(lines)
        os.remove(self.filename.split("/")[-1])
        os.remove(pdb_id+"_pre.pdb")
        return True


    def geo_center(self,coords):
        xi = 0
        yi = 0
        zi = 0
        for (x,y,z) in coords:
            xi = xi + x
            yi = yi + y
            zi = zi + z

        return [numpy.mean(x), numpy.mean(y), numpy.mean(z)]


        ############## Overlap ###############

    def turnPos(self, v):
        if v < 0:
            return -v
        else:
            return v

    def get_overlap(self, mode1, mode2):
        a = self.turnPos(numpy.dot(mode1, mode2))
        b = numpy.linalg.norm(mode1)*numpy.linalg.norm(mode2)
        return a/b


    def runTest(self):

        print(self.filename)

        if self.pdbid is not None:
            print("Downloading structure\n")
            self.download()
        elif os.path.exists(self.filename):
            print("\n",self.filename,"file exists!\n")
        else:
            sys.exit("Must either specify pdb id or provide a real input file")

        self.pull_pockets()
        print("\nComputing dynamics\n")
        self.prep_mol()

        try:
            self.active_site_index = [self.ids.index(x) for x in self.active_site if x[-1]]
        except:
            self.active_site_index = []

        model_normal = GNM(self.coords, self.calpha, self.tip, self.coord_tips, cutoff=self.cutoff)
        eigen_vectors_l, eigen_values_l = model_normal.calculate_modes()

        P_dict = dict()
        store_ratios = dict()
        norm_dict = dict()
        comp_dict = dict()
        psse_dict = dict()
        print("Collecting pocket dynamics data\n")
        toolbar_width = len(self.pocket_data)
        for numi, (key, val) in enumerate(self.pocket_data.items()):
            t = str((float(numi+1)/toolbar_width)*100)[:4]
            sys.stdout.write("\r%s%%" %t)
            sys.stdout.flush()

            try:
                #pocket_index = [self.ids.index(res) for res in val]
                pocket_index = val
            except:
                print("pocket ", key, "ddnt work\n")

            psse_dict[key] = self.psse_poc(eigen_values_l, eigen_vectors_l, pocket_index)

        print("\n")

        psse_values = list(psse_dict.values())
        for numval, (key, val) in enumerate(psse_dict.items()):
            psse_dict[key] = (val - numpy.mean(psse_values))/numpy.std(psse_values)

        return psse_dict

    def pull_pocket_features(self):
        fdata = open(self.filename[:-4] + "_out/" + self.filename.split("/")[-1][:-4] + "_info.txt")
        flag = False
        store_features = dict()
        for lines in fdata:
            if lines.strip() == "":
                flag = False

            if flag == True:
                feature = float(lines.strip().split()[-1])
                store_features["pocket" + pocket_name].append(feature)

            if "Pocket" in lines:
                flag = True
                pocket_name = lines.strip().split(" ")[1]
                store_features["pocket" + pocket_name] = list()

        return store_features

    def get_pockets(self):
        store_psse = self.runTest()
        all_pocket_features = self.pull_pocket_features()
        hydro_dict = {key: val[7] for key, val in all_pocket_features.items()}
        hydro_zscore = {key:(val-numpy.mean(list(hydro_dict.values())))/numpy.std(list(hydro_dict.values())) for key, val in hydro_dict.items()}
        store_means = dict()
        for key, val in store_psse.items():
            store_means[key] = numpy.mean([store_psse[key], hydro_zscore[key.split("/")[-1].rstrip("_atm.pdb")]])

        top_pockets_psse = sorted(store_means, key=store_means.get, reverse=True)
        store_pockets = dict()
        for i, pocket_id in enumerate(top_pockets_psse):
            mol = molecule.load_structure(pocket_id)
            element = mol[0].get_residues()
            store_pockets["Rank="+str(i+1)] = (pocket_id.split("/")[-1], store_means[pocket_id], list(set([str(x.get_id()) + x.get_parent().get_id() for x in element])))

        return store_pockets


    def psse_poc(self, eigen_values_l, eigen_vectors_l, pocket_index):

        model_pert = GNM(self.coords, self.calpha, self.tip, self.coord_tips, restrain=pocket_index, cutoff=self.cutoff)
        eigen_vectors2, eigen_values2 = model_pert.calculate_modes()

        eigen_vectors2 = list(eigen_vectors2)
        bestPair = dict()
        for i, mode1 in enumerate(numpy.array(eigen_vectors_l).T[1:10]):
            rem_mode = None
            bestOverlap = 0
            try:
                g=eigen_vectors2.pop(rem_mode)
            except TypeError:
                pass

            for j, mode2 in enumerate(numpy.array(eigen_vectors2).T[1:10]):
                O = self.get_overlap(mode1, mode2)
                if O > bestOverlap:
                    rem_mode = j
                    bestPair[i+1] = j+1
                    bestOverlap = O

            eigen_diff = dict()
            for key, val in bestPair.items():
                if key > 5:
                    break
                eigen_diff[key] = ((eigen_values2[val] - eigen_values_l[key])*100)/eigen_values_l[key]

            eigen_means = numpy.mean(list(eigen_diff.values()))

        return eigen_means








