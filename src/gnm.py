
from packman import molecule
from functools import reduce
import numpy
from numpy import around
import os
import glob
import sys
import shutil
import urllib


class GNM():
    def __init__(self, coord, calpha, tip, coord_tip, cutoff=10.0, restrain=[]):
        self.calpha = calpha
        self.cutoff = cutoff
        self.coords = numpy.array(coord)
        self.restrain = restrain
        #print(self.restrain)
        
    def calculate_hessian(self, pf=False):
        
        gamma = 1.0
        n_atoms=len(self.coords)
        hessian=numpy.zeros((n_atoms, n_atoms), float)
        DMAT=numpy.ones((n_atoms, n_atoms), float)
        gamma_normal = 1.0
        gamma_poc = 10.0
        for i in range(len(self.coords)):
            diff = self.coords[i+1:, :] - self.coords[i]
            squared_diff = diff**2
            res_i = str(self.calpha[i].get_parent().get_id())+str(self.calpha[i].get_parent().get_parent().get_id())
            for j, s_ij in enumerate(squared_diff.sum(1)):
                
                res_j = str(self.calpha[j + i + 1].get_parent().get_id())+str(self.calpha[j + i + 1].get_parent().get_parent().get_id())
                
                if s_ij <= self.cutoff**2 or (res_i in self.restrain and res_j in self.restrain):
                    
                    if self.restrain != []:
                        if res_i in self.restrain and res_j in self.restrain:
                            use_gamma = gamma_poc
                        else:
                            use_gamma = gamma_normal
                    else:
                        use_gamma = gamma_normal
                        
                    diff_coords = diff[j]
                    j = j + i + 1
                    
                    hessian[i, j] = -use_gamma
                    hessian[j, i] = -use_gamma
                    hessian[i, i] = hessian[i, i] + use_gamma
                    hessian[j, j] = hessian[j, j] + use_gamma
                    DMAT[i, j] = s_ij
                    DMAT[j, i] = s_ij
        
        #if pf == True:
        #    hessian = numpy.divide(hessian, DMAT)
            
        self.hessian=hessian
        return True
        
    def calculate_decomposition(self): return numpy.linalg.eigh(self.hessian)

    def calculate_modes(self, n_modes=10):
        self.calculate_hessian()
        self.eigen_values, self.eigen_vectors = self.calculate_decomposition()
        return self.eigen_vectors, self.eigen_values
        

        

                

                


        
    
