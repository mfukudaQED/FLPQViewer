import os
import pprint
import copy
import numpy as np

from parameters import BohrR

class Data_cube:
    """ class of data from md file """
    def __init__(self,
                 atomnum = None,
                 numbers = None,   # int
                 valence = None,   # int
                 positions_bohr = None,  # bohr, coordinate of atoms
                 positions_ang = None,  # Ang, coordinate of atoms
                 indices = None,
                 origin = None, # Ang
                 num_mesh = None,
                 mesh_vec = None, # Ang primitive vector for mesh
                 mesh_grid = None, # Ang
                 val = None,
                 small_mesh_grid = None, # Ang
                 small_val = None):

        pass

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"

    def read(self,filename):
        with open(filename, 'r') as file:

            # Skip the header lines
            for _ in range(2):
                next(file)

            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            #self.origin = [float(coord) for coord in line[1:4]]
            self.origin = [float(coord)*BohrR for coord in line[1:4]]


            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]


            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            self.indices = []
            for i in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])
                self.indices.append(i+1)

                          ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)

            self.val = np.empty((self.num_mesh[0]*self.num_mesh[1]*self.num_mesh[2])).tolist()

            # Read the file using the next function until the end of the file
            idx = 0
            while True:
                try:
                    line = next(file).strip().split()
                    for value in line:
                        self.val[idx] = float(value)
                        #print(self.val[idx])
                        idx += 1

                except StopIteration:
                    # End of file reached, break out of the loop
                    break

            ### Data structure ###
            #idx = 0
            #for i in range(self.num_mesh[0]):
            #    x = self.origin[0] + i*self.mesh_vec[0]
            #    for j in range(self.num_mesh[1]):
            #        y = self.origin[1] + j*self.mesh_vec[1]
            #        for k in range(num_mesh[2]):
            #            z = self.origin[2] + k*self.mesh_vec[2]
            #            print(x,y,z,self.val[idx])
            #            idx += 1
        return

    def read_header(self,filename):
        with open(filename, 'r') as file:

            # Skip the header lines
            for _ in range(2):
                next(file)

            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            self.origin = [float(coord)*BohrR for coord in line[1:4]]

            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]

            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            self.indices = []
            for i in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])
                self.indices.append(i+1)

            ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)

        return

    def write(self,filename):
        with open(filename, 'w') as file:
            # Write header lines
            file.write("Generated by Data_cube class\n")
            file.write("Data cube file format\n")

            # Write number of atoms and origin
            file.write(f"{self.atomnum} {self.origin[0]/BohrR:.6f} {self.origin[1]/BohrR:.6f} {self.origin[2]/BohrR:.6f}\n")

            # Write mesh info
            for i in range(3):
                num_mesh_out = self.num_mesh[i] if self.num_mesh[i]>0 else -self.num_mesh[i]
                file.write(f"{num_mesh_out} {self.mesh_vec[i][0]/BohrR:.6f} {self.mesh_vec[i][1]/BohrR:.6f} {self.mesh_vec[i][2]/BohrR:.6f}\n")

            # Write atomic symbols and positions_bohr
            for i in range(self.atomnum):
                file.write(f"{self.numbers[i]} {self.valence[i]} {self.positions_bohr[i][0]:.6f} {self.positions_bohr[i][1]:.6f} {self.positions_bohr[i][2]:.6f}\n")

            # Write data values
            for idx, value in enumerate(self.val):
                file.write(f"{value:.6e} ")
                if (idx + 1) % 6 == 0:
                    file.write("\n")
            if len(self.val) % 6 != 0:
                file.write("\n")
        return


# Add two cube files
def add_cube(cube1, cube2):
    """ Add two cube files """
    if cube1.num_mesh != cube2.num_mesh:
        raise ValueError("Number of mesh points do not match.")
    if cube1.origin != cube2.origin:
        raise ValueError("Origins do not match.")
    if cube1.mesh_vec != cube2.mesh_vec:
        raise ValueError("Mesh vectors do not match.")

    result_cube = copy.deepcopy(cube1)
    result_cube.val = [v1 + v2 for v1, v2 in zip(cube1.val, cube2.val)]
    return result_cube

# Subtract two cube files
def subtract_cube(cube1, cube2):
    """ Subtract two cube files """
    if cube1.num_mesh != cube2.num_mesh:
        raise ValueError("Number of mesh points do not match.")
    if cube1.origin != cube2.origin:
        raise ValueError("Origins do not match.")
    if cube1.mesh_vec != cube2.mesh_vec:
        raise ValueError("Mesh vectors do not match.")

    result_cube = copy.deepcopy(cube1)
    result_cube.val = [v1 - v2 for v1, v2 in zip(cube1.val, cube2.val)]
    return result_cube


################################
#if __name__ == '__main__':
#   # Example usage
#   filename_cube1=sys.argv[1]
#   filename_cube2=sys.argv[2]
#
#   cube1 = Data_cube()
#   cube1.read(filename_cube1)
#
#   cube2 = Data_cube()
#   cube2.read(filename_cube2)
#
#   added_cube = add_cube(cube1, cube2)
#   added_cube.write('added_cube.cube')
#
#   subtracted_cube = subtract_cube(cube1, cube2)
#   subtracted_cube.write('subtracted_cube.cube')
