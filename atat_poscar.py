#!/opt/local/bin/python
import os               # Miscellaneous operating system interfaces (e.g. for manipulating directories).    https://docs.python.org/3/library/os.html
import sys              # System-specific parameters and functions (e.g. for exiting after error).          https://docs.python.org/3/library/sys.html
import numpy as np      # Fundamental package for scientific computing.                                     https://numpy.org/doc/stable/
import argparse         # Parser for command-line options, arguments and sub-commands.                      https://docs.python.org/3/library/argparse.html#module-argparse

class Structure:
    def __init__(self, data, atom_list):
        try:                                            
            self.data = data
            pos = 0
            # The first three lines describe the coordinate system used, the vectors (defined by the lines) make up the rows of the matrix:
            self.coordinate_system = np.matrix([values for vectors in self.data[pos:pos+3] for values in vectors.split()], dtype=float).reshape((3, 3))
            pos += 3
            # The following three lines describe the lattice vectors of the unit cell:
            # The unit cell is expanded in most cases, therefore this generally not the identity matrix.
            self.lattice_vectors_org = np.matrix([values for vectors in self.data[pos:pos+3] for values in vectors.split()], dtype=float).reshape((3, 3))
            pos += 3
            self.lattice_vectors_car = self.lattice_vectors_org*self.coordinate_system
            
            # Determine the atom types present.
            self.atom_types = list(set([line.split()[-1] for line in self.data[pos:-1]]))  # beware: the slicing to -1 is to not read 'end'
            if atom_list is not None:
                if len(self.atom_types) == len(atom_list):
                    self.atom_types = atom_list
                else: raise ValueError("Length of order specification is different from number of atom types found in input file.")

            self.atom_numbers = [0 for atom_type in self.atom_types]
            atoms_unsorted = []
            # Read in all the atoms.
            for i in range(len(self.data[pos:-1])):
                atom_info = self.data[pos+i].split()
                atom_type = atom_info.pop()
                self.atom_numbers[self.atom_types.index(atom_type)] += 1
                coords_org = np.transpose(np.matrix([float(val) for val in atom_info]))
                atoms_unsorted.append(Atom(self, atom_type, coords_org))

            self.atoms = []
            # Getting the list of atoms in the right order (VASP needs this):
            for atom_type in self.atom_types:
                index=0
                for multiplicity in range(int(self.atom_numbers[self.atom_types.index(atom_type)])):
                    atom_found = False
                    while not atom_found:
                        if atoms_unsorted[index].atom_type==atom_type:
                            self.atoms.append(atoms_unsorted[index])
                            atom_found=True
                        index += 1

            # Output some structural information.
            if sqs_count == 1:
                print("\nStructural information (assuming all sqs are of the same composition):")
                info_list = [None]*(2*len(self.atom_types))
                info_list[::2] = self.atom_types
                info_list[1::2] = [str(num) for num in self.atom_numbers]
                print("Atom types (quantity):"+("  {} ({})"*len(self.atom_types)).format(*[info_list[i] for i in range(2*len(self.atom_types))])+"\n")
            
        except OSError as oserr:
            print("\nError while reading SQS {}: {}".format(sqs_count, oserr))
            sys.exit(oserr.errno)
        except ValueError as valerr:
            print("\nError while reading SQS {}: {}".format(sqs_count, valerr))
            sys.exit(-1)


class Atom:
    def __init__(self, struc, atom_type, coords_org):
        self.atom_type = atom_type
        self.coords_org = coords_org
        self.struc = struc
    
    # Transformation of the original coordinates to cartesian or direct.
    # One may insure oneself about these transformations by manually calculating the matrix operations.
    def org_to_car(self):
        transformation = np.transpose(self.struc.coordinate_system)
        return np.dot(transformation, self.coords_org)

    def org_to_dir(self):
        transformation = np.linalg.inv(np.transpose(self.struc.lattice_vectors_car))
        return np.dot(transformation, self.org_to_car())


def printPOSCAR(structure, name, representation):
    # For the exact format, check the VASP manual on POSCAR files: https://www.vasp.at/wiki/index.php/POSCAR
    poscar_dir = output_dir+"/poscar_"+str(sqs_count)
    os.mkdir(poscar_dir)
    os.chdir(poscar_dir)
    output_file = open("POSCAR", 'w')
    
    output_file.write(name+" (SQS number "+str(sqs_count)+")\n")                                                        # Name of structure (comment).
    output_file.write(" 1.0\n")                                                                                         # Scaling factor.
    for i in range(3):
        output_file.write(" "+("  {:12.8f}"*3).format(*[structure.lattice_vectors_car[i,j] for j in range(3)])+"\n")    # Three vectors describing the lattice.
    output_file.write(" "+("  {:<4s}"*len(structure.atom_types)).format(*structure.atom_types)+"\n")                    # Atom types in order of appearance.
    output_file.write(" "+("  {:>4d}"*len(structure.atom_numbers)).format(*structure.atom_numbers)+"\n")                # Corresponding quantities.

    if representation in ['car', 'cartesian']:
        output_file.write("Cartesian\n")        # With this tag present the atom positions get read in cartesian coordinates.
        for i in range(len(structure.atoms)):
            output_file.write(" "+("  {:20.16f}"*3).format(*np.transpose(structure.atoms[i].org_to_car()).tolist()[0])+"\n")
    else:
        output_file.write("Direct\n")           # With this tag present the atom positions get read in the coordinates defined by the lattice vectors.
        for i in range(len(structure.atoms)):
            output_file.write(" "+("  {:20.16f}"*3).format(*np.transpose(structure.atoms[i].org_to_dir()).tolist()[0])+"\n")
    output_file.close()
    
    os.chdir(output_dir)
    return

def main():
    # Necessities for command line usage:
    print("\nNote: Reading of the coordinate system as \'a, b, c, alpha, beta, gamma\' has not been implemented yet.")

    parser = argparse.ArgumentParser(description='Converts SQS output file (e.g. sqs.out) to VASP POSCAR files.')
    parser.add_argument('-i', '--ifile', dest='file_name', type=str, help='Name of input file.', required=False, default='best_sqs.out')
    parser.add_argument('-n', '--name', dest='name', type=str, help='Name of structure (first line of POSCAR file)', required=False, default='Comment (name of structure).')
    parser.add_argument('-o', '--order', nargs='*', dest='atom_list', help='Atom types in desired order (VASP calculations require matching order of atoms in POSCAR and POTCAR files).', required=False, default=['Li', 'Nb', 'Ta', 'O'])
    parser.add_argument('-r', '--repr', default='cartesian', type=str, choices=['car', 'cartesian', 'dir', 'direct'], required=False, help='Choose the representation of atom positions in the POSCAR file.', dest='representation')
    
    args = parser.parse_args()
    file_name = args.file_name
    name = args.name
    atom_list = args.atom_list
    representation = args.representation
    print("\nInput file name: ", file_name)

    global output_dir
    output_dir = os.getcwd()+"/output_files"
    global input_dir
    input_dir = os.getcwd()

    # Delete and/or rename previous results:
    if "output_files" in os.listdir(os.getcwd()):
        print("")
        if "output_files~" in os.listdir(os.getcwd()):
            for root, dirs, files in os.walk(os.getcwd()+"/output_files~", topdown=False):
                for file_to_del in files:
                    os.remove(os.path.join(root, file_to_del))
                for dir_to_del in dirs:
                    os.rmdir(os.path.join(root, dir_to_del))
            os.rmdir(os.getcwd()+"/output_files~")
            print("Deleted second to last results stored within \'/output_files~\'.")
        os.rename("output_files", "output_files~")
        print("Renamed previous results. \'/output_files\' is now \'/output_files~\'.")
    os.mkdir(output_dir)

    input_file = open(file_name, 'r')
    content = [line.rstrip() for line in input_file.readlines()]        # Create list with lines, truncated space at end
    input_file.close()
    os.chdir(output_dir)

    pos = 0
    global sqs_count
    sqs_count = 0
    data = []

    # Creating a separate structure / POSCAR file for each SQS found in the input file.
    num_backspace=0
    while pos in range(len(content)):
        line = content[pos]
        data.append(line)
        if line == 'end':
            sqs_count += 1
            printPOSCAR(Structure(data, atom_list), name, representation)
            data = []
            pos += 1
            # Live count of the number of converted structures:
            if sqs_count > 1:
                num_backspace=int(np.floor(np.log10(sqs_count-1))+1)
            sys.stdout.write("\b"*num_backspace+"%s" %sqs_count)
            sys.stdout.flush()
        pos += 1

        
    print(" SQS successfully converted to POSCAR file(s).\n")

if __name__ == "__main__":
    main()

"""
--------------------------
Manuals and documentaries:
--------------------------

VASP POSCAR files:      https://www.vasp.at/wiki/index.php/POSCAR
argparse documentary:   https://docs.python.org/3/library/argparse.html
mcsqs documentary:      https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/manual/node47.html
"""
