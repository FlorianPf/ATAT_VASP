#!/opt/local/bin/python
import os               # Miscellaneous operating system interfaces (e.g. for manipulating directories).    https://docs.python.org/3/library/os.html
import sys              # System-specific parameters and functions (e.g. for exiting after error).          https://docs.python.org/3/library/sys.html
import numpy as np      # Fundamental package for scientific computing.                                     https://numpy.org/doc/stable/
import argparse         # Parser for command-line options, arguments and sub-commands.                      https://docs.python.org/3/library/argparse.html#module-argparse

def my_length(vec):
    return np.sqrt(sum([value**2 for value in vec]))

def main():
    parser = argparse.ArgumentParser(description='Determine average lattice constant of random Li(NbTa)O3 alloy.')

    parser.add_argument('-i', '--ifile', dest='file_name', type=str, help='Name of input file (some POSCAR or CONTCAR etc.). Defaults to \'CONTCAR\' if omitted.', required=False, default='CONTCAR')
    parser.add_argument('-t', '--tol', dest='tolerance', type=float, help='Maximum mismatch between the lengths of the relaxed and unrelaxed lattice vectors. Defaults to 0.1 Angstroem if omitted.', required=False, default=0.5)
    
    args = parser.parse_args()

    tol = args.tolerance
    input_content = [line.rstrip() for line in open(args.file_name, 'r').readlines()]

    conc = float(input_content[6].split()[2])/24
    # print(conc)

    scale = float(input_content[1])
    lattice_vectors = [[scale*float(value) for value in line.split()] for line in input_content[2:5]]

    lengths = [my_length(vector) for vector in lattice_vectors]

    # print(lattice_vectors)
    # print()
    # print(lengths)
    # print()

    # index = lengths.index(min(lengths))
    # lengths.pop(index)
    # a = lattice_vectors.pop(index)
    # index = lengths.index(min(lengths))
    # lengths.pop(index)
    # b = lattice_vectors.pop(index)
    # c = lattice_vectors.pop(index)

    a, b, c = [lattice_vectors[lengths.index(length)] for length in sorted(lengths)]
    lengths = [my_length(vector) for vector in [a, b, c]]

    # print(a)
    # print(b)
    # print(c)
    # print()

    struc_not_found = "Structure not found, check tolerance or structure."

    if my_length(a)>5.1616-tol and my_length(a)<5.1616+tol:
        if my_length(b)>8.9401-tol and my_length(b)<8.9401+tol:
            if my_length(c)>27.8036-tol and my_length(c)<27.8036+tol:
                struc = 'long'
            else: print(struc_not_found)
        elif my_length(b)>13.9018-tol and my_length(b)<13.9018+tol:
            if my_length(c)>17.8802-tol and my_length(c)<17.8802+tol:
                struc = 'medium'
            else: print(struc_not_found)
        else: print(struc_not_found)
    elif my_length(a)>8.9401-tol and my_length(a)<8.9401+tol:
        if my_length(b)>10.3232-tol and my_length(b)<10.3232+tol and my_length(c)>13.9018-tol and my_length(c)<13.9018+tol:
            struc = 'short'
        else: print(struc_not_found)
    else: print(struc_not_found)

    # print("Structure is {}.".format(struc))
    # print()

    if struc == 'short':
        d_1 = my_length(b)/2
        d_2 = my_length([a[i]+b[i]/2 for i in range(3)])/2
        d_3 = my_length([a[i]-b[i]/2 for i in range(3)])/2
        # print(d_1, d_2, d_3)
        # print()

    if struc == 'medium':
        d_1 = my_length(a)
        d_2 = my_length([a[i]+c[i]/2 for i in range(3)])/2
        d_3 = my_length([a[i]-c[i]/2 for i in range(3)])/2
        # print(d_1, d_2, d_3)
        # print()

    if struc == 'long':
        d_1 = my_length(a)
        d_2 = my_length([a[i]+b[i] for i in range(3)])/2
        d_3 = my_length([a[i]-b[i] for i in range(3)])/2
        # print(d_1, d_2, d_3)
        # print()

    print(str(conc)+"\t"+str((d_1+d_2+d_3)/3))


    # atoms = [[float(val) for val in input_content[8+i].split()] for i in range(num_atoms)]

    # char=str(input_content[7])[0]
    # if char in ['D', 'd']:
    #     lattice = np.matrix([values for vectors in input_content[2:5] for values in vectors.split()], dtype=float).reshape((3, 3))
    #     # print(input_content[2:5])
    #     # print(lattice)
    #     # quit()
    #     positions_test = [np.dot(lattice, atom).tolist()[0] for atom in atoms]
    #     # print(positions)
    # else:
    #     positions_test = atoms
    #     # print(positions)
    # # print(atoms)

    # for i in range(num_atoms):
    #     print(np.sqrt(sum([(positions[i][k]-positions_test[i][k])**2 for k in range(3)])))



if __name__ == "__main__":
    main()


"""
Old version:

    args = parser.parse_args()

    input_content = [line.rstrip() for line in open(args.file_name, 'r').readlines()]

    # print(input_content[6])

    num_atoms=sum([int(val) for val in input_content[6].split()[:-1]])
    # print(input_content[6].split()[:-1])
    # print(num_atoms)
    # for i in range(num_atoms):
    #     print("atom nr. "+str(i)+":\t"+input_content[8+i])
    
    atoms = [[float(val) for val in input_content[8+i].split()] for i in range(num_atoms)]

    char=str(input_content[7])[0]
    if char in ['D', 'd']:
        lattice = np.matrix([values for vectors in input_content[2:5] for values in vectors.split()], dtype=float).reshape((3, 3))
        # print(input_content[2:5])
        # print(lattice)
        # quit()
        positions = [np.dot(lattice, atom).tolist()[0] for atom in atoms]
        # print(positions)
    else:
        positions = atoms
        # print(positions)
    # print(atoms)

    test=[]
    for i in range(num_atoms):
        for j in range(i+1,num_atoms):
            dist=np.sqrt(sum([(positions[i][k]-positions[j][k])**2 for k in range(3)]))
            if dist<6.5 and dist>3.5:
                test.append(dist)

    # print(np.mean(test))
    print(test)
    print(len(test))

"""