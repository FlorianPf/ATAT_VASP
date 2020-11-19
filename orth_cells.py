#!/opt/local/bin/python
import numpy as np      # Fundamental package for scientific computing.                                     https://numpy.org/doc/stable/

def main():
    global tol
    tol = 10**(-8)

    lattice_file = [line.rstrip() for line in open("lat.in", 'r').readlines()]
    coordinate_system = np.matrix([values for vectors in lattice_file[:3] for values in vectors.split()], dtype=float).reshape((3, 3))

    sqscell_file = [line.rstrip() for line in open("sqscell.out", 'r').readlines()]
    num_sqscells = int(sqscell_file.pop(0))
    
    num_orth = 0
    
    for i in range(num_sqscells):

        lattice_vectors_org = np.matrix([values for vectors in sqscell_file[4*i+1:4*i+4] for values in vectors.split()], dtype=float).reshape((3, 3))
        lattice_vectors_car = lattice_vectors_org*coordinate_system

        a = lattice_vectors_car[0].reshape((3,1))
        b = lattice_vectors_car[1].reshape((3,1))
        c = lattice_vectors_car[2].reshape((3,1))

        if abs(np.transpose(a)*b)<tol and abs(np.transpose(b)*c)<tol and abs(np.transpose(c)*a)<tol:
            print(i)

if __name__ == "__main__":
    main()
