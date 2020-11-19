#!/opt/local/bin/python
import sys                  # System-specific parameters and functions (e.g. for exiting after error).          https://docs.python.org/3/library/sys.html
import numpy as np          # Fundamental package for scientific computing.                                     https://numpy.org/doc/stable/
import argparse             # Parser for command-line options, arguments and sub-commands.                      https://docs.python.org/3/library/argparse.html#module-argparse

class Cluster:
    def __init__(self, data, coordinate_system):
        self.multiplicity = int(data[0])
        self.longest_length = float(data[1])
        self.num_atoms = int(data[2])

        self.atoms = []
        for i in range(self.num_atoms):
            coords_org = np.transpose(np.matrix([float(val) for val in data[3+i].split()[:3]]))
            coords_car = np.dot(np.transpose(coordinate_system), coords_org)
            self.atoms.append(np.transpose(coords_car).tolist()[0])
    
    def mean_distance(self):
        if self.num_atoms == 1:
            return None
        else:
            diff_sum = 0
            for i in range(self.num_atoms):
                for j in range(i+1, self.num_atoms):
                    diff_sum += np.sqrt(sum([(self.atoms[i][xyz]-self.atoms[j][xyz])**2 for xyz in range(3)]))
            return diff_sum*2/(self.num_atoms-1)/self.num_atoms     # Comprehend why the factor is indeed N!

def print_errors(output_file_name, rcorr, scorr, clusters):
    output_file = open(output_file_name, 'w')
    error_list = []
    num_backspace = 0
    for i in range(num_sqs):
        if i > 0:
            num_backspace = int(np.floor(np.log10(i))+1)
        sys.stdout.write("\b"*num_backspace+"%s" %(i+1))
        sys.stdout.flush()

        error = sum([clusters[j].multiplicity/(clusters[j].num_atoms*clusters[j].mean_distance())**damping*abs(scorr[i,j]-rcorr[i,j]) for j in range(1, num_clus)]) # The first cluster (index 0) is excluded, because it contains only one vertex. A correlation function for a single atom is not of interest.
        error_list.append(error)
        output_file.write(str(error)+"\n")
    print(" error functions corresponding to the input sqs where saved to \'"+output_file_name+"\'.\n")
    return error_list

def select_sqs(num_best, error_list):
    content = [line.rstrip() for line in open(args.sqs_file, 'r').readlines()]
    sqs_list = []
    single_sqs = []
    for line in content:
        single_sqs.append(line)
        if line == '':
            sqs_list.append(single_sqs)
            single_sqs = []

    output_file = open("best_sqs.out", 'w')
    for i in range(num_best):
        index = error_list.index(min(error_list))
        error_list.pop(index)
        sqs_temp = sqs_list.pop(index)
        for line in sqs_temp:
            output_file.write(line+"\n")
    output_file.close()

    print("The best {} sqs have been saved to \'best_sqs.out\'.\n".format(num_best))

def make_bin_chart(error_list):
    plt.bar(error_list)
    plt.savefig()

def restricted_float(x):
    if float(x) < 0.0:
        raise argparse.ArgumentTypeError("%r smaller than zero. Check -h for explanation of usage."%(x,))
    return float(x)

# p = argparse.ArgumentParser()
# p.add_argument("--arg", type=restricted_float)

def main():
    parser = argparse.ArgumentParser(description='Calculate error function.')

    parser.add_argument('-r', '--rcorr', dest='random_correlation_file', type=str, default='tcorr_finalRND.out',
        help='Name of file containing random correlations. Defaults to \'tcorr_finalRND.out\' if omitted.')
    parser.add_argument('-s', '--scorr', dest='sqs_correlation_file', type=str, default='tcorr_final.out',
        help='Name of file containing sqs correlations. Defaults to \'tcorr_final.out\' if omitted.')
    parser.add_argument('-c', '--clus', dest='clusters_file', type=str, default='clusters.out',
        help='Name of file containing clusters. Defaults to \'clusters.out\' if omitted.')
    parser.add_argument('-st', '--struc', dest='structure_file', type=str, default='lat.in',
        help='Name of the structure file. First three lines need to contain the coordinate system used by the calculations. Used for calculation of the cartesian atom positions. defaults to \'lat.in\'.')
    parser.add_argument('-d', '--damp', dest='damping_constant', type=float, default=2.0,
        help='Damping constant for weighting the impact of small and large clusters. Choose between 1 and 5. The higher, the less big clusters are taken into account. Defaults to 2.')
    parser.add_argument('-o', '--ofile', dest='output_file_name', type=str, default="errors.out",
        help='Name of output file to write the errors to. Defaults to \'errors.out\'.')
    parser.add_argument('-b', '--best', dest='num_best', type=restricted_float, default=0.,
        help='If this value is given as positive float, the best sqs get written to \'best_sqs.out\'. For num_best < 1 the best 100*num_best %% of sqs are chosen, for num_best >= 1 the round(num_best) best sqs get selected.')
    parser.add_argument('-sq', '--sqs', dest='sqs_file', type=str, default='sqs.out',
        help='File to take sqs from when using -b option. Defaults to \'sqs.out\'.')

    global args
    args = parser.parse_args()

    try:
        random_correlation_file = open(args.random_correlation_file, 'r')
    except:
        print("\nFile {} not found.".format(args.random_correlation_file))
        sys.exit()
    lines = random_correlation_file.readlines()
    global num_sqs
    num_sqs = len(lines)
    print("\nNumber of sqs: \t\t", num_sqs)
    global num_clus
    num_clus = len(lines[0].split())
    print("Number of clusters: \t{} (including the point cluster)\n".format(num_clus))

    rcorr = np.matrix([values for line in lines for values in line.split()], dtype=float).reshape((num_sqs, num_clus))

    try:
        sqs_correlation_file = open(args.sqs_correlation_file, 'r')
    except:
        print("\nFile {} not found.".format(args.sqs_correlation_file))
        sys.exit()
    lines = sqs_correlation_file.readlines()

    scorr = np.matrix([values for line in lines for values in line.split()], dtype=float).reshape((num_sqs, num_clus))

    try:
        structure_file = open(args.structure_file, 'r')
    except:
        print("\nFile {} not found.".format(args.structure_file))
        sys.exit()
    lines = structure_file.readlines()
    coordinate_system = np.matrix([values for vectors in lines[:3] for values in vectors.split()], dtype=float).reshape((3, 3))

    try:
        clusters_file = open(args.clusters_file, 'r')
    except:
        print("\nFile {} not found.".format(args.clusters_file))
        sys.exit()
    lines = clusters_file.readlines()
    clusters = []
    data = []
    for line in lines:
        if line != '\n':
            data.append(line)
        else:
            clusters.append(Cluster(data, coordinate_system))
            data = []

    global damping
    damping = args.damping_constant

    output_file_name = args.output_file_name

    error_list = print_errors(output_file_name, rcorr, scorr, clusters)

    num_best = args.num_best
    if num_best != 0:
        if num_best >= 1:
            select_sqs(round(num_best), error_list)
        else:
            select_sqs(round(num_best*num_sqs), error_list)

if __name__ == "__main__":
    main()

# Idea for improvement: Optional argument (number) for writing the best sqs to an additional output file.
# If number < 1: Take best number percent of sqs.
# If number >= 1: Take best round(number) sqs.

""" From the corrdump documentary:
    
    Cluster file format (clusters.out)

    for each cluster:
    [multiplicity]
    [length of the longest pair within the cluster]
    [number of points in cluster]
    [coordinates of 1st point] [number of possible species-2] [cluster function]
    [coordinates of 2nd point] [number of possible species-2] [cluster function]
    etc.

"""
