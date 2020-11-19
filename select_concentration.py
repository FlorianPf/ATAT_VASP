#!/opt/local/bin/python
import argparse             # Parser for command-line options, arguments and sub-commands.                      https://docs.python.org/3/library/argparse.html#module-argparse
import sys                  # System-specific parameters and functions (e.g. for exiting after error).          https://docs.python.org/3/library/sys.html
import numpy as np          # Fundamental package for scientific computing.                                     https://numpy.org/doc/stable/

def main():
    parser = argparse.ArgumentParser(description='Select and output sqs of desired concentration.')

    parser.add_argument('-i', dest='input_file', type=str, default='sqs.out', required=False,
        help='Name of input file. Defaults to \'sqs.out\' if omitted.')
    parser.add_argument('-Nb', dest='num_nb', type=int, required=True,
        help='Desired number of Nb atoms in sqs.')
    parser.add_argument('-Ta', dest='num_ta', type=int, required=True,
        help='Desired number of Ta atoms in sqs.')
    parser.add_argument('-o', dest='output_file', type=str, default='sqs_sel.out', required=False,
        help='Name of output file. Defaults to \'sqs_sel.out\' if omitted.')

    args = parser.parse_args()

    num_nb = args.num_nb
    num_ta = args.num_ta

    input_file = [line.rstrip() for line in open(args.input_file, 'r').readlines()]
    output_file = open(args.output_file, 'w')

    count_sqs = 0
    count_matches = 0

    print("\n")

    structure = []
    nb_count = 0
    ta_count = 0
    num_backspace = 1
    sys.stdout.write("%s" %count_matches)
    sys.stdout.flush()
    for line in input_file:
        nb_count += line.count('Nb')
        ta_count += line.count('Ta')
        structure.append(line)
        if line == '':
            count_sqs += 1
            if nb_count == num_nb and ta_count == num_ta:
                count_matches += 1
                for content in structure:
                    output_file.write(content+"\n")
                if count_matches > 1:
                    num_backspace=int(np.floor(np.log10(count_matches-1))+1)
                sys.stdout.write("\b"*num_backspace+"%s" %count_matches)
                sys.stdout.flush()
            structure = []
            nb_count = 0
            ta_count = 0
    output_file.close()

    print(" of the {} provided structures had {} Nb and {} Ta atoms. They were saved to \'{}\'.\n".format(count_sqs, num_nb, num_ta, args.output_file))

if __name__ == "__main__":
    main()