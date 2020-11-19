#!/opt/local/bin/python
import argparse             # Parser for command-line options, arguments and sub-commands.                      https://docs.python.org/3/library/argparse.html#module-argparse
import numpy as np          # Fundamental package for scientific computing.                                     https://numpy.org/doc/stable/

def restricted_float(x):
    if float(x) < 0.0:
        raise argparse.ArgumentTypeError("%r smaller than zero. Check -h for explanation of usage."%(x,))
    return float(x)

def main():
    parser = argparse.ArgumentParser(description='Select the best sqs when given the errors and structures.')

    parser.add_argument('-e', dest='error_file', type=str, default='errors.out', required=False,
        help='Name of file containing the errors. Defaults to \'errors.out\' if omitted.')
    parser.add_argument('-b', dest='num_best', type=restricted_float, required=True,
        help='If this value is given as positive float, the best sqs get written to \'best_sqs.out\'. For num_best < 1 the best 100*num_best %% of sqs are chosen, for num_best >= 1 the round(num_best) best sqs get selected.')
    parser.add_argument('-s', dest='sqs_file', type=str, default='sqs.out',
        help='File to take sqs from when using -b option. Defaults to \'sqs.out\' if omitted.')
    parser.add_argument('-o', dest='output_file', type=str, default='best_sqs.out', required=False,
        help='Name of output file to write best sqs to. Defaults to \'best_sqs.out\' if omitted.')

    args = parser.parse_args()

    num_best = 0
    arg_best = args.num_best
    if arg_best != 0:
        if arg_best >= 1:
            num_best = round(arg_best)
        else:
            num_best = round(num_best*num_sqs)

    error_list = open(args.error_file, 'r').readlines()
    sqs_list = []
    single_sqs = []
    for line in open(args.sqs_file, 'r').readlines():
        single_sqs.append(line)
        if line.rstrip() == '':
            sqs_list.append(single_sqs)
            single_sqs = []

    output_file = open(args.output_file, 'w')
    for i in range(num_best):
        index = error_list.index(min(error_list))
        error_list.pop(index)
        sqs_temp = sqs_list.pop(index)
        for line in sqs_temp:
            output_file.write(line)
    output_file.close()

    print("\nThe best {} sqs have been saved to \'{}\'.\n".format(num_best, args.output_file))

if __name__ == "__main__":
    main()