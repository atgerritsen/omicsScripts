import sys
import csv
import argparse
from collections import defaultdict

# read in tab-delimited file from R as arg1, fasta reference as arg2, and desired output as arg3


parser = argparse.ArgumentParser(description='A script to read a file of search targets and find the corresponding sequences.')
parser.add_argument('-i', '--input_path', type=str, nargs='?', dest="input",
                    help='path to tab-separated file of search targets')
parser.add_argument('-r', '--reference', type=str, nargs='?', dest="ref",
                    help='path to reference fasta file')
parser.add_argument('-o', '--output_path', type=str, nargs='?', dest="output",
                    help='output file name')


args = parser.parse_args()

columns = defaultdict(list)  # each value in each column is appended to a list

with open(args.input, 'r') as f:
    reader = csv.reader(f, delimiter="\t")  # read rows into a dictionary format
    reader.next()
    for row in reader:  # read a row as {column1: value1, column2: value2,...}
        for (i,v) in enumerate(row):  # go over each column name and value 
            columns[i].append(v)  # append the value into the appropriate list

geneID = columns[0]
# print(geneID)

f1 = open(args.ref, 'r')
f3 = open(args.output, 'w')


skip = 0
for line in f1:
    if line[0] == '>':
        _splitline = line.split(' ')
        accessorIDWithArrow = _splitline[0]
        accessorID = accessorIDWithArrow.lstrip('>')
        # print accessorID
        if accessorID in geneID:  # match DE genes in reference fasta and write if match
            f3.write(line)
            skip = 0
        else:
            skip = 1
    else:
        if not skip:
            f3.write(line)

print("Sequences extracted.")

f1.close()
f3.close()
