"""
This script is used for adding energy functions according to user-defined base pairs,
Programmed by Jun-Jie Zhu and Zheng-Xin Li,
Energy functions fitted by Zheng-Xin Li and Ge Song,
Feel free to contact us Nucleotide AlchemiSts (NASs) if you have any problems!
"""

import argparse
import re

# parameters to define
parser = argparse.ArgumentParser()
parser.add_argument('--gro', '-g', default='2a43.gro',
                    help='name (path) for input gromacs file')
parser.add_argument('--pair', '-p', default='G3-C16,C4-G15,G5-C14,G6-C13,G11-C26,C10-G27,C9-G28',
                    help='user-defined base pairs')
parser.add_argument('--out', '-o', default='fix.txt',
                    help='name (path) for output file')
args = parser.parse_args()

# atom pair info. for specific base pairs
AU_pair = [['N1', 'H61'], ['H3', 'O4']]
GC_pair = [['H21', 'H1', 'O6'], ['O2', 'N3', 'H41']]
GU_pair = [['O6', 'H1'], ['H3', 'O2']]


#####
# energy functions for base pairs
AU_energy = ['(418.4*(-0.475)*(0.19-x)^2)*step(x-0.19)*step(0.49-x)+(418.4*(-0.475)*0.3*(2*x-0.68))*step(x-0.49)',
             '(418.4*(0.475)*10*(0.19-x)^2)*step(x-0.19)*step(0.22-x)+(418.4*(0.475)*10*0.03*(2*x-0.41))*step(x-0.22)',
             ]
             
AU_energy = ['(-198.74*(0.19-x)^2)*step(x-0.19)*step(0.49-x)+(-59.622*(2*x-0.68))*step(x-0.49)',
             '(1987.4*(0.19-x)^2)*step(x-0.19)*step(0.22-x)+(59.622*(2*x-0.41))*step(x-0.22)',
             ]         

GC_energy = ['(418.4*(-0.32)*(0.2-x)^2)*step(x-0.2)*step(0.5-x)+(418.4*(-0.32)*0.3*(2*x-0.7))*step(x-0.5)',
             '(418.4*0.32*2*(0.2-x)^2)*step(x-0.2)*step(0.35-x)+(418.4*0.32*2*0.15*(2*x-0.55))*step(x-0.35)'
             ]
             
GC_energy = ['(-133.888*(0.2-x)^2)*step(x-0.2)*step(0.5-x)+(-40.1664*(2*x-0.7))*step(x-0.5)',
             '(267.776*(0.2-x)^2)*step(x-0.2)*step(0.35-x)+(40.1664*(2*x-0.55))*step(x-0.35)'
             ]             

GU_energy = ['(418.4*(-0.45)*(0.2-x)^2)*step(x-0.2)*step(0.45-x)+(418.4*(-0.45)*0.25*(2*x-0.65))*step(x-0.45)',
             '(418.4*0.45*2.5*(0.2-x)^2)*step(x-0.2)*step(0.3-x)+(418.4*0.45*2.5*0.1*(2*x-0.5))*step(x-0.3)'
             ]

GU_energy = ['(-188.28*(0.2-x)^2)*step(x-0.2)*step(0.45-x)+(-47.07*(2*x-0.65))*step(x-0.45)',
             '(470.7*(0.2-x)^2)*step(x-0.2)*step(0.3-x)+(47.07*(2*x-0.5))*step(x-0.3)'
             ]
# process input base pair info.
pair_list = []
pair_info = []
for items in args.pair.split(','):
    tmp_pair_list, tmp_info_list = [], []
    for bps in items.split('-'):
        # identify the input format
        matches = re.match(r"([A-Za-z]+)(\d+)", bps)
        reverse_matches = re.match(r"(\d+)([A-Za-z]+)", bps)
        if matches:
            alphabetic_part = matches.group(1)
            numeric_part = matches.group(2)

            tmp_pair_list.append(alphabetic_part)
            tmp_info_list.append(numeric_part + alphabetic_part)
        elif reverse_matches:
            numeric_part = matches.group(1)
            alphabetic_part = matches.group(2)

            tmp_pair_list.append(alphabetic_part)
            tmp_info_list.append(numeric_part + alphabetic_part)
        else:
            raise ValueError('The format of given base pair info. is unsupported now.')

    if (tmp_pair_list[0] == 'G' and tmp_pair_list[1] == 'C') \
            or (tmp_pair_list[0] == 'G' and tmp_pair_list[1] == 'U') \
            or (tmp_pair_list[0] == 'A' and tmp_pair_list[1] == 'U'):
        pair_list.append(tmp_info_list)
        pair_info.append(tmp_pair_list[0] + tmp_pair_list[1])
    elif (tmp_pair_list[0] == 'C' and tmp_pair_list[1] == 'G') \
            or (tmp_pair_list[0] == 'U' and tmp_pair_list[1] == 'G') \
            or (tmp_pair_list[0] == 'U' and tmp_pair_list[1] == 'A'):
        pair_list.append(tmp_info_list[::-1])
        pair_info.append(tmp_pair_list[1] + tmp_pair_list[0])
    else:
        raise ValueError('The base pair you give is unsupported now.')

# find atom pairs according to the input base pair
atom_pair_info = []
for i in range(len(pair_info)):
    former = []
    latter = []
    with open(args.gro, 'r') as f:
        for lines in f.readlines():
            data = lines.split()
            if data[0] == pair_list[i][0] and data[1] in globals()[pair_info[i] + '_pair'][0]:
                former.append(data[2])
            elif data[0] == pair_list[i][1] and data[1] in globals()[pair_info[i] + '_pair'][1]:
                latter.append(data[2])
    atom_pair_info.append([former, latter])

# check if all atom pairs exist
for bps in range(len(atom_pair_info)):
    for base in range(2):
        if len(atom_pair_info[bps][base]) == 0:
            raise ValueError('The base pair %s-%s you give is not found in the given gromacs file.' %
                             (pair_list[bps][0], pair_list[bps][1]))

index = 0
with open(args.out, 'w') as f:
    for i in range(len(atom_pair_info)):
        for j in range(len(atom_pair_info[i][0])):
            index += 1
            f.write('f%d: DISTANCE ATOMS=%s,%s\n' % (index, atom_pair_info[i][0][j], atom_pair_info[i][1][j]))
            f.write('f_f%da: MATHEVAL ARG=f%d FUNC=%s PERIODIC=NO\n' %
                    (index, index, globals()[pair_info[i] + '_energy'][0]))
            f.write('f_f%db: MATHEVAL ARG=f%d FUNC=%s PERIODIC=NO\n' %
                    (index, index, globals()[pair_info[i] + '_energy'][1]))
    f.write('BIASVALUE ARG=')
    for i in range(index):
        f.write('f_f%da,f_f%db' % (i + 1, i + 1))
        if i != index - 1:
            f.write(',')
        else:
            f.write(' LABEL=b_f')
