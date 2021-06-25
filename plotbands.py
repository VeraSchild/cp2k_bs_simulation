"""
Created Sep 28 2020
@author: Vera Schild
Based on the work of Alex McCab and Bram van der Linden

Creates a bandplot based on a .out file and .bs file. These files can be specified,
if not the newest/last modified files in the current folder will be used.

How to run example:
python plot_bands.py
- will use newest out- and bsfile
python plot_bands.py -o cp2k.out -b Sibulk.bs
- will use specified files, order in which they are called does not matter, you
- can also specify only one file.
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import re

parser = argparse.ArgumentParser(description="optional aguments")
parser.add_argument("--outfile", "-o")
parser.add_argument("--bandstructurefile", "-b")
args = parser.parse_args()

def main():
    # Check if files where specified, differently get newest out and/or bs file
    if args.outfile:
        cp2koutfile = args.outfile
    else:
        cp2koutfile = max(filter(lambda f: f.endswith('k.out'), os.listdir(os.path.abspath(os.getcwd()))), key=os.path.getctime)
    if args.bandstructurefile:
        bsfile = args.bandstructurefile
    else:
        bsfile = max(filter(lambda f: f.endswith('.bs'), os.listdir(os.path.abspath(os.getcwd()))), key=os.path.getctime)

    fermi_energy, special_kpoints, dist_special_kpoints, total_kpoints, project_name = readoutfile(cp2koutfile)
    energies, n_bands = readbs(bsfile)

    print("Plotting the bandstructure of {} for the special kpoints {}\nWith a fermi energy of {:.3f}eV and {} bands.".format(project_name, ' '.join([str(kpoint) for kpoint in special_kpoints]), fermi_energy, n_bands))

    if len(total_kpoints) > len(energies):
        total_kpoints = total_kpoints[:len(energies)]
    elif len(total_kpoints) < len(energies):
        energies = energies[:len(total_kpoints)]

    plotBands(energies, n_bands, dist_special_kpoints, special_kpoints, fermi_energy, total_kpoints, project_name)


# reading the needed information out of the outfile by searching for a line which
# has some of the words in it of which we know are on the same line as the information
def readoutfile(outfile):
    special_kpoints = []
    dist_special_kpoints = []
    skpoint_pos = []
    total_kpoints = 0

    with open(outfile) as f:
        for line in f:
            if "PROGRAM STARTED IN" in line:
                project_name = line.split()[-1].split("/")[-1]
            if "Number of K-Points in Set" in line:
                total_kpoints = float(line.split()[-1]) - 1
            if "Fermi energy:" in line:
                fermi_energy = float(line.split()[-1]) * 27.2113838565563
            if "Special K-Point" in line:
                special_kpoints.append(line.split()[4])
                skpoint_pos.append([float(i) for i in line.split()[5:]])

                # getting the distince from the kx, ky, kz of 2 points
                if len(dist_special_kpoints) == 0:
                    dist_special_kpoints.append(0)
                else:
                    dist_special_kpoints.append(np.sqrt((skpoint_pos[-1][0]-skpoint_pos[-2][0])**2 + (skpoint_pos[-1][1]-skpoint_pos[-2][1])**2 + (skpoint_pos[-1][2]-skpoint_pos[-2][2])**2))

    special_kpoints = ["$\Gamma$" if KPOINT == "GAMMA" else KPOINT for KPOINT in special_kpoints]

    # go from ex. dist special kpoints [0, 0.5, 0.75] to [0, 0.5, 1.25]
    for i in range(1, len(dist_special_kpoints)):
        dist_special_kpoints[i] = dist_special_kpoints[i] * (total_kpoints/(len(special_kpoints)-1)) + dist_special_kpoints[i-1]

    return fermi_energy, special_kpoints, dist_special_kpoints, total_kpoints, project_name

def restart(line):
    try:
        if int(line[1]) == 1:
            return True
    except:
        return False
# read the bsfile: reading out the energies and number of bands from the correct lines
def readbs(bsfile):
    energy = []
    energies = []

    with open(bsfile) as f:
        for line in f:

            line = line.split()

            if restart(line):
                energies = []

            if not any(c.isalpha() for c in line):
                if len(line) == 1:
                    if len(energy) != 0:
                        energies.append(energy)
                    else:
                        number_of_bands = int(line[0])
                    energy = get_energy(f, number_of_bands)

    return energies, number_of_bands

def get_energy(file, number_of_bands):
    lines_to_read = number_of_bands//4 + 1
    energy = []
    for lines in range(lines_to_read):
        energy += file.readline().split()

    return [float(e) for e in energy]

# Creating a list for x and y for the fermi line in the plot
def createFermiLine(x_axis_length,fermi_energy_ev ):
    x = np.arange(0,x_axis_length,0.01)
    y = np.empty(len(x))
    y.fill(fermi_energy_ev)

    return x, y

# Creating a scaled x axis
def createXAxis(k_paths, n_specialkpoints, total_kpoints):
    x = []
    for i in range(len(k_paths)-1):
        x.append(np.linspace(k_paths[i], k_paths[i+1], total_kpoints/(n_specialkpoints-1)))

    return np.array(x).reshape(1, len(x)*len(x[0]))[0]

def plotBands(energies, number_of_bands, kpoint_dis, special_kpoints, Ef, total_kpoints, project_name):
    ymin = 0
    ymax = 0

    x_axis = createXAxis(kpoint_dis, len(special_kpoints), total_kpoints)
    x_fermi, y_fermi = createFermiLine(max(x_axis), 0)

    # plot the line of every band
    for n in range(number_of_bands):
        band_energies = [energy[n] for energy in energies]
        re_aligned_energies = [energy - Ef for energy in band_energies]
        plt.plot(x_axis,re_aligned_energies, color='b', linewidth=1)

        if min(re_aligned_energies) < ymin:
            ymin = min(re_aligned_energies)
        if max(re_aligned_energies) > ymax:
            ymax = max(re_aligned_energies)

    plt.title(project_name)
    plt.xticks(kpoint_dis, special_kpoints)
    plt.ylabel('Energy (eV)')
    plt.rc('axes', labelsize=20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim(min(x_axis), max(x_axis))
    plt.ylim(ymin, ymax)
    plt.plot(x_fermi,y_fermi,'k--', linewidth=1)
    plt.savefig('Bandplot_%s.png'%project_name)
    plt.show()


main()
