"""
Plots electronic PDOS
"""

import argparse as ap
import os
import sys
from xml.etree import ElementTree as ET

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ticker
import numpy as np
import pylab as pyl
import re
from readrawdata import readrawdata
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

def option_parser():
    """
    Parse command line inputs

    Parse:
        directory
        files
        eboundary
        assign_type
        phonon
        no_legend
        scale_box
        dos_boundary
        eunit
        funit
        title
        no_title
        legend_position
        scale_box
        reverse_colors
        reverse_plots
        no_fill
        no_reverse_spin
        max_ticks_x
        max_ticks_y
        legend_label
        grid
        show

    :return input_options: Dictionary of parsed command line arguments
    """
    p = ap.ArgumentParser(description=\
                'Plot single and multiple electronic/phonon density of sataes.')

    help_directory = 'List of the directories in which the data to be plotted have to be found.    If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'

    help_files = 'List of DOS files to plot. If this argument is not specified, defaults to TDOS. Useful for plots of projected DOS.'

    help_eboundary = 'One or two floats corresponding to the minimum and maximum energy in the plot of the electronic DOS in the units specified by --eunit, respectively. If the argument --phonon is present, the option values specify the minimum and maximum frequency in cm-1 appearing in the phonon-dispersion plot. If either the argument --eboundary is not present or no floats are given, the energy (frequency) boundaries are chosen correspondingly to the maximum and minimum of the full data to be plotted.'

    help_dos_boundary = 'One or two floats corresponding to the minimum and maximum density in the plot of the density of states (DOS) according to the units specified by --eunit, respectively. If the argument --phonon is present, the option values specify the minimum and maximum DOS, according to the units specified by --funit. If either the argument --dos_boundary is not present or no floats are given, the DOS boundaries are chosen correspondingly to the maximum and minimum of the full data to be plotted. Please notice that in the case of spin-polarized system the electronic spin-down DOS is represented by default by negative values.'

    help_assign_type = "List of the description keys for each plot of the electronic DOS. Possible choices are 'KS' (standard Kohn-Sham calculation), 'GW' (G0W0 calculation), and 'WA' (data are interpolated by using Wannier functions). If not present, the option 'KS' is assumed for all plots. Not used if the argument --phonon is present."

    help_eunit = "Set the units of the energy appearing in the plot of the electronic density of states. Possible choices are 'eV' (electronvolt, default) and 'Ha' (Hartree)."

    help_funit = "Set the units of the frequency appearing in the plot of the phonon density of states. Possible choices are 'icm' (inverse centimeter, cm^-1, default), 'meV' (millielectronvolt), and 'THz' (terahertz)."

    help_legend_position = "The location of the legend. The strings 'upper left', 'upper right', 'lower left', 'lower right' place the legend at the corresponding corner of the axes/figure. The strings 'upper center', 'lower center', 'center left', 'center right' place the legend at the center of the corresponding edge of the axes/figure. The string 'center' places the legend at the center of the axes/figure. The string 'best' places the legend at the location, among the nine locations defined so far, with the minimum overlap with other drawn artists. This option can be quite slow for plots with large amounts of data; your plotting speed may benefit from providing a specific location. For back-compatibility, 'center right' (but no other location) can also be spelled 'right', and each string locations can also be given as the corresponding numeric value."

    help_phonon = 'If present, it tags the plotting of the phonon DOS. If absent, the electronic DOS is plotted.'

    help_no_legend = 'If present, it disables the plotting of the legend.'

    help_no_fill = 'If present, plots are not filled.'

    help_title = "Used as --title 'String as a title' assign a title to the plot."

    help_no_title = 'If present, it disables the writing of the title.'

    help_scale_box = "One or two floats corresponding to the scaling factor in the horizontal and vertical size of the plot appearence, respectively."

    help_reverse_colors = "If present, the order of the sequence of colors of the plots is reversed."

    help_reverse_plots = "If present, the order of appearance of the plots is reversed."

    help_no_reverse_spin = "If present and a spin-polarized density of state is plotted, the spin-down DOS is represented by positive values instead of negative ones."

    help_max_ticks_x = "Specifies the maximum number of ticks along the x-axis in the plot."

    help_max_ticks_y = "Specifies the maximum number of ticks along the y-axis in the plot."

    help_legend_label = "Specifies the labels to appear in the legend for each plot."

    help_grid = 'If present, a grid is plotted in correspondence to the position of the major ticks.'

    help_show = "Opens plot in new window"

    help_norbitals = "Number of orbital contributions to plot for QP weights. Defaults to 3 (s,p,d) and must be at least 1."
    #---------------------------------------------------------------------------

    p.add_argument('-d','--directory',
                   nargs = '*', default = ["."],
                   type = str, help = help_directory)

    p.add_argument('-e','--eboundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_eboundary)

    p.add_argument('-f', '--files',
                    nargs='*', default=None,
                    type=str, help=help_files)

    p.add_argument('-db','--dos_boundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_dos_boundary)

    p.add_argument('-s','--scale_box',
                   nargs = '*', default = [1.0, 1,0],
                   type = float, help = help_scale_box)

    p.add_argument('-a','--assign_type',
                   nargs = '*', default = ['KS'],
                   choices = ['KS', 'GW', 'WA'],
                   type = str, help = help_assign_type)

    p.add_argument('-p','--phonon', action='store_true', help = help_phonon)

    p.add_argument('-nf','--no_fill', action='store_true', help = help_no_fill)

    p.add_argument('-rc','--reverse_colors', action='store_true', help = help_reverse_colors)

    p.add_argument('-rp','--reverse_plots', action='store_true', help = help_reverse_plots)

    p.add_argument('-nl','--no_legend', action='store_true', help = help_no_legend)

    p.add_argument('-nrs','--no_reverse_spin', action='store_true', help = help_no_reverse_spin)

    p.add_argument('-eu','--eunit',
                   type = str, help = help_eunit,
                   choices = ['eV', 'Ha'], default = 'eV')

    p.add_argument('-fu','--funit',
                   type = str, help = help_funit,
                   choices = ['icm', 'meV', 'THz'], default = 'icm')

    p.add_argument('-t','--title',
                   type = str, default = None, help = help_title)

    p.add_argument('-nt','--no_title', action='store_true', help = help_no_title)

    p.add_argument('-mtx','--max_ticks_x',
                   type = int, default = None, help = help_max_ticks_x)

    p.add_argument('-mty','--max_ticks_y',
                   type = int, default = None, help = help_max_ticks_y)

    p.add_argument('-ll','--legend_label',
                   nargs = '*', default = [],
                   type = str, help = help_legend_label)

    p.add_argument('-l','--legend_position',
                   type = str, help = help_legend_position,
                   choices = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                              'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'],
                   default = 'best')

    p.add_argument('-g','--grid', action='store_true', help = help_grid)

    p.add_argument('-sh', '--show',
                   action='store_true',
                   help=help_show)
    #p.add_argument('-shx','--share_x', default="col", help = help_sharex)
    #p.add_argument('-shy','--share_y', default="col", help = help_sharey)
    p.add_argument('-no', '--norbitals', type=int, default=3, help=help_norbitals)

    #---------------------------------------------------------------------------

    args = p.parse_args()
    input_options = {}

    input_options['directory'] = args.directory
    input_options['files'] = args.files
    input_options['atype'] = args.assign_type

    input_options['emin'] = None
    if len(args.eboundary) >= 1: input_options['emin'] = args.eboundary[0]
    input_options['emax'] = None
    if len(args.eboundary) >= 2: input_options['emax'] = args.eboundary[1]

    input_options['dmin'] = None
    if len(args.dos_boundary) >= 1: input_options['dmin'] = args.dos_boundary[0]
    input_options['dmax'] = None
    if len(args.dos_boundary) >= 2: input_options['dmax'] = args.dos_boundary[1]

    input_options['sx'] = 1.0
    if len(args.scale_box) >= 1: input_options['sx'] = args.scale_box[0]
    input_options['sy'] = 1.0
    if len(args.scale_box) >= 2: input_options['sy'] = args.scale_box[1]

    input_options['phonon'] = args.phonon

    input_options['title'] = args.title
    input_options['no_title'] = args.no_title

    input_options['eunit'] = args.eunit
    input_options['funit'] = args.funit

    input_options['maxticksx'] = args.max_ticks_x
    input_options['maxticksy'] = args.max_ticks_y

    input_options['legend'] = args.legend_label

    input_options['no_legend'] = args.no_legend
    input_options['leg_pos'] = args.legend_position
    if len(args.legend_position)<=2:
        input_options['leg_pos'] = int(args.legend_position)

    input_options['no_reverse_spin'] = args.no_reverse_spin
    input_options['reverse_plots'] = args.reverse_plots
    input_options['reverse_colors'] = args.reverse_colors
    input_options['no_fill'] = args.no_fill
    input_options['grid'] = args.grid

    input_options['show'] = args.show
    input_options['norbitals'] = args.norbitals

    return input_options


def find_steps(infile):
    inf = open(infile,"r")
    steps = 0
    while True:
       line=inf.readline().strip().split()
       if len(line)==0: break
       steps += 1
    inf.close()
    return steps
def read_evalcore() -> dict:
    """
    Extracts the species and their numbers from EVALCORE.OUT
    """
    #all_species_file = os.listdir(f"{os.environ['EXCITINGROOT']}/species")
    all_species = []
    try:
        lines = []
        species = []
        species_index = {}
        with open('EVALCORE.OUT') as f:
            for line in f:
                if not line.startswith('Species'):
                    continue
                lines.append(line)
        f.close()
        for line in lines:
            information = line.split(" ").remove(":")
            #if any(species in line for species in all_species):
            species_index[re.search(r'\d+',line).group()] =line[line.find("(")+1:line.find(")")]
        return species_index

    except OSError:
        print('Cannot open EVALCORE.OUT.')


def determine_species(filename: str) -> str:
    """
    Determines the species of PDOS_Sss_Aaaaa.OUT files for QP contribution
    """
    file_name = filename.strip(".OUT").split("_")
    species = read_evalcore()
    species_index = file_name[1].strip("S")
    # Indices in evalcore have no leading zero, convert to int then str
    return species[str(int(species_index))]

def get_orbital_contribution(infile: list[str], norbitals: int):
    """
    Input:
        infile: list of strings of filenames with DOS of the same species.

    Description:
    - Determine DOS orbital contribution from given files of same species.
    - Given PDOS_Sss_Aaaaa.OUT, each block of nwdos rows determine PDOS contributions from norbitals orbitals by summing up the corresponding blocks from each Aaaaa file for a given Sss
    - First block is assumed s. Next 3 blocks are 2px,2py,2pz. Next 5 are the five 3d orbitals etc
    - Columns of the form Energy | PDOS

    Returns:
    Tuple of np array containing total PDOS contribution from each orbital of a given species.
    """
    norbitals_per_subshell = [1,3,5,7,9]
    nblocks = sum(norbitals_per_subshell[:norbitals])
    orbital_contributions = []
    for ifile in infile:
        dos = readrawdata(ifile)[0]
       # print(len(readrawdata(ifile)))
        # only 2 columns
        # need to extract the number of blocks according to norbitals
        #block_size = find_steps(ifile)
        orbital_data = dos[:nblocks]
        #print("wtf", orbital_data)
        orbital_contributions.append(orbital_data[:,1])
    # print(orbital_contributions)
    # for i in range(len(orbital_contributions[0])):
    #     for row in orbital_contributions[0]:
    s_contrib = np.zeros_like(orbital_contributions[0][0])
    p_contrib = np.zeros_like(orbital_contributions[0][0])
    d_contrib = np.zeros_like(orbital_contributions[0][0])
    for contrib in orbital_contributions:
        s_contrib += contrib[0]
        p_contrib += contrib[1]
        p_contrib += contrib[2]
        p_contrib += contrib[3]
        d_contrib += contrib[4]
        d_contrib += contrib[5]
        d_contrib += contrib[6]
        d_contrib += contrib[7]
        d_contrib += contrib[8]
    # total_contrib = np.sum(orbital_contributions, axis=0)
    #print(total_contrib)
    #print("total_contrib: ", total_contrib)
    return [s_contrib, p_contrib, d_contrib]

def main(input_options):
    '''
    input:
    :input_options: dictionary that holds the input options parsed
                    from the command line arguments
    '''
    directory = input_options['directory']
    files = input_options['files']
    atype = input_options['atype']
    phonon = input_options['phonon']
    emin = input_options['emin']
    emax = input_options['emax']
    dmin = input_options['dmin']
    dmax = input_options['dmax']
    eunit = input_options['eunit']
    funit = input_options['funit']
    sx = input_options['sx']
    sy = input_options['sy']
    title = input_options['title']
    no_title = input_options['no_title']
    leg_pos = input_options['leg_pos']
    no_leg = input_options['no_legend']
    no_fill = input_options['no_fill']
    reverse_colors = input_options['reverse_colors']
    reverse_plots = input_options['reverse_plots']
    no_reverse_spin = input_options['no_reverse_spin']
    maxticksx = input_options['maxticksx']
    maxticksy = input_options['maxticksy']
    legend = input_options['legend']
    grid = input_options['grid']
    show = input_options['show']
    norbitals = input_options['norbitals']
    ha2ev = 27.211396132
    #print(files)
    sorted_files = []
    read_species = read_evalcore()
    for i in range(0,len(read_species)):
        sorted_files.append([])
    #print(files)
    for key in read_species:
        temp = ""
        if len(key) == 1:
            temp = "S0"+key
        else:
            temp = "S"+key
        for j in range(0,len(files)):
            if re.search(rf'{temp}',files[j]):
                sorted_files[int(key) - 1].append(files[j])
            else:
                continue
    dos_data = []
    print("wtf", sorted_files)
    ordered_species = []
    for filenames in sorted_files:
        # Get correct order of species to label plots
        if len(filenames) > 0:
            ordered_species.append(determine_species(filenames[0]))
   # #
   #for atomfiles in sorted_files:
    #    data = get_orbital_contribution(atomfiles)
     #   extract_rel_orbitals = data[:,:norbitals]
        # Should be length of nspecies
        # qp_data dimensions: N_s x  N_b x (norbitals) x N
      #  qp_data.append(extract_rel_orbitals)


    if (len(files) > norbitals * len(dos_data)):
        print("WARNING: Number of plots greater than available subplot space. Please change the number of files you are plotting, or norbitals.")
    relevant_data = []
    for files in sorted_files:
        data = get_orbital_contribution(files, norbitals)
        print("??", data)
        relevant_data.append(data)

    #print("len: ", len(relevant_data))
    size_title = "40"

    elab = 'Energy$-E_F$ ['+eunit+']'
    if phonon:
        elab = 'Frequency [cm$^{-1}$]'
        if funit!= 'icm': elab = 'Frequency [' + funit + ']'

    dlab = 'DOS [states/'+eunit+'/unit cell]'
    if phonon:
        dlab = 'Phonon DOS [states/cm$^{-1}$]'
        if funit!= 'icm':
            dlab = 'Phonon DOS [states/'+funit+']'
            if funit== 'THz': dlab = 'Phonon DOS [states/$\,$' + funit + ']'

    line_thickness = 3.0
    sline_thickness = 2.0
    axes_thickness = 4.0
    leg_size = 30

    dpi = 300

    figcolor = 'white'

    line_color = ["lightcoral", "yellow", "blue", "lightgreen", "green", "turquoise", "darkslategrey", "cyan", "darkgoldenrod"]
    fill_color = ["cornflowerblue", "lightsalmon", "lightgreen", "moccasin"]
    fig = plt.figure(figsize=(10,5),dpi=dpi)
    fig.patch.set_edgecolor(figcolor)
    fig.patch.set_facecolor(figcolor)

    plt.rcParams['axes.linewidth']  = axes_thickness # set the value globally
    plt.rcParams['grid.linewidth']  = 1.5
    plt.rcParams['xtick.labelsize'] = 24
    plt.rcParams['ytick.labelsize'] = 24
    plt.rcParams['axes.edgecolor']  = 'black'
    plt.rcParams['axes.labelsize']  = 18      # fontsize of the x any y labels
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.axisbelow']  = 'True'  # axis gridlines and ticks are below
                                              # the axes elements (lines, text, etc)
    plt.rcParams['legend.fontsize'] = leg_size
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10
    
    ax = plt.gca()
    ax.set_xlim([emin, emax])
    ax.set_ylim([dmin, dmax])
    ax.set_ylabel(dlab)
    ax.set_xlabel(elab)
    ax.xaxis.grid(True, which='minor')
    ax.yaxis.grid(True, which='minor')
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.tick_params(axis='y', which='minor', bottom=True)
    orbital_label=["s", "p", "d", "f"]
    plt.rcParams.update({'mathtext.default':'regular'})
    counter = 0
    #xaxis = readrawdata(species[0])[0][0][0] * ha2ev
    for i,species in enumerate(sorted_files):
        xaxis = readrawdata(species[0])[0][0][0] * ha2ev
        data = get_orbital_contribution(species, norbitals)
        for j,orbitals in enumerate(data):
            #print("orb shape: ", orbitals)
            #print()
            plt.plot(xaxis, orbitals/ha2ev, color=line_color[counter],label=f"{ordered_species[i]} {orbital_label[j]}" )
            counter += 1

    ax.legend(fontsize=10, loc="upper right", ncol = 3)
    output_formats=["png", "pdf", "eps"]
    for format in output_formats:
        if format == "png":
            fig.savefig(f"PLOT.{format}", format="png", dpi=300, bbox_inches='tight')
        else:
            fig.savefig(f"PLOT.{format}",format=format, bbox_inches="tight")
    sys.exit()
if __name__ == '__main__':
    input_options = option_parser()
    main(input_options)
