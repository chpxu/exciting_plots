
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import sys
import matplotlib         as mpl
import matplotlib.pyplot  as plt
import numpy              as np
import os
import argparse as ap
import re
from readrawdata import readrawdata

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

ha2ev = 27.211396132 ## conversion hartree -> eV
bo2an = 0.52917721067



def option_parser():
    """
    Parse command line inputs for exciton and QP weights

    Parse:
        directory
        files
        scale
        ymin
        ymax
        title
        grid
        path
        qp
        eunit
    :return input_options: Dictionary of parsed command line arguments
    """
    p = ap.ArgumentParser(description=\
                'Plot exciton weight data from single and multiple files as subplots.')

    help_directory = 'List of the directories in which the data to be plotted have to be found. If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'

    help_files = 'List of file names containing data to plot. At least one file must be specified.'

    help_title = "Used as --title 'String as a title' assign a title to the plot."
    help_no_title = 'If present, it disables the writing of the title.'

    help_grid = 'If present, a grid is plotted in correspondence to the position of the major ticks.'

    help_scale = "Array of float corresponding to the scaling factor in the horizontal and vertical size of the exciton circles for each plot."

    help_ymin = "Minimum y-value to output from."
    help_ymax = "Maximum y-value to output from."

    #help_nrows = "Number of rows in subplot."
    #help_ncols = "Number of columns in subplot."
    help_sharex = "Whether the x-axes in subplot should share the same domain."
    help_sharey = "Whether the y-axes in subplot should share the same domain."
    help_path = "Path through Brillouin Zone to print on the x-axis of bandstructure. Default is 'hexagonal' which corresponds to the tuple ('W','L','$\Gamma$','X','W','K')."
    help_eunit = "Energy units to use in the y-axis. Defaults to eV, which will autoscale QP weights if -qp is set."
    help_norbitals = "Number of orbital contributions to plot for QP weights. Defaults to 3 (s,p,d) and must be at least 1."
    #---------------------------------------------------------------------------

    p.add_argument('-d','--directory',
                   nargs = '*', default = ["."],
                   type = str, help = help_directory)

    p.add_argument('-f','--files',
                   nargs = '*', default = [],
                   type = str, help = help_files)

    p.add_argument('-y1','--ymin', default = [None, None],
                   type = float, help = help_ymin)

    p.add_argument('-y2','--ymax', default = [None, None],
                   type = float, help = help_ymax)

    p.add_argument('-s','--scale', default = [1.0], nargs="*",
                   type = float, help = help_scale)

    p.add_argument('-t','--title',
                   type = str, default = None, help = help_title)

    p.add_argument('-nt','--no_title', action='store_true', help = help_no_title)

    #p.add_argument('-nr','--nrows', type=int, help = help_nrows)
    #p.add_argument('-nc','--ncols', type=int,  help = help_ncols)
    p.add_argument('-p', '--path', type=str, help=help_path, default='hexagonal',choices=['hexagonal', 'cubic'])

    p.add_argument('-g','--grid', action='store_true', help = help_grid)
    p.add_argument('-shx','--share_x', default="col", help = help_sharex)
    p.add_argument('-shy','--share_y', default="col", help = help_sharey)
    p.add_argument('-eu', '--eunit', default="eV", choices=['eV', 'Ha'], help=help_eunit)
    p.add_argument('-no', '--norbitals', type=int, default=3, help=help_norbitals)
    #---------------------------------------------------------------------------

    args = p.parse_args()
    input_options = {}

    input_options['directory'] = args.directory

    if ( len(args.files)==0 ):
        sys.exit("\n ERROR: At least a filename must be specified:\n\n"
                 +"        PLOT-exciton-weights.py -f [FILES [FILES ...]]\n")
    input_options['files'] = args.files

    input_options['ymin'] = args.ymin
    #if ( len(args.yboundary) >= 1 ): input_options['ymin'] = args.yboundary[0]
    input_options['ymax'] = args.ymax
    #if ( len(args.yboundary) >= 2 ): input_options['ymax'] = args.yboundary[1]

    input_options['scale'] = args.scale
    input_options['title'] = args.title
    input_options['no_title'] = args.no_title
    input_options['grid'] = args.grid
    #input_options['nrows'] = args.nrows
    #input_options['ncols'] = args.ncols
    input_options['share_x'] = args.share_x
    input_options['share_y'] = args.share_y
    input_options['path'] = args.path
    input_options['eunit'] = args.eunit
    input_options['norbitals'] = args.norbitals
    return input_options

# Remove or add image output formats
output_formats: list = ["png", "eps", "pdf"]
def read_input(filename: str):
        """
        Reads KPATH and BANDLINE data."""
        bandwgt, bandwgtdim = readrawdata(filename)
        if os.path.isfile('../BANDLINES.OUT'):
            bandlin, bandlindim = readrawdata( '../BANDLINES.OUT')
        else:
            bandlin, bandlindim = readrawdata( './BANDLINES.OUT')
        bandlines = bandlin[:,0,0]

        return bandwgt, bandwgtdim, bandlin, bandlindim, bandlines
def read_evalcore_species() -> dict:
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
            #information = line.split(" ").remove(":")
            #if any(species in line for species in all_species):
            species_index[re.search(r'\d+',line).group()] =line[line.find("(")+1:line.find(")")]
        return species_index

    except OSError:
        print('Cannot open EVALCORE.OUT.')


def determine_species(filename: str) -> str:
    """
    Determines the species of BAND_Sss_Aaaaa.OUT files for contribution
    """
    file_name = filename.strip(".OUT").split("_")
    species = read_evalcore_species()
    species_index = file_name[1].strip("S")
    # Indices in evalcore have no leading zero, convert to int then str 
    return species[str(int(species_index))]

def get_orbital_contribution(infile: list[str]):
    """
    Input:
        infile: list of strings of filenames with QP weights of the same species.

    Description:
    - Determine orbital contribution from given files of same species.
    - Given BAND_Sss_Aaaaa.OUT, 4th 5th, 6th columns determine QP contributions from s,p,d orbitals
    - Possible for larger atoms that there are more columns - this will handle it.
    - This assumes the first 2 columns of the file are KPATH,ENERGY and the rest are ORBITAL WEIGHT.

    Returns:
    Tuple of np array containing total QP contribution from each orbital of a given species.
    """
    orbital_contributions = []
    for ifile in infile:
        bandwgt = read_input(ifile)[0]
        # columns 3,4,5,6 etc => indices 2 3 4 5 etc
        # need to determine number of orbitals (columns above 3)
        # first 3 columns always kpath, energy, total weight so ignore columns 0,1 2
        # Will be array of N_b x (c - 3) x N
        orbital_data = bandwgt[:,3:] * ha2ev
        orbital_contributions.append(orbital_data)
    s_contrib = np.zeros((1,len(orbital_contributions[0]), len(orbital_contributions[0][0][0])))
    p_contrib = np.zeros((1,len(orbital_contributions[0]),  len(orbital_contributions[0][0][0])))
    d_contrib = np.zeros((1,len(orbital_contributions[0]),  len(orbital_contributions[0][0][0])))
    for contrib in orbital_contributions:
        for index,block in enumerate(contrib):
            s_contrib[0][index] += block[0]
            p_contrib[0][index] += block[1]
            d_contrib[0][index] += block[2]
    return [s_contrib, p_contrib, d_contrib]
        
def l2norm(x,y):
    denom = np.max(y) - np.min(y)
    if denom == 0:
        return y
    else:
        normed = (y - np.min(y)) / (denom)
        return normed
def main(input_options):
    '''
    input:
    :input_options: dictionary that holds the input options parsed 
                    from the command line arguments
    '''
    directory = input_options['directory']
    files: list[str] = input_options['files']
    ymin: float = input_options['ymin'] 
    ymax: float = input_options['ymax'] 
    scale: list[float] = input_options['scale'] 
    title: str = input_options['title']
    no_title: bool = input_options['no_title']
    grid: bool = input_options['grid']
    #ncols: int = input_options['ncols']
    #nrows: int = input_options['nrows']
    share_x: bool = input_options['share_x']
    share_y: bool = input_options['share_y']
    path: tuple = input_options['path']
    eunit: str = input_options['eunit']
    norbitals: int = input_options['norbitals']
    path_BZ = ('A','H','$\Gamma$','M','L','K') if path == 'hexagonal' else ('W','L','$\Gamma$','X','W','K')
    
    
    #########################
    # Settings for the plot #
    #########################

    figcolor = 'white'

    mpl.rcParams['axes.linewidth'] = 1.0 # set the value globally
    mpl.rcParams['grid.linewidth'] = 1.0
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.labelsize'] = 10    # fontsize of the x and y labels
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['axes.axisbelow'] = 'True'   # whether axis gridlines and ticks are below
                                            # the axes elements (lines, text, etc)
    mpl.rcParams['legend.fontsize'] = 25
    plt.rcParams['xtick.major.pad'] = 5
    plt.rcParams['ytick.major.pad'] = 5

    plt.rcParams.update({'mathtext.default':'regular'})
    subshells = ["s", "p", "d", "f", "g", "h"]
    colors = ["c", "m", "y"]
    # Extract Data and Plot
    counter = 0
    # Sort QP files by species
    sorted_files = []
    read_species = read_evalcore_species()
    for i in range(0,len(read_species)):
        sorted_files.append([])
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
    qp_data = []
    print(sorted_files)
    if len(scale) != 1 and len(scale) != len(sorted_files) * norbitals:
        print("Numebr of scale values provided is incorrect. Either provide a single float or len(files) * norbitals float values")
        sys.exit(1)
    ordered_species = []
    for files in sorted_files:
        if len(files) == 0:
            continue
        # Get correct order of species to label plots
        ordered_species.append(determine_species(files[0]))
    
    for atomfiles in sorted_files:
        data = get_orbital_contribution(atomfiles)
        qp_data.append(data)
    print("len", len(qp_data[0]))
    

    if (len(files) > norbitals * len(qp_data)):
        print("WARNING: Number of plots greater than available subplot space. Please change the number of files you are plotting, or norbitals.")
        sys.exit(1)
    # Automatically make subplot grid based on species and orbital contributions rather than ncols,nrows
    fig, axs = plt.subplots(len(sorted_files), norbitals, sharex=share_x, sharey=share_y, figsize=(10,10))
    fig.patch.set_edgecolor(figcolor)
    fig.patch.set_facecolor(figcolor)
    species_index = 0
    orbital_index = 0
    counter2 = 0
    flat = []
    for xs in sorted_files:
        for x in xs:
            flat.append(x)

    print("flat", flat)
    # Read QP bandstructure and plot it instead of KS BS
    for row in axs:
        row[0].set_ylabel( 'Energy [eV]')
    for i, ax in enumerate(fig.axes):
        if i % 3 == 0 and i != 0: 
            species_index += 1
            orbital_index = 0
            counter2 = 0
        if orbital_index >= norbitals and i >= len(axs.flatten()): break
        if species_index >= len(ordered_species): break
        col2 = ax.twinx()
        output = read_input(flat[0])
        # if ( eunit == "eV"):
        #     qp_data[species_index] *= 1
 
        ax.xaxis.grid( True, which='major', color='k', linestyle='-', linewidth=1)
        ax.xaxis.set_label_position('bottom')
        ax.set_xticks(output[4])
        ax.set_xticklabels(path_BZ)
        for line in ax.get_xticklines() + ax.get_yticklines():
            line.set_markersize(5)
            line.set_markeredgewidth(1)
        xmin = np.amin( output[0][:,0,:])
        xmax = np.amax( output[0][:,0,:])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        print(colors[species_index])
        col2.text(0.2,0.2,f"{ordered_species[species_index]} {subshells[orbital_index]}",zorder=50, backgroundcolor="w", color=colors[species_index], fontsize=12, bbox=dict(boxstyle= "square", edgecolor=colors[species_index], facecolor="w"))
        col2.text(0,6,f"{scale[counter]}",zorder=50, backgroundcolor="w", color=colors[species_index], fontsize=12, bbox=dict(boxstyle= "square", edgecolor =colors[species_index], facecolor="w"))
        for j in range( output[1][2]):
            ax.plot( output[0][j,0,:], output[0][j,1,:]*ha2ev, 'k', lw=1.0, zorder=10)
        if len(scale) > 1:
            for j in range( output[1][2]):
                #normedcircles = l2norm(np.arange(0,1,1/len(qp_data[species_index][j,orbital_index,:])),qp_data[species_index][j,orbital_index,:])
                ax.scatter( output[0][j,0,:], output[0][j,1,:]*ha2ev, s=(scale[counter]*qp_data[species_index][orbital_index][0][j])**2, lw=0.35, edgecolor=colors[species_index], facecolor='none', zorder=11)
        else:
            for j in range( output[1][2]):
                ax.scatter( output[0][j,0,:], output[0][j,1,:]*ha2ev, s=(scale[0]*qp_data[species_index][orbital_index][0][j])**2, lw=0.35, edgecolor=colors[species_index], facecolor='none', zorder=11)
        
    
        
                #Fermi level
        ax.plot( [xmin, xmax], [0, 0], 'k', lw=2.0, ls='-')
            
        col2.set_ylim( ax.get_ylim())
        col2.set_yticks( [0, 0])
        col2.set_yticklabels( ('$\\mathregular{E_{F}}$', ''))

        ax.grid( True)
        plt.title( title, fontsize=mpl.rcParams['ytick.labelsize'], y=1.03)

        counter2 += 1
        counter +=1
        orbital_index += 1

    #save file format
    for format in output_formats:
        if format == "png":
            fig.savefig(f"PLOT.{format}", format="png", dpi=600, bbox_inches='tight')
        else:
            fig.savefig(f"PLOT.{format}",format=format, bbox_inches="tight")
    sys.exit()

if __name__ == "__main__":
    input_options = option_parser()
    main(input_options)
