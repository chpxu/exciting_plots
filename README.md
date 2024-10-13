# exciting_plots
Repository hosting modified plotting files I have made for data produced from exciting. The following sections are documentation for each file. 

The prescence of `readrawdata.py` is a necessity for PDOS and PBS to ensure maximum compatibility with previous versions of exciting.

Requires lxml and matplotlib. `EVALCORE.OUT` should be present in the directory where you are trying to run these files.

## PDOS
To plot PDOS, please add the attribute `lmirep = "true"` to the `<dos>` element and run `exciting`. Then, for every `<atom />` declared inside a `<species>` element, you will find files of the format `PDOS_Sxx_Ayyyy.OUT` where `Sxx` is the species number (indexed starting from `xx=01`, in the order that they were declared inside `input.xml`) and `Ayyyy` is the atom number (indexed starting from `yyyy=0001` in the order they were declared for its associated `<species>`). These will all be outputted in the directory where `exciting` was ran, typically in the same directory as `input.xml`.

### Interpreting the PDOS OUT files
If you look inside `$EXCITINGROOT/species/*.xml` files, you will see all the different electronic states for your atoms: the 1s, 2s, 2p_x, 2p_y etc. In general, there are one s states, 3-p states, 5 d states, 7 f states and so on (i.e., all the s-states are considered together and so on).

The PDOS files are split into blocks of size `nwdos` (declared on the `<dos>` element) where the number of blocks equals the sum of the number of orbitals. If for example, the available states of any of your atoms goes up to and includes f states, and you have `nwdos=1000`, then you will have `(1 + 3 + 5 + 7)` blocks of `1000` data points. Every PDOS file has 2 columns, the first column is the energy in the range $[e_1, e_2]$ (declared by `winddos="e_1 e_2"`) and the 2nd column is the PDOS at that point in states/Ha/unit cell. Therefore:
- The total s-contribution from an atom of a species is given by the first block
- The total p-contribution from an atom of a species is given by the sum of the 2nd,3rd,4th blocks
- The total d-contribution from an atom of a species is given by the sum of the 5th-9th blocks inclusive.
- And so on

For multiple atoms of the same species, the blocks are then added together to give the total orbital contribution from a species.

### Using `pdos.py`
The file to call is `pdos.py` and should work for Python `>= 3.7`, NumPy `< 2`. Not tested with NumPy `2.x`. Example usage:
```sh
python3 pdos.py -f [files] -e -2 6 -db 0 2.5
```

- `-f` is the list of PDOS files. Order does not matter - it sorts them In the future, can just make a function to automatically find PDOS files.
- `-e` is the energy range, i.e. x-axis range
- `-db` is the y-axis (DOS) range. Note that the data is automatically converted to states/eV/unit cell


## Projected Bandstructure (PBS)
To plot PBS, please add the attribute `character="true"` to the `<bandstructure>` element and run `exciting`. Then, for every `<atom />` declared inside a `<species>` element, you will find files of the format `BAND_Sxx_Ayyyy.OUT` where `Sxx` is the species number (indexed starting from `xx=01`, in the order that they were declared inside `input.xml`) and `Ayyyy` is the atom number (indexed starting from `yyyy=0001` in the order they were declared for its associated `<species>`). These will all be outputted in the directory where `exciting` was ran, typically in the same directory as `input.xml`.

### Interpreting the BAND files
If you look inside `$EXCITINGROOT/species/*.xml` files, you will see all the different electronic states for your atoms: the 1s, 2s, 2p\_x, 2p\_y etc. In general, there are one s states, 3-p states, 5 d states, 7 f states and so on (i.e., all the s-states are considered together and so on).

The BAND files are split into blocks of size `steps` (declared on the `<path>` element) where the number of blocks equals the number of bands. Every BAND file has $3 + N\_s$  columns ($N\_s$ is the number of species), the first column is the KPATH, the 2nd column is the energy in the range $[e_1, e_2]$  and the 3rd column is the sum of the contributions for each k-point. Therefore:
- The total s-contribution at a k-point from an atom of a species is given by the 4th column
- The total p-contribution at a k-point from an atom of a species is given by the sum of the 5th column
- The total d-contribution at a k-point from an atom of a species is given by the sum of the 6th column.
- And so on

For multiple atoms of the same species, the blocks are then added together to give the total orbital contribution from a species.

### Using `qp_weights.py`
The file to call is `qp_weights.py` and should work for Python `>= 3.7`, NumPy `< 2`. Not tested with NumPy `2.x`. Example usage:
```sh
python3 qp_weights.py -f [files] -y1 e_1 -y2 e_2 -s float | [float] -p "cubic" | "hexagonal" -shx "none" | "all" -shy "none" | "all" --norbitals 3 --scissor 0.0 | float
```

- `-f` is the list of BAND files. Order does not matter - it sorts them In the future, can just make a function to automatically find BAND files.
- `-y1 -y2` is the energy range, i.e. y-axis range. Energy is automatically converted to eV.
- `-p "cubic" | "hexagonal"` controls the x-axis tick labels for the path. Can be manually adjusted to your liking.
- `-norbitals` is the number of orbitals to plot. By default, it is 3 (s,p,d).
- `-s`  is the scale factor. A single float scales the weights for all plots. Otherwise must be a list of floats of size $N\_s \times \text{norbitals}$
- `-shx` is "none" or "all" and controls whether each plot has its own labelled x-axis ("none") or shared.
- `-shy` is "none" or "all" and controls whether each plot has its own labelled y-axis ("none") or shared.
- `--scissor` is the scissor shift to apply to the conduction bands in eV.

This file will output bandstructures in a $N_s \times \text{norbitals}$ grid with a box in the top left containing the scale factor, and each plot labelled with species and orbital.

