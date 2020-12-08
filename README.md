# Si-H_GAP_training
Code and scripts used to iteratively train the Si-H GAP potential. Full version is hosted on Dropbox, to more effectively host the large numbers of data files. This repo contains some sample data, but a minimal number in order to reduce push/ pull times. 

Most relevant files and folders:
- Espresso-to-ExtendedXYZ.py, the code used to take QE DFT outputs and convert them into extended xyz format, the format used by QUIP and GAP
- ExtendedXYZ-to-Espresso.py, the code used to convert extended xyz format files into DFT inputs. Used at the beginning of the project to construct the training database, and now used to run DFT on trial MD simulations, and feed poor fits back into the GAP training database
- "Pair Potentials", the folder which contains the method by which we created 2-body potentials to add on to GAP (to help it fit the high slope potential regions more efficiently)
- Iterative_Training/, the folder which contains the scripts used to create new files to run trial MD simulations on. Typical simulations are optimization, annealing, and heating/ quenching
- "Iterative_Training_DFT/, the folder containing sample DFT inputs and outputs from this iterative training process. Most data has been scrubbed to preserve a small repo size
- GAP_fitting/, the folder containing all scripts and files for the actual training procedure. It also includes scripts for visualizing and comparing GAP results to DFT. Most data has been scrubbed apart from a few examples and images, to preserve a small repo size
