# strling-denovo
Uses pedigree information to filter STRling output to find de novo variants
Will need a ped file as well as STRling output. For peddy library, must be in python version 3.7 or earlier.

#strling-denovo-threshold150
The same as above, except in the STRling output the threshold for being classified as a true novel amplification is a difference greater than 150 compared to both parents.

This repository also holds various pieces of code I have used for data analysis and visualization.

## Install
`conda env create environment.yml`
`conda activate strling-denovo`
