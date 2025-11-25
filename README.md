# PHYNeXFASTA
This project was a challenge that involved creating a Python script to convert files between three different formats (FASTA, NeXus, and PHYLIP) without using any external packages or modules.

# Requisities 
-Python >= v3.6 (tested with Python v3.12.2)

# Instalation
`git clone https://github.com/BSantos04/PHYNeXFASTA.git` 

# Usage
## Basic usage
`python3 PHYNeXFASTA.py --infile {path/to/input/file} --format {fasta; nexus; phylip; fa; nex; phy}`

### Example 
`python3 PHYNeXFASTA.py --infile example.nex --format fasta`

## Other flags
`--outname:` Specifies the output filename. If not specified, the script will use the same filename as the input file by default.

`--outgroup:` An useful option when converting to the NeXus format. Since the script also creates a MrBayes block with some default features, you can also specify the outgroup species in it using that flag.

### More detailed example
`python3 PHYNeXFASTA.py --infile example.fasta --format nex --outname test --outgroup Podarcis`

# License 
GNU General Public License V3.0
