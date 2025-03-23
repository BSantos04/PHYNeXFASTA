# PHYNeXFASTA
This project consists solely as a challenge to create a Python script that converts files within three different formats (FASTA, NeXus and PHYLIP) without any external package/module.

# Requisities 
-Python3 3.12.2 (but I guess any Python 3 version can run it without any trouble)

# Instalation
`git clone https://github.com/BSantos04/PHYNeXFASTA.git` 

# Usage
## Basic usage
`python3 PHYNeXFASTA.py --infile {path/to/input/file} --format {fasta; nexus; phylip; fa; nex; phy}`

### Example 
`python3 PHYNeXFASTA.py --infile example.nex --format fasta`

## Other flags
`--outname:` Specifies the output filename. If not specified, the script will use the same filename as the input file by default.

`--outgroup:` A useful option when converting to the NeXus format. Since the script also creates a MrBayes block with some default features, you can also specify the outgroup species in it using that flag.

### A more detailed example
`python3 PHYNeXFASTA.py --infile example.fasta --format nex --outname test --outgroup Podarcis`

# License 
GNU General Public License V3.0
