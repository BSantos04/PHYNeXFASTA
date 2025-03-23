import argparse
import sys

def fasta_to_dict(file):
    """
    Summary:
        This function converts a FASTA file into a dictionary containing the sequences and sequence names.
    
    Parameters:
        file: Input FASTA file.
    
    Returns:
        fasta_dict: Dictionary with the sequence names and respective sequences.
    """
    try:
        # Create a dictionary object to store the sequences and respective names
        fasta_dict = {}
        
        # Open input file
        with open(file) as infile:
            # Remove every empty space on each line of the file
            for line in infile:
                line = line.strip()
                # If the line is "empty", skip that line
                if not line:
                    continue
                # If it starts with a header, stores in the dictionary as a key (sequence ID)
                if line.startswith(">"):
                    header = line[1:]
                    fasta_dict[header] = ""
                # All the other lines are stored as values (sequences)
                fasta_dict[header] = line
        return fasta_dict
    except Exception as e:
        # If something goes wrong, it will display a message and terminate the process
        print(f"An error has ocurred: {e}")
        sys.exit(1)
    

def phylip_to_dict(file):
    """
    Summary:
        This function converts a PHYLIP file into a dictionary containing the sequences and sequence names.
    
    Parameters:
        file: Input PHYLIP file.
    
    Returns:
        dict(phylip_matrix): Dictionary with the sequence names and respective sequences.
    """
    try:
        # Create two empty lists
        phylip_list = []
        phylip_matrix = []
        
        # Open input file
        with open(file) as infile:
            for line in infile:
                # Remove every empty space on each line of the file
                line = line.strip()
                # If the line is "empty", skip that line
                if not line:
                    continue
                # Add every line to the first list
                phylip_list.append(line)
            # Remove the first item of the list (the number and length of sequences)
            phylip_list = phylip_list[1:]
            # Create a matrix with sequences and respective IDs
            for item in phylip_list:
                phylip_matrix.append(item.split(maxsplit=1))
        # Convert the matrix into a dictionary
        return dict(phylip_matrix)
    except Exception as e:
        # If something goes wrong, it will display a message and terminate the process
        print(f"An error has ocurred: {e}")
        sys.exit(1)

def nexus_to_dict(file):
    """
    Summary:
        This function converts a NeXus file into a dictionary containing the sequences and sequence names.
    
    Parameters:
        file: Input NeXus file.
    
    Returns:
        dict(nexus_matrix): Dictionary with the sequence names and respective sequences.
    """
    try:
        # Create two empty lists
        nexus_list = []
        nexus_matrix = []
        
        # Open input file
        with open(file) as infile:
            for line in infile:
                # Remove every empty space on each line of the file
                line = line.strip().upper()
                # If the line is "empty", skip that line
                if not line:
                    continue
                # Skip every non, sequence line
                if line.startswith(("#NEXUS", "DIMENSIONS", "FORMAT", "MATRIX", "BEGIN DATA")):
                    continue
                # Stop adding lines to the list if it encounters a ";"
                if line == ";":
                    break
                # Add sequence lines to the list
                nexus_list.append(line)
            # Create a matrix with every sequence and their respective IDs
            for item in nexus_list:
                nexus_matrix.append(item.split(maxsplit=1))
        # Convert the matrix into a dictionary
        return dict(nexus_matrix)
    except Exception as e:
        # If something goes wrong, it will display a message and terminate the process
        print(f"An error has ocurred: {e}")
        sys.exit(1)

def dict_to_fasta(file_dict, outname):
    """
    Summary:
        This function converts a NeXus file into a dictionary containing the sequences and sequence names.
    
    Parameters:
        file: Input NeXus file.
    """
    try:
        # Create te output file based on the given dictionary and output name
        with open(f"{outname}.fasta", "w") as outfile:
            # Write the sequences into the output file in FASTA format
            for seqid, seq in file_dict.items():
                sequence = f">{seqid.split()[0]}\n{seq}\n"
                outfile.write(sequence)
    except Exception as e:
        # If something goes wrong, it will display a message and terminate the process
        print(f"An error has ocurred: {e}")
        sys.exit(1)

def dict_to_phylip(file_dict, outname):
    """
    Summary:
        This function converts a dictionary with sequences and respective sequence names into a PHYLIP file.
    
    Parameters:
        file_dict: Dictionary with sequences.
    """
    try:
        # Create te output file based on the given dictionary and output name
        with open(f"{outname}.phy", "w") as outfile:
            # Write the first line of the file with the number of sequences and their length
            values = list(file_dict.values())
            first_line = f"{len(file_dict)} {len(values[0])}\n"
            outfile.write(first_line)
            # Write the sequences into the output file in PHYLIP format
            for seqid, seq in file_dict.items():
                sequence = f"{seqid.split()[0]}   {seq}\n"
                outfile.write(sequence)
    except Exception as e:
        # If something goes wrong, it will display a message and terminate the process
        print(f"An error has ocurred: {e}")
        sys.exit(1)

def dict_to_nexus(file_dict, outname, outgroup = None):
    """
    Summary:
        This function converts a NeXus file into a dictionary containing the sequences and sequence names.
    
    Parameters:
        file: Input NeXus file.
    """
    try:
        # Create te output file based on the given dictionary and output name
        with open(f"{outname}.nex", "w") as outfile:
            # Write NeXus file header
            header = f"#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX={len(file_dict)} NCHAR={len(next(iter(file_dict.values())))};\nFORMAT DATATYPE=DNA MISSING=N GAP=-;\nMATRIX\n"
            outfile.write(header)
            # Write the sequences into the output file in PHYLIP format
            for seqid, seq in file_dict.items():
                sequence = f"    {seqid.split()[0]}  {seq}\n"
                outfile.write(sequence)
            # Write MrBayes block based on given outgroup
            end = "  ;\nEND;\n"
            mrbayes = f"\nbegin mrbayes;\n  set autoclose=yes;\n  outgroup {outgroup};\n  mcmcp ngen=200000 printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename={outname};\n  mcmc;\n  sumt filename={outname}\nend;\n"
            outfile.write(end+mrbayes)
    except Exception as e:
        # If something goes wrong, it will display a message and terminate the process
        print(f"An error has ocurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Define input arguments
    parser = argparse.ArgumentParser(description="Converts an input file within these three formats: FASTA, NeXus and PHYLIP.")
    parser.add_argument("--infile", type=str, help="Input file (.fasta, .fa, .nex or .phy)")
    parser.add_argument("--format", type=str, help="Output file format (FASTA, NeXus or PHYLIP.")
    parser.add_argument("--outname", type=str, help="Name of the output file (same as the input if not specified)", required=False)
    parser.add_argument("--outgroup", type=str, help="Name of the outgroup", required=False)
    args = parser.parse_args()
    
    # If the input format isn't one of the options in the array, it will raise an error message and exit
    if args.format.upper() not in {"NEXUS", "FASTA", "PHYLIP", "NEX", "FA", "PHY"}:
        parser.error("File format not recognized!\nTry FASTA, NeXus or PHYLIP!")
    
    # Create an empty list to store every line of the input file
    identifier = []
    
    # Add every line of the input file
    with open(args.infile) as infile:
        for line in infile:
            line = line.strip()
            identifier.append(line)
    
    # Recognize the format of the input file based on it's first line
    if identifier[0].startswith(">"):
        file_dict = fasta_to_dict(args.infile)
    elif identifier[0].strip() == "#NEXUS":
        file_dict = nexus_to_dict(args.infile)
    elif identifier[0].isdigit():
        file_dict = phylip_to_dict(args.infile)
    else:
        # Exit the program if it doesn't recognize the file format
        print("Input file format not recognized!!!")
        sys.exit(1)
    
    # Ensure that the input file is aligned based on it's sequences length, otherwise it will dispay an error message
    same_len = all(len(i) == len(next(iter(file_dict.values()))) for i in file_dict.values())
    if same_len == True:
        # Based on the selected output format, it will convert the dictionary into a file of the specified format
        # If a filename is specified, it will use that filename
        if args.outname:
            if args.format.upper() == "FASTA" or args.format.upper() == "FA":
                dict_to_fasta(file_dict, outname=args.outname)
            if args.format.upper() == "NEXUS" or args.format.upper() == "NEX":
                dict_to_nexus(file_dict, outname=args.outname, outgroup = args.outgroup)
            if args.format.upper() == "PHYLIP" or args.format.upper() == "PHY":
                dict_to_phylip(file_dict, outname=args.outname)
        else:
            filename = args.infile.split(".")[0]
            if args.format.upper() == "FASTA" or args.format.upper() == "FA":
                dict_to_fasta(file_dict, outname=filename)
            if args.format.upper() == "NEXUS" or args.format.upper() == "NEX":
                dict_to_nexus(file_dict, outname=filename, outgroup = args.outgroup)
            if args.format.upper() == "PHYLIP" or args.format.upper() == "PHY":
                dict_to_phylip(file_dict, outname=filename)
    else:
        print("The sequences must be aligned!!!")
        sys.exit(1)
