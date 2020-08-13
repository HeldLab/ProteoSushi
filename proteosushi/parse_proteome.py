"""parse_proteome.py: parses the proteome FASTA file to make sure it is in an acceptable format
and grabs the species ID"""

import os

def fasta_entry_check(header: str, sequence: str) -> str:
    """Checks whether the header has organism info and returns as tuple

    Arguments:
        header {str} -- header line
        sequence {str} -- protein sequence
    Returns:
        str -- species ID
    """
    unpid = "N/A"
    # Check organism
    if 'OX=' in header:
        organism = int(header.split('OX=')[1].split()[0])
    elif "OS=" in header:
        organism = header.split('OS=')[1].split(" GN=")[0]
    else:
        return None  # skip altogether if no organism info in line
    # Attempts to grab the Uniprot ID from the header
    if header.count('|') >= 2:
        unpid = header.split('|')[1]

    if unpid == "N/A":
        return "ERROR"

    return str(organism)

def parse_proteome(filename: str) -> list:
    """checks FASTA file to make sure it is in the correct format and grabs the species ID

    Arguments:
        filename {str} -- the filepath+name of the FASTA file
    Returns:
        bool -- whether an error occurred during parsing (e.g. incorrect format)
        str -- the species ID from the fasta file
    """
    species_IDs = list()
    if os.path.exists(filename):
        with open(filename) as input_file:
            print(f"Checking proteome file: {filename}")
            header = None
            seq = ''
            for line in input_file:
                if line.startswith('>') and seq:
                    species_ID = fasta_entry_check(header, seq)
                    if species_ID == "ERROR":
                        return True, None
                    elif species_ID is None:
                        pass
                    else:
                        species_IDs.append(species_ID)
                    header = line.strip()
                    seq = ''
                elif line.startswith('>'):  # first entry of file
                    header = line.strip()
                else:
                    seq += line.strip()
            fasta_entry_check(header, seq)  # For last entry
            if species_ID == "ERROR":
                return True, None
            elif species_ID is None:
                pass
            else:
                species_IDs.append(species_ID)

    # If there are no species listed, it will send back None
    if not species_IDs:
        return False, None
    # Grabs the most common of the species IDs and uses that
    return False, max(set(species_IDs), key=species_IDs.count)
#EOF