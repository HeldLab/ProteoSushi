"""ps_utilities.py: helper functions used by parsers in proteoSushi"""

__author__ = "Arshag Mooradian, Rob Seymour"
__email__ = "rseymour@wustl.edu"

import pickle
import csv
from collections import defaultdict
import os
from ps_digest import digest, fasta_producer, cleave_rule_determination


def parse_maxquant_summary(infile: str) -> list:
    """Pulls out information needed from the MaxQuant files
    Arguments:
        infile {str} -- the name of the MaxQuant summary file that will be used
    Returns:
        str -- the name of the enzyme used (usually trypsin)
        int -- the maximum allowed number of missed cleavages
    """
    mqfile = open(infile, 'r')
    #outfile = open("mq_parsed.csv", 'w', newline='')
    #assert infile.endswith(".txt")
    tsvReader = csv.reader(mqfile, delimiter='\t', quotechar='"')
    #csvWriter = csv.writer(outfile)
    enzyme = ""
    PTMs = []
    maxMissedCleaves = 0
    # From proteinGroups: Need Intensity[DM], PTM data[GS-GX], Potential Contaminant(?)[GK]
    # From summary: Enzyme[C], Variable Modifications[H][I-J?], Max. missed cleavages[O]
    # From parameters: Modifications included in protein quantification,
    header = next(tsvReader)  # Skips the header line
    row = next(tsvReader)
    enzyme = row[header.index("Enzyme")]
    assert enzyme != ""
    #PTMs = row[header.index("Variable modifications")].split(';')  # I might just do this in the evidence file
    maxMissedCleaves = int(row[header.index("Max. missed cleavages")])
    mqfile.close()
    return enzyme.lower(), maxMissedCleaves

def parse_mascot(in_file: str) -> list:
    """Pulls out enzyme, quantitative range, and variable modifications, along
    with the quant data from a Mascot file.
    Arguments:
        in_file {str} -- a string for the name of the mascot file
    Returns:
        list -- enzyme {str}, quant_range {tuple}, var_mod_map {dict},
        missed_cleaves {int}
    """
# We need to parse Mascot file first, because raw Mascot output is difficult to work with
    print(f"Parsing file {in_file}.")
    with open(in_file) as f, \
        open('quant_parsed.csv', 'w', newline='') as out_file:
        # Quantitative section of Mascot file
        assert in_file.endswith('.csv'), "File must be a csv"
        file_reader = csv.reader(f)
        out_writer = csv.writer(out_file)
        enzyme = ""
        missed_cleaves = -1
        for row in file_reader:
            if row:  # lots of blank rows
                # Pull data we need
                if row[0] == "Enzyme":
                    enzyme = row[1].lower()
                elif row[0] == "Variable modifications" and '--' not in row[1]:
                    var_mods = row[1].split(',')
                    var_mod_map = {val.lower():
                                   str(pos+1) for pos, val in enumerate(var_mods)}
                elif row[0] == "Maximum Missed Cleavages":
                    missed_cleaves = int(row[1])
                elif row[0] == "prot_hit_num" or row[0] == "prot_acc":
                    # This should be last elif; start of quantitative info
                    header = row
                    second_row = next(file_reader)  # need to normalize headers
                    quant_range = (len(header), len(second_row)-1)  #range where there will be quant; appears to skip some
                    header.extend(['SampleQuant' 
                                   for i in range(len(second_row) - len(header))])
                    out_writer.writerows([header, second_row])
                    for row in file_reader:
                        out_writer.writerow(row)
                    break  # Exit loop as we have all of the header information
        print('done parsing')
        return enzyme, quant_range, var_mod_map, missed_cleaves


def load_pepdict(proteome_fasta_filepath: str, enzyme: str, missed_cleaves: int) -> dict:
    """Looks for pep_dict; if it doesn't find, it looks for a FASTA to make one.
    Fasta would be in current directory but this can be changed.

    Arguments:
        proteome_fasta_filepath {str} -- the filepath for the proteome fasta file
        enzyme {str} -- protease used in the sample digestion
        missed_cleaves {int} -- the max number of allowed missed cleavages
    Returns:
        dict -- a dictionary that 
    """
    pep_dict = None
    wd = os.getcwd()
    pep_dict_file = f"{proteome_fasta_filepath.split('.')[0]}.pepdict"
    if not any([x == pep_dict_file for x in os.listdir(wd)]):
        print("Pepdict not found in current directory. Trying to generate from FASTA")
        digest_if_needed(proteome_fasta_filepath, enzyme, missed_cleaves)
    '''
    pep_dicts = 0
    for f in os.listdir(wd):
        if f.endswith('.pepdict'):
            pep_dicts += 1
            pep_dict_file = f
    assert pep_dicts == 1, \
        f"Found {pep_dicts} pepdict files. Need exactly 1 in {wd}."
        '''
    
    print(f"Loading {pep_dict_file} into memory")
    with open(pep_dict_file, 'rb') as f:
        pep_dict = pickle.load(f)
    print("Pept Dictionary Loaded.")
    return pep_dict


def digest_if_needed(proteome_fasta_filepath: str, enzyme: str, mcleave: bool, minlen=5, maxlen=55, 
                     m_excise="Both"):
    """Runs digest and dumps pickle file (dict)"""
    seq_list = []  # (gene, organism, sequence, unpid)
    seq_list = fasta_producer(proteome_fasta_filepath, seq_list)
    rule = cleave_rule_determination(enzyme)
    pep_dict = defaultdict(set)  # NOTE: this won't work for large fasta/sequence databases as it's in ram; ~1GB for human
    for i in seq_list:
        gene, _organism, seq, unpid = i
        if m_excise == "Both":
            # Run digest results twice and combine; inefficient but gets job done
            digest_results_true = digest(seq, minlen, maxlen, mcleave, True,
                                         li_swap=True, rule_to_use=rule,
                                         reverse=False)
            digest_results_false = digest(seq, minlen, maxlen, mcleave, False,
                                          li_swap=True, rule_to_use=rule,
                                          reverse=False)
            digest_results = set().union(digest_results_true,
                                         digest_results_false)
        else:
            digest_results = digest(seq, minlen, maxlen, mcleave, mexcise,
                                    li_swap=True, rule_to_use=rule,
                                    reverse=False)
        for result in digest_results:
            peptide, position = result
            peptide = peptide.replace("L", "I")
            pep_dict[peptide].add((gene, position, unpid))  # careful with site specificity here
    with open(f"{proteome_fasta_filepath.split('.')[0]}.pepdict", "wb") as w1:
        pickle.dump(pep_dict, w1)
    del pep_dict  # Free memory for later
