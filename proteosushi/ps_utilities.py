"""ps_utilities.py: helper functions used by parsers in proteoSushi"""

__author__ = "Arshag Mooradian, Rob Seymour"
__email__ = "rseymour@wustl.edu"

import pickle
import csv
from collections import defaultdict
import ntpath
import os
from re import finditer

from .ps_digest import digest, fasta_producer, cleave_rule_determination


def clean_pep_seq(rule: tuple, pep_mod_seq: str, user_PTMs: list, old_pep_seq: str) -> str:
    """Removes the missed cleavages without selected PTMs

    Arguments:
        rule {tuple} -- contains info about the protease
        pep_mod_seq {str} -- the peptide sequence with modifications
        user_PTMs {list} -- a list of PTMs that the user has selected for Analysis
    Returns:
        str -- the cleaned peptide
    """
    #print(user_PTMs)
    # Puts break points for unmodified sequence in a list
    breaks = finditer(rule[0], old_pep_seq)
    cut_unmod_sites = []
    cut_unmod_peptides = []
    for breakp in breaks:
        cut_unmod_sites.append(breakp.start())
    # Adds the unmodified cut peptides in one at a time
    last_site = 0
    if rule[1].lower() == 'c':
        for i in cut_unmod_sites:
            cut_unmod_peptides.append(old_pep_seq[last_site:i+1])
            last_site = i + 1
        cut_unmod_peptides.append(old_pep_seq[last_site:])
    elif rule[1].lower() == 'n':
        for i in cut_unmod_sites:
            cut_unmod_peptides.append(old_pep_seq[last_site:i])
            last_site = i
        cut_unmod_peptides.append(old_pep_seq[last_site:])

    # Puts the break points for modified sequence in a list
    breaks = finditer(rule[0], pep_mod_seq)
    cut_sites = []
    cut_peptides = []
    for breakp in breaks:
        cut_sites.append(breakp.start())
    # Adds the cut peptides in one at a time
    last_site = 0
    if rule[1].lower() == 'c':
        for i in cut_sites:
            cut_peptides.append(pep_mod_seq[last_site:i+1])
            last_site = i + 1
        cut_peptides.append(pep_mod_seq[last_site:])
    elif rule[1].lower() == 'n':
        for i in cut_sites:
            cut_peptides.append(pep_mod_seq[last_site:i])
            last_site = i
        cut_peptides.append(pep_mod_seq[last_site:])
    #print(cut_peptides)

    # This gets the range of the first and last pep_slice with a PTM
    start = 0
    while start < len(cut_peptides):
        if any(user_PTM in cut_peptides[start] for user_PTM in user_PTMs):
            break
        else:
            start += 1
    
    end = len(cut_peptides) - 1
    while end > start:
        if any(user_PTM in cut_peptides[end] for user_PTM in user_PTMs):
            break
        else:
            end -= 1

    if start == end:
        new_pep_mod_seq = cut_peptides[start]
        new_pep_seq = cut_unmod_peptides[start]
    else:
        new_pep_mod_seq = ''.join(cut_peptides[start:end+1])
        new_pep_seq = ''.join(cut_unmod_peptides[start:end+1])
    return new_pep_mod_seq, new_pep_seq, len(''.join(cut_unmod_peptides[0:start]))  #Check this isn't modded instead

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
        quant_range = None
        var_mod_map = dict()
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
                    header = [i for i in header if i]
                    second_row = next(file_reader)  # need to normalize headers
                    #print(header)
                    #print(second_row)
                    quant_range = (len(header), len(second_row)-1)  #range where there will be quant; appears to skip some
                    header.extend(['SampleQuant' 
                                   for i in range(len(second_row) - len(header))])
                    out_writer.writerows([header, second_row])
                    for row in file_reader:
                        out_writer.writerow(row)
                    break  # Exit loop as we have all of the header information
        print("Finished Parsing")
        if quant_range is None:
            return -5, None, None, None
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
    #wd = os.getcwd()
    fasta_dir = os.path.dirname(os.path.abspath(proteome_fasta_filepath))
    pep_dict_file = f"{proteome_fasta_filepath.split('.')[0]}.pepdict"
    #print(os.listdir(fasta_dir))
    if not any([x == ntpath.basename(pep_dict_file) for x in os.listdir(fasta_dir)]):
        print("Pepdict not found in the FASTA directory. Trying to generate from the FASTA file")
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
    print("\033[95m {}\033[00m".format("Peptide Dictionary Loaded."))
    return pep_dict


def digest_if_needed(proteome_fasta_filepath: str, enzyme: str, mcleave: bool, minlen=5, maxlen=55, 
                     m_excise="Both"):
    """Runs digest and dumps pickle file (dict)"""
    seq_list = []  # (gene, organism, sequence, unpid)
    seq_list = fasta_producer(proteome_fasta_filepath, seq_list)
    rule = cleave_rule_determination(enzyme)
    pep_dict = defaultdict(set)  # NOTE: this won't work for large fasta/sequence databases as it's in ram; ~1GB for human
    for i in seq_list:
        gene, _organism, seq, unpid, protein_name = i
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
            digest_results = digest(seq, minlen, maxlen, mcleave, m_excise,
                                    li_swap=True, rule_to_use=rule,
                                    reverse=False)
        for result in digest_results:
            peptide, position = result
            peptide = peptide.replace("L", "I")
            pep_dict[peptide].add((gene, position, unpid, protein_name))  # careful with site specificity here
    with open(f"{proteome_fasta_filepath.split('.')[0]}.pepdict", "wb") as w1:
        pickle.dump(pep_dict, w1)
    del pep_dict  # Free memory for later
