"""ps_utilities.py: helper functions used by parsers in proteoSushi"""

__author__ = "Arshag Mooradian, Rob Seymour"
__email__ = "rseymour@wustl.edu"

import pickle
import csv
from collections import defaultdict
import ntpath
import os
from re import finditer
try:
    from .ps_digest import digest, fasta_producer, cleave_rule_determination
except ImportError:
    from ps_digest import digest, fasta_producer, cleave_rule_determination


def clean_localization_pep_seq(rule: tuple, pep_mod_seqs: list, user_PTMs: list, old_pep_seq: str) -> tuple:
    """Removes the missed cleavages without selected PTMs

    Arguments:
        rule {tuple} -- contains info about the protease
        pep_mod_seq {str} -- the peptide sequence with modifications
        user_PTMs {list} -- a list of PTMs that the user has selected for Analysis
    Returns:
        list -- the cleaned PTM localization peptide
        str -- the cleaned peptide
        int -- length of the front part of peptide that was removed
    """
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
    #if old_pep_seq == "KACLNPASPIVK":
    #    input("peptide split " + str(cut_unmod_peptides))

    PTM_cut_peptides = []
    for pep_mod_seq in pep_mod_seqs:
        if pep_mod_seq == '':
            PTM_cut_peptides.append('')
            continue
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
        PTM_cut_peptides.append(cut_peptides)
    #if old_pep_seq == "KACLNPASPIVK":
    #    input("peptide split " + str(PTM_cut_peptides))
    try:
        cut_peptides = [s for s in PTM_cut_peptides if s != ""][0]
    except IndexError:
        return None, None, None
    #if old_pep_seq == "KACLNPASPIVK":
    #    print(cut_peptides)
    first_start = len(cut_peptides)
    last_end = 0
    for cut_peptides in [s for s in PTM_cut_peptides if s != ""]:
    # This gets the range of the first and last pep_slice with a PTM
        start = 0
        while start < len(cut_peptides):
            if '(' in cut_peptides[start] or '[' in cut_peptides[start]:
                break
            else:
                start += 1
        if start < first_start:
            first_start = start

        end = len(cut_peptides) - 1
        while end > start:
            if '(' in cut_peptides[end] or '[' in cut_peptides[end]:
                break
            else:
                end -= 1
        if end > last_end:
            last_end = end
    #if old_pep_seq == "KACLNPASPIVK":
    #    input("start end " + str(first_start) + ' ' + str(last_end))
    if first_start == last_end:
        new_pep_mod_seqs = [peptide[first_start] for peptide in PTM_cut_peptides if peptide != '']
        new_pep_seq = cut_unmod_peptides[first_start]
    else:
        new_pep_mod_seqs = [''.join(peptide[first_start:last_end+1]) for peptide in PTM_cut_peptides if peptide != '']
        new_pep_seq = ''.join(cut_unmod_peptides[first_start:last_end+1])
    return new_pep_mod_seqs, new_pep_seq, len(''.join(cut_unmod_peptides[0:first_start]))

def clean_pep_seq(rule: tuple, pep_mod_seq: str, user_PTMs: list, old_pep_seq: str, is_MQ = False) -> str:
    """Removes the missed cleavages without selected PTMs

    Arguments:
        rule {tuple} -- contains info about the protease
        pep_mod_seq {str} -- the peptide sequence with modifications
        user_PTMs {list} -- a list of PTMs that the user has selected for Analysis
        old_pep_seq {str} -- the unmodified peptide sequence
        is_MQ {bool} -- whether the search engine is MQ
    Returns:
        str -- the cleaned peptide
    """
    # Puts break points for unmodified sequence in a list
    breaks = finditer(rule[0], old_pep_seq)
    cut_unmod_sites = []
    cut_unmod_peptides = []
    for breakp in breaks:
        cut_unmod_sites.append(breakp.start())
    # Adds the unmodified cut peptides in the list one at a time
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
    # Adds the cut peptides in the list one at a time
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

    # This gets the range of the first and last pep_slice with a PTM
    start = 0
    while start < len(cut_peptides):
        if is_MQ:
            if any(user_PTM[:2].lower() in cut_peptides[start] for user_PTM in user_PTMs):
                break
            else:
                start += 1
        else:
            if any(user_PTM in cut_peptides[start] for user_PTM in user_PTMs):
                break
            else:
                start += 1
    
    end = len(cut_peptides) - 1
    while end > start:
        if is_MQ:
            if any(user_PTM[:2].lower() in cut_peptides[end] for user_PTM in user_PTMs):
                break
            else:
                end -= 1
        else:
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
    new_start = start
    while len(new_pep_seq) < 6:
        if new_start <= 0:
            end += 1
            new_pep_mod_seq = ''.join(cut_peptides[new_start:end+1])
            new_pep_seq = ''.join(cut_unmod_peptides[new_start:end+1])
        elif end >= len(cut_unmod_peptides) - 1:
            new_start -= 1
            new_pep_mod_seq = ''.join(cut_peptides[new_start:end+1])
            new_pep_seq = ''.join(cut_unmod_peptides[new_start:end+1])
        elif len(cut_unmod_peptides[new_start-1]) >= len(cut_unmod_peptides[end+1]):
            new_start -= 1
            new_pep_mod_seq = ''.join(cut_peptides[new_start:end+1])
            new_pep_seq = ''.join(cut_unmod_peptides[new_start:end+1])
        else:
            end += 1
            new_pep_mod_seq = ''.join(cut_peptides[new_start:end+1])
            new_pep_seq = ''.join(cut_unmod_peptides[new_start:end+1])
    #input(f"old {pep_mod_seq}, new {new_pep_mod_seq}")
    return new_pep_mod_seq, new_pep_seq, len(''.join(cut_unmod_peptides[0:start]))

def parse_maxquant_summary(infile: str) -> list:
    """Pulls out information needed from the MaxQuant files
    Arguments:
        infile {str} -- the name of the MaxQuant summary file that will be used
    Returns:
        str -- the name of the protease used (usually trypsin)
        int -- the maximum allowed number of missed cleavages
    """
    mqfile = open(infile, 'r')
    tsvReader = csv.reader(mqfile, delimiter='\t', quotechar='"')
    protease = ""
    PTMs = []
    maxMissedCleaves = 0
    # From proteinGroups: Need Intensity[DM], PTM data[GS-GX], Potential Contaminant(?)[GK]
    # From summary: protease[C], Variable Modifications[H][I-J?], Max. missed cleavages[O]
    # From parameters: Modifications included in protein quantification,
    header = next(tsvReader)  # Skips the header line
    row = next(tsvReader)
    protease = row[header.index("Enzyme")]
    assert protease != ""
    maxMissedCleaves = int(row[header.index("Max. missed cleavages")])
    mqfile.close()
    return protease.lower(), maxMissedCleaves

def parse_mascot(in_file: str) -> list:
    """Pulls out protease, quantitative range, and variable modifications, along
    with the quant data from a Mascot file.
    Arguments:
        in_file {str} -- a string for the name of the mascot file
    Returns:
        list -- protease {str}, quant_range {tuple}, var_mod_map {dict},
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
        protease = ""
        missed_cleaves = -1
        quant_range = None
        var_mod_map = dict()
        for row in file_reader:
            if row:  # lots of blank rows
                # Pull data we need
                if row[0] == "Enzyme":
                    protease = row[1].lower()
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
        return protease, quant_range, var_mod_map, missed_cleaves


def load_pepdict(proteome_fasta_filepath: str, protease: str, missed_cleavages: int) -> dict:
    """Looks for pep_dict; if it doesn't find, it looks for a FASTA to make one.
    Fasta would be in current directory but this can be changed.

    Arguments:
        proteome_fasta_filepath {str} -- the filepath for the proteome fasta file
        protease {str} -- protease used in the sample digestion
        missed_cleavages {int} -- the max number of allowed missed cleavages
    Returns:
        dict -- a dictionary that connects peptide sequence to gene, position, protein name and unpid
    """
    pep_dict = None
    #wd = os.getcwd()
    fasta_dir = os.path.dirname(os.path.abspath(proteome_fasta_filepath))
    pep_dict_file = f"{proteome_fasta_filepath.split('.')[0]}_{protease.replace('/','^')}.pepdict"
    #print(os.listdir(fasta_dir))
    if not any([x == ntpath.basename(pep_dict_file) for x in os.listdir(fasta_dir)]):
        print("Pepdict not found in the FASTA directory. Trying to generate from the FASTA file")
        digest_if_needed(proteome_fasta_filepath, protease, missed_cleavages)
    
    print(f"Loading {pep_dict_file} into memory")
    with open(pep_dict_file, 'rb') as f:
        pep_dict = pickle.load(f)
    print("\033[95m {}\033[00m".format("Peptide Dictionary Loaded."))
    return pep_dict


def digest_if_needed(proteome_fasta_filepath: str, protease: str, missed_cleavages: int, minlen=6, maxlen=55, 
                     m_excise="Both"):
    """Runs digest and dumps pickle file (dict)"""
    seq_list = []  # (gene, organism, sequence, unpid)
    seq_list = fasta_producer(proteome_fasta_filepath, seq_list)
    rule = cleave_rule_determination(protease)
    pep_dict = defaultdict(set)  # NOTE: this won't work for large fasta/sequence databases as it's in ram; ~1GB for human
    for i in seq_list:
        gene, _organism, seq, unpid, protein_name = i
        if m_excise == "Both":
            # Run digest results twice and combine; inefficient but gets job done
            digest_results_true = digest(seq, minlen, maxlen, missed_cleavages, True,
                                         li_swap=True, rule_to_use=rule,
                                         reverse=False)
            digest_results_false = digest(seq, minlen, maxlen, missed_cleavages, False,
                                          li_swap=True, rule_to_use=rule,
                                          reverse=False)
            digest_results = set().union(digest_results_true,
                                         digest_results_false)
        else:
            digest_results = digest(seq, minlen, maxlen, missed_cleavages, m_excise,
                                    li_swap=True, rule_to_use=rule,
                                    reverse=False)
        for result in digest_results:
            peptide, position = result
            peptide = peptide.replace("L", "I")
            pep_dict[peptide].add((gene, position, unpid, protein_name))  # careful with site specificity here
    with open(f"{proteome_fasta_filepath.split('.')[0]}_{protease.replace('/','^')}.pepdict", "wb") as w1:
        pickle.dump(pep_dict, w1)
    del pep_dict  # Free memory for later
