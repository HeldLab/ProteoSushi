#!/usr/bin/env python

"""parse_MaxQuant.py: using helper functions, the script parses MaxQuant result files"""

from collections import defaultdict
import csv
import os.path
from re import findall, finditer, match
from time import sleep
try:
    from .ps_utilities import clean_localization_pep_seq, clean_pep_seq
    from .ps_utilities import load_pepdict, parse_maxquant_summary
except ImportError:
    from ps_utilities import clean_localization_pep_seq, clean_pep_seq
    from ps_utilities import load_pepdict, parse_maxquant_summary


def get_MQ_PTMs(filename: str) -> list:
    """Parses evidence file for just the PTMs available for selection
    Arguments:
        filename {str} -- name of the evidence file with the location
    Returns:
        list -- a list of the PTMs we can use for the analysis
    """
    PTMs = list()
    with open(filename, 'r') as evidence_file:
        tsv_reader = csv.reader(evidence_file, delimiter='\t', quotechar='"')
        header = next(tsv_reader)
        header_lower = [s.lower() for s in header]
        start_index = header_lower.index([ptm for ptm in header_lower if "score diffs" in ptm][-1]) + 1
        end_index = header_lower.index("missed cleavages")
        PTMs = header[start_index:end_index]
    return PTMs


def parse_evidence_localization(filename: str, user_PTMs: list, cleave_rule: tuple, localization_threshold: float) -> dict:
    """parses the evidence file and builds a dictionary for PTMs
    Arguments:
        filename {str} -- name of the evidence file with the location
        user_PTMs {list} -- list of user-chosen PTMs for analysis
        cleave_rule {tuple} -- information on the protease for cleaving in silico
        localization_threshold {float} -- the user-chosen threshold for localization probability
    Returns:
        dict -- a dictionary of the PTMs with sequence as key and a list of tuples with 
        PTM/location as value
        list -- a list of the PTMs we can use for the analysis
    """
    mod_dict = dict()
    PTMs = list()
    with open(filename, 'r') as evidence_file:
        tsv_reader = csv.reader(evidence_file, delimiter='\t', quotechar='"')
        header = next(tsv_reader)
        header_lower = [s.lower() for s in header]
        seq_index = header_lower.index("sequence")
        mod_indices = [i for i, s in enumerate(header) if "Probabilities" in s]
        modified_index = header_lower.index("modified sequence")
        PTMs = [ptm.split(" Probabilities")[0] for ptm in header if "Probabilities" in ptm]
        # Go through each line of the evidence file
        for row in tsv_reader:
            mod_seq = row[modified_index].strip('_')
            if len(row[seq_index]) < 6:
                continue
            loc_seqs = [row[i] if header[i].split(" Probabilities")[0] in user_PTMs else "" for i in mod_indices]
            new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, 
                                                                        mod_seq, 
                                                                        user_PTMs, 
                                                                        row[seq_index],
                                                                        True)
            new_loc_seqs, new_pep_seq, missed_cleave_fix = clean_localization_pep_seq(cleave_rule, 
                                                                                      loc_seqs, 
                                                                                      user_PTMs, 
                                                                                      row[seq_index])
            
            if new_loc_seqs is None:
                continue
            # Go through each of the PTM localization sequences
            index = 0
            while index < len(new_loc_seqs):
                if new_loc_seqs[index] == "":
                    index += 1
                    continue
                loc_iter = finditer(r"(\[.+?\])|(\(.+?\))", new_loc_seqs[index])
                loc_indices = []
                for loc_i in loc_iter:
                    loc_indices.append(loc_i.start())
                
                try:
                    loc_probs = findall(r"(\[.+?\])|(\(.+?\))", new_loc_seqs[index]).remove('')
                except ValueError:
                    loc_probs = findall(r"(\[.+?\])|(\(.+?\))", new_loc_seqs[index])
                
                index_correction = 0
                j = 0
                # Go through each of the PTM sites
                while j < len(loc_probs):
                    # Skips sites below the threshold
                    if float(loc_probs[j][1].replace('(','').replace(')','').replace('[','').replace(']','')) < localization_threshold:
                        index_correction += len(loc_probs[j][1])
                        j += 1
                        continue
                    try:  # The - 1 is a hotfix as I'm not sure why it is just 1 off.
                        mod_dict[new_mod_seq].append(tuple((PTMs[index], loc_indices[j] - index_correction - 1)))
                        index_correction += len(loc_probs[j][1])
                        j += 1
                    except KeyError:
                        mod_dict[new_mod_seq] = [tuple((PTMs[index], loc_indices[j] - index_correction - 1))]
                        index_correction += len(loc_probs[j][1])
                        j += 1
                index += 1
    return mod_dict, PTMs

def parse_evidence(filename: str, user_PTMs: list, cleave_rule: tuple) -> dict:
    """parses the evidence file and builds a dictionary for PTMs
    Arguments:
        filename {str} -- name of the evidence file with the location
    Returns:
        dict -- a dictionary of the PTMs with sequence as key and a list of tuples with 
        PTM/location as value
        list -- a list of the PTMs we can use for the analysis
    """
    mod_dict = dict()
    PTMs = list()
    with open(filename, 'r') as evidence_file:
        tsv_reader = csv.reader(evidence_file, delimiter='\t', quotechar='"')
        header = next(tsv_reader)
        header_lower = [s.lower() for s in header]
        seq_index = header_lower.index("sequence")
        modified_index = header_lower.index("modified sequence")
        mods_index = header_lower.index("modifications")
        start_index = header_lower.index([ptm for ptm in header_lower if "score diffs" in ptm][-1]) + 1
        end_index = header_lower.index("missed cleavages")
        PTMs = header[start_index:end_index]
        for row in tsv_reader:
            pep_mods = dict()
            modified_sequence = row[modified_index].strip('_')
            new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, modified_sequence, user_PTMs, row[seq_index])
            abr_ptms = [PTM[:2].lower() for PTM in user_PTMs]

            mods = findall(r"(\[.+?\])|(\(.+?\))", new_mod_seq)
            breaks = finditer(r"(\[.+?\])|(\(.+?\))", new_mod_seq)
            cut_sites = []
            for breakp in breaks:
                cut_sites.append(breakp.start())
            # This fixes the indices of each mod to connect with the unmodified sequence
            correction = 0
            fixed_indices = []
            num_of_mods = 0
            for site in cut_sites:
                fixed_indices.append(site - correction - 1)  # TODO: Check this!
                correction += len(mods[num_of_mods][1])
                num_of_mods += 1

            #if row[seq_index] == "AASCVLLHTGQK":
            #    input(str(mods))
            i = 0
            while i < len(cut_sites):
                if not mods[i][1][1:-1] in abr_ptms:
                    #print(f"{mods[i][1][1:-1]} not in {abr_ptms}")
                    i += 1
                    continue
                try:
                    mod_dict[new_mod_seq].append(tuple(([ptm for ptm in user_PTMs if mods[i][1][1:-1] == ptm[:2].lower()][0], 
                                                        fixed_indices[i])))
                except KeyError:
                    mod_dict[new_mod_seq] = [tuple(([ptm for ptm in user_PTMs if mods[i][1][1:-1] == ptm[:2].lower()][0], 
                                                    fixed_indices[i]))]
                i += 1

            '''
            for PTM in PTMs:
                modified_sequence = row[modified_index].strip('_')
                new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, modified_sequence, user_PTMs, row[seq_index])
                mods_B4_mod = 0
                # Grabs the first two letters of PTM that maxquant uses in mod_seq
                abr_ptm = PTM[:2].lower()  # NOTE: I could possibly imagine this not working in certain circumstances...
                # while this PTM is in the sequence we are looking at and listed among the PTMs for this peptide
                while PTM in row[mods_index] and abr_ptm in new_mod_seq:  # This seems a little sketchy. What if there are two different PTMs with the same first 2 letters?
                    other_PTMs = PTMs.copy()
                    other_PTMs.remove(PTM)
                    # Gets the index of the PTM within the sequence (that is right here XX(|ptm)XX)
                    for other_ptm in other_PTMs:
                        mod_index = new_mod_seq.index(abr_ptm)
                        abr_other_ptm = other_ptm[:2].lower()
                        while (abr_other_ptm in new_mod_seq and 
                               new_mod_seq.index(abr_other_ptm) < mod_index):
                            #This takes the ptms, one by one, until each before is removed
                            mods_B4_mod += 1
                            ptm_index = new_mod_seq.index(abr_other_ptm)
                            new_mod_seq = (new_mod_seq[:ptm_index]
                                               + new_mod_seq[ptm_index + len(abr_other_ptm):])  # NOTE: removes the mod, but not the parenthesis
                    mod_index = new_mod_seq.index(abr_ptm)
                    pep_mods[new_mod_seq.index(abr_ptm) - (mods_B4_mod * 2)] = PTM
                    new_mod_seq = (new_mod_seq[:mod_index]
                                       + new_mod_seq[mod_index + len(abr_ptm):])
                    mods_B4_mod += 1
            
            for loc in sorted(pep_mods):
                try:
                    if not tuple((pep_mods[loc], loc - 2)) in mod_dict[new_mod_seq]:
                        mod_dict[new_mod_seq].append(tuple((pep_mods[loc],
                                                     loc - 2 - missed_cleave_fix)))
                except KeyError:
                    mod_dict[new_mod_seq] = [tuple((pep_mods[loc],
                                             loc - 2 - missed_cleave_fix))]
            '''
    return mod_dict, PTMs

def compile_localization_data_maxquant(search_engine_filepath: str, user_PTMs: list, 
                                       cleave_rule: tuple, localization_threshold: float) -> list:
    """Takes the lists and dictionaries needed to parse files
    Arguments:
        search_engine_filepath {str} -- filepath to the SE directory
        user_PTMs {list} -- a list of user-selected PTMs for analysis
        cleave_rule {tuple} -- the protease cleave rule
        localization_threshold {float} -- user-chosen threshold for localization probability
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        list -- list of indices of modified sequence w/ localization prob
        dict -- mod_dict
        int -- intensity_start (maybe just skyline)
        str -- evidence_filename
    """
    MQ_dir = search_engine_filepath
    sum_file = os.path.join(MQ_dir, "summary.txt")
    protease, max_missed_cleaves = parse_maxquant_summary(sum_file)
    evidence_filename = os.path.join(MQ_dir, "evidence.txt")
    mod_dict, PTMs = parse_evidence_localization(evidence_filename, user_PTMs, cleave_rule, localization_threshold)  # This grabs the PTMs from the evidence file
    PTMs = user_PTMs
    evidence_file = open(evidence_filename, 'r')
    tsv_reader = csv.reader(evidence_file, delimiter='\t', quotechar='"')
    header = next(tsv_reader)
    header_lower = [s.lower() for s in header]
    sequence = header_lower.index("sequence")
    mod_seq = header_lower.index("modified sequence")
    loc_seqs = [i for i, h in enumerate(header) if " probabilities" in h.lower()]
    reporter_intensities = [i for i, h in enumerate(header) if "intensit" in h.lower() and 
                                                               not "max intensity m/z" in h.lower()]
    evidence_file.close()
    return sequence, mod_seq, loc_seqs, mod_dict, reporter_intensities, evidence_filename

def compile_data_maxquant(search_engine_filepath: str, user_PTMs: list, cleave_rule: tuple) -> list:
    """Takes the lists and dictionaries needed to parse files
    Arguments:
        search_engine_filepath {str} -- filepath to the SE directory
        user_PTMs {list} -- a list of user-selected PTMs for analysis
        cleave_rule {tuple} -- the protease cleave rule
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        dict -- mod_dict
        int -- intensity_start (maybe just skyline)
        str -- evidence_filename
        list -- PTMs the user has chosen to be used
    """
    MQ_dir = search_engine_filepath
    sum_file = os.path.join(MQ_dir, "summary.txt")
    protease, max_missed_cleaves = parse_maxquant_summary(sum_file)
    evidence_filename = os.path.join(MQ_dir, "evidence.txt")
    mod_dict, PTMs = parse_evidence(evidence_filename, user_PTMs, cleave_rule)  # This grabs the PTMs from the evidence file
    PTMs = user_PTMs
    evidence_file = open(evidence_filename, 'r')
    tsv_reader = csv.reader(evidence_file, delimiter='\t', quotechar='"')
    header = next(tsv_reader)
    header_lower = [s.lower() for s in header]
    sequence = header_lower.index("sequence")
    mod_seq = header_lower.index("modified sequence")
    reporter_intensities = [i for i, h in enumerate(header) if "intensit" in h.lower() and 
                                                               not "max intensity m/z" in h.lower()]
    evidence_file.close()
    return sequence, mod_seq, mod_dict, reporter_intensities, evidence_filename, PTMs
#EOF