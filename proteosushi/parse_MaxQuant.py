#!/usr/bin/env python

"""parse_MaxQuant.py: using helper functions, the script parses MaxQuant result files"""

from collections import defaultdict
import csv
import os.path
from re import match
from time import sleep

from .ps_utilities import clean_pep_seq
from .ps_utilities import load_pepdict, parse_maxquant_summary


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
        #modIndices = [i for i, s in enumerate(header) if "Probabilities" in s]
        modified_index = header_lower.index("modified sequence")
        mods_index = header_lower.index("modifications")
        #PTMs = [ptm.split(" Probabilities")[0] for ptm in header if "Probabilities" in ptm]
        start_index = header_lower.index([ptm for ptm in header_lower if "score diffs" in ptm][-1]) + 1
        end_index = header_lower.index("missed cleavages")
        PTMs = header[start_index:end_index]
        for row in tsv_reader:
            pep_mods = dict()
            #sequence = row[seq_index].replace('L', 'I')
            #if sequence == "MAPACQIIR":
            #    print("peptide found")
            for PTM in PTMs:
                modified_sequence = row[modified_index].strip('_')
                new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, modified_sequence, user_PTMs, row[seq_index])
                mods_B4_mod = 0
                # Grabs the first two letters of PTM that maxquant uses in mod_seq
                abr_ptm = PTM[:2].lower()  # NOTE: I could possibly imagine this not working in certain circumstances...
                # while this PTM is in the sequence we are looking at and listed among the PTMs for this peptide
                while PTM in row[mods_index] and abr_ptm in modified_sequence:  # This seems a little sketchy. What if there are two different PTMs with the same first 2 letters?
                    other_PTMs = PTMs.copy()
                    other_PTMs.remove(PTM)
                    # Gets the index of the PTM within the sequence (that is right here XX(|ptm)XX)
                    for other_ptm in other_PTMs:
                        mod_index = modified_sequence.index(abr_ptm)
                        abr_other_ptm = other_ptm[:2].lower()
                        while (abr_other_ptm in modified_sequence and 
                               modified_sequence.index(abr_other_ptm) < mod_index):
                            #This takes the ptms, one by one, until each before is removed
                            mods_B4_mod += 1
                            ptm_index = modified_sequence.index(abr_other_ptm)
                            modified_sequence = (modified_sequence[:ptm_index]
                                               + modified_sequence[ptm_index + len(abr_other_ptm):])  # NOTE: removes the mod, but not the parenthesis
                    mod_index = modified_sequence.index(abr_ptm)
                    pep_mods[modified_sequence.index(abr_ptm) - (mods_B4_mod * 2)] = PTM
                    modified_sequence = (modified_sequence[:mod_index]
                                       + modified_sequence[mod_index + len(abr_ptm):])
                    mods_B4_mod += 1
            for loc in sorted(pep_mods):
                try:
                    if not tuple((pep_mods[loc], loc - 2)) in mod_dict[row[modified_index].strip('_')]:
                        mod_dict[row[modified_index].strip('_')].append(tuple((pep_mods[loc],
                                                        loc - 2)))
                except KeyError:
                    mod_dict[row[modified_index].strip('_')] = [tuple((pep_mods[loc],
                                                loc - 2))]
    return mod_dict, PTMs
'''
def __prompt_dir() -> str:
    """Prompts the user for the directory with the MaxQuant files
    Returns:
        str -- the path to the directory
    """
    print("\033[96m {}\033[00m".format("Please provide the path to the MaxQuant files."))
    print("\033[96m {}\033[00m".format("For example: C:/experiment1/Combined/combined/txt/, or"))
    MQ_dir = input("\033[96m {}\033[00m".format("/home/[user]/Documents/experiment1/txt/\n"))
    if not os.path.exists(f"{MQ_dir}evidence.txt"):
        print("\033[91m {}\033[00m".format("Invalid path, try again!\n"))
        sleep(.5)
        return __prompt_dir()
    return MQ_dir

def __prompt_PTMs(PTMs: list) -> list:
    """prompts the user for which PTMs should be used
    Arguments:
        PTMs {list} -- a list of PTM names
    Returns:
        list -- a list of PTMs specified by the user
    """
    print("\033[96m {}\033[00m".format("Possible PTMs for analysis include:"))
    print("\033[93m {}\033[00m".format("\n ".join([f"[{i}] {ptm}" for i, ptm in enumerate(PTMs)])))
    print("\033[96m {}\033[00m".format("Enter the number for each of the PTMs of interest separated by a comma."))
    ptm_indices = input("\033[96m {}\033[00m".format("For example: 0,1,2\n"))
    if match(r"{^0-9|,}", ptm_indices):  # r"[0-9](?:,[0-9])*"
        print("\033[91m {}\033[00m".format("Invalid input, try again!\n"))
        sleep(.5)
        return __prompt_PTMs(PTMs)  # NOTE: the end of this recursion is getting it right
    try:
        mod_PTMs = [PTMs[int(x.strip())] for x in ptm_indices.strip().split(',')]
    except (ValueError, IndexError):
        print("\033[91m {}\033[00m".format("Invalid input, try again!\n"))
        sleep(.5)
        return __prompt_PTMs(PTMs)
    return mod_PTMs
'''
def compile_data_maxquant(search_engine_filepath: str, user_PTMs: list, cleave_rule: tuple) -> list:
    """Takes the lists and dictionaries needed to parse files
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        dict -- mod_dict
        int -- intensity_start (maybe just skyline)
        str -- evidence_filename
        list -- PTMs the user has chosen to be used
    """
    MQ_dir = search_engine_filepath #__prompt_dir()
    sum_file = os.path.join(MQ_dir, "summary.txt")
    enzyme, max_missed_cleaves = parse_maxquant_summary(sum_file)
    evidence_filename = os.path.join(MQ_dir, "evidence.txt")
    mod_dict, PTMs = parse_evidence(evidence_filename, user_PTMs, cleave_rule)  # This grabs the PTMs from the evidence file
    PTMs = user_PTMs#__prompt_PTMs(PTMs)
    evidence_file = open(evidence_filename, 'r')
    tsv_reader = csv.reader(evidence_file, delimiter='\t', quotechar='"')
    header = next(tsv_reader)
    header_lower = [s.lower() for s in header]
    sequence = header_lower.index("sequence")
    mod_seq = header_lower.index("modified sequence")
    reporter_intensities = [i for i, h in enumerate(header) if "intensit" in h.lower() and 
                                                               not "max intensity m/z" in h.lower()]
    #reporter_intensity = header.index("Intensity")
    evidence_file.close()
    return sequence, mod_seq, mod_dict, reporter_intensities, evidence_filename, PTMs
#EOF