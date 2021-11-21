#!/usr/bin/env python

"""parse_generic.py: parses and rolls up the results from any search engine"""

from collections import defaultdict
import csv
import logging
import os.path
from re import findall, finditer, match
import sys
from time import sleep

try:
    from .ps_utilities import clean_pep_seq
    from .proteoSushi_constants import cleave_rules
    from .ps_utilities import load_pepdict
except ImportError:
    from ps_utilities import clean_pep_seq, clean_localization_pep_seq
    from proteoSushi_constants import cleave_rules
    from ps_utilities import load_pepdict


def get_PTMs(sky_filename: str) -> list:
    """Gets the list of PTMs available through parsing the generic output

    Arguments:
        sky_filename {str} -- the name of the generic output file
    Returns:
        list -- a list of PTMs in the output file
    """
    with open(sky_filename, 'r') as sky_output:
        tsv_reader = csv.reader(sky_output, quotechar='"')
        header = next(tsv_reader)
        header_lower = [s.lower() for s in header]
        try:
            seq_index = header_lower.index("peptide sequence")
        except ValueError:
            try:
                seq_index = header_lower.index("peptide")
            except ValueError:
                return -4
        mod_index = -1
        try:
            mod_index = header_lower.index("peptide modified sequence") #NOTE: apparently, this changes, so check here.
        except ValueError:
            try:
                mod_index = header_lower.index("modified sequence")
            except ValueError:
                try:
                    mod_index = header_lower.index("modified peptide")
                except ValueError:
                    print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the generic output file."))
                    print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                    return -4
        
        mod_list = []
        for row in tsv_reader:
            # This grabs the PTMs in each sequence and builds a list of all of them
            sequence = row[seq_index].replace("L","I")
            mod_seq = row[mod_index]
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", mod_seq)

            #If there are no PTMs in the sequence
            if len(mods) < 1:
                return -3  # This will be there error listed below, but handled in the GUI
                raise ValueError("A sequence in the Peptide Modified Sequence column is missing PTMs")

            for mod in mods[0]:
                if mod == '':
                    continue
                stripped_mod = mod
                if not stripped_mod in mod_list:
                    mod_list.append(stripped_mod)
    return mod_list

def create_loc_mod_dict(sky_filename: str, user_PTMs: list, cleave_rule: tuple, 
                        localization_threshold: float) -> dict:
    """Parses the generic output to retrieve the modification info

    Arguments:
        sky_filename {str} -- the name of the generic output file
        user_PTMs {list} -- a list of PTMs that the user wants to analyze
        cleave_rule {tuple} -- information needed to cleave peptides
        localization_threshold {float} -- the threshold for the lowest allowed localization
        probability; all localization prob is above this number
    Returns:
        dict -- a dictionary of sequences with the type and location of modifications
    """
    mod_dict = {}
    with open(sky_filename, 'r') as sky_output:
        tsv_reader = csv.reader(sky_output, quotechar='"')
        header = next(tsv_reader)
        header_lower = [s.lower() for s in header]
        try:
            seq_index = header_lower.index("peptide sequence")
        except ValueError:
            try:
                seq_index = header_lower.index("peptide")
            except ValueError:
                return -4
        mod_indices = [i for i, s in enumerate(header_lower) if "probabilities" in s]
        PTMs = [ptm.split(" Probabilities")[0] for ptm in header if "Probabilities" in ptm]
        mod_index = -1
        try:
            mod_index = header_lower.index("peptide modified sequence") #NOTE: apparently, this changes, so check here.
        except ValueError:
            try:
                mod_index = header_lower.index("modified sequence")
            except ValueError:
                try:
                    mod_index = header_lower.index("modified peptide")
                except ValueError:
                    logging.error("no generic peptide modified sequence")
                    print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the generic output file."))
                    print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                    return -4
        
        for row in tsv_reader:
            #print(row)
            # This grabs the PTMs in each sequence and builds a list of all of them
            #sequence = row[seq_index].replace("L","I")
            #if sequence == "FACAVVCIQK":
            #    print("here")
            mod_seq = row[mod_index]
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
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", mod_seq)
            if len(mods) < 1:
                raise ValueError("A sequence in the Peptide Modified Sequence column is missing PTMs")
            
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
                    #if row[seq_index] == "KACLNPASPIVK":
                    #    input("numbers: "+  str(loc_indices[j])+" - "+str(index_correction)+" - 1 - "+str(missed_cleave_fix))
                    # Skips sites below the threshold
                    if float(loc_probs[j][1].replace('(','').replace(')','').replace('[','').replace(']','')) < localization_threshold:
                        index_correction += len(loc_probs[j][1])
                        j += 1
                        continue
                    try:  # The - 1 is a hotfix as I'm not sure why it is just 1 off.
                        mod_dict[new_mod_seq].append(tuple((PTMs[index], loc_indices[j] - index_correction - 1)))# - missed_cleave_fix)))
                        index_correction += len(loc_probs[j][1])
                        j += 1
                    except KeyError:
                        mod_dict[new_mod_seq] = [tuple((PTMs[index], loc_indices[j] - index_correction - 1))]# - missed_cleave_fix))]
                        index_correction += len(loc_probs[j][1])
                        j += 1
                index += 1
    return mod_dict

def create_mod_dict(sky_filename: str, user_PTMs: list, cleave_rule: tuple) -> dict:
    """Parses the generic output to retrieve the modification info

    Arguments:
        sky_filename {str} -- the name of the generic output file
        user_PTMs {list} -- a list of PTMs that the user wants to analyze
        cleave_rule {tuple} -- information needed to cleave peptides
    Returns:
        dict -- a dictionary of sequences with the type and location of modifications
    """
    mod_dict = {}
    with open(sky_filename, 'r') as sky_output:
        logging.debug("Grabbing header info from search engine file")
        tsv_reader = csv.reader(sky_output, quotechar='"')
        header = next(tsv_reader)
        header_lower = [s.lower() for s in header]
        try:
            seq_index = header_lower.index("peptide sequence")
        except ValueError:
            try:
                seq_index = header_lower.index("peptide")
            except ValueError:
                return -4
        mod_index = -1
        try:
            mod_index = header_lower.index("peptide modified sequence") #NOTE: apparently, this changes, so check here.
        except ValueError:
            try:
                mod_index = header_lower.index("modified sequence")
            except ValueError:
                try:
                    mod_index = header_lower.index("modified peptide")
                except ValueError:
                    logging.error("generic file had no peptide modified sequence column")
                    print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the generic output file."))
                    print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                    return -4
        
        logging.debug("building a list with all PTMs in each sequence")
        for row in tsv_reader:
            # This grabs the PTMs in each sequence and builds a list of all of them
            mod_seq = row[mod_index]
            if len(row[seq_index]) < 6:
                continue
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", mod_seq)
            if len(mods) < 1:
                raise ValueError("A sequence in the Peptide Modified Sequence column is missing PTMs")
            
            new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, mod_seq, user_PTMs, row[seq_index])
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", new_mod_seq)  # NOTE that these grab the AA as well
            breaks = finditer(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", new_mod_seq)
            logging.debug(f"Mods are: {mods};\nbreaks are: {breaks}")
            cut_sites = []
            for breakp in breaks:
                cut_sites.append(breakp.start())
            # This fixes the indices of each mod to connect with the unmodified sequence
            logging.debug(f"fixing the indices of {mod_seq}")
            correction = 0
            fixed_indices = []
            num_of_mods = 0
            for site in cut_sites:
                fixed_indices.append(site - correction)  # TODO: Check this!
                correction += len(mods[num_of_mods][0]) - 1  # The -1 is because of the AA connected
                num_of_mods += 1
            logging.debug("Building the mod dict")
            i = 0
            while i < len(cut_sites):
                if not mods[i][0] in user_PTMs:
                    i += 1
                    continue
                try:
                    mod_dict[new_mod_seq].append(tuple((mods[i][0], fixed_indices[i])))  # TODO: Check This!
                except KeyError:
                    mod_dict[new_mod_seq] = [tuple((mods[i][0], fixed_indices[i]))]
                i += 1
    logging.debug("Finished building the mod_dict")
    return mod_dict



def compile_localization_data_generic(search_engine_filepath: str, user_PTMs: list, 
                                      cleave_rule: tuple, localization_threshold: float) -> list:
    """Takes the lists and dictionaries needed to parse files

    Arguments:
        search_engine_filepath {str} -- the filepath for the generic search engine output
        user_PTMs {list} -- the list of PTMs that the user will use for the analysis
        cleave_rule {tuple} -- information for how to cleave peptides
        localization_threshold {float} -- the threshold for the lowest allowed localization
        probability; all localization prob is above this number
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        dict -- pep_dict
        dict -- mod_dict
        int -- intensity_values
        file -- mito_genes (file)
        str -- sky_filename
        dict -- mascot var mod dict (blank)
    """
    sky_filename = search_engine_filepath
    logging.debug("Reading in the mod dict")
    mod_dict = create_loc_mod_dict(sky_filename, user_PTMs, cleave_rule, localization_threshold)
    logging.debug("Opening the search engine file")
    skyFile = open(sky_filename, 'r')
    tsv_reader = csv.reader(skyFile, quotechar='"')
    header = next(tsv_reader)
    header_lower = [s.lower() for s in header]
    try:
        sequence = header_lower.index("peptide sequence")
    except ValueError:
        try:
            sequence = header_lower.index("peptide")
        except ValueError:
            return -4
    pms = -1
    try:
        pms = header_lower.index("peptide modified sequence") #NOTE: apparently, this changes, so check here.
    except ValueError:
        try:
            pms = header_lower.index("modified sequence")
        except ValueError:
            try:
                pms = header_lower.index("modified peptide")
            except ValueError:
                logging.error("generic file had no peptide modified sequence")
                print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the generic output file."))
                print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                return -4
    logging.debug("Reading in localization indices and intensity values")
    localization_indices = [i for i, h in enumerate(header) if " probabilities" in h.lower()]
    intensity_values = [i for i, h in enumerate(header) if "intensit" in h.lower()]  # NOTE: This makes the assumption that the intensity columns are directly after the charge column
    if len(intensity_values) == 0:
        intensity_values = None
    return sequence, pms, localization_indices, mod_dict, intensity_values, sky_filename, None

def compile_data_generic(search_engine_filepath: str, user_PTMs: list, cleave_rule: tuple) -> list:
    """Takes the lists and dictionaries needed to parse files

    Arguments:
        search_engine_filepath {str} -- the filepath for the generic search engine output
        user_PTMs {list} -- the list of PTMs that the user will use for the analysis
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        dict -- pep_dict
        dict -- mod_dict
        int -- intensity_values
        file -- mito_genes (file)
        str -- sky_filename
        dict -- mascot var mod dict (blank)
    """
    sky_filename = search_engine_filepath
    logging.debug("Creating mod dict")
    mod_dict = create_mod_dict(sky_filename, user_PTMs, cleave_rule)
    logging.debug("opening search engine file")
    skyFile = open(sky_filename, 'r')
    tsv_reader = csv.reader(skyFile, quotechar='"')
    header = next(tsv_reader)
    header_lower = [s.lower() for s in header]
    try:
        sequence = header_lower.index("peptide sequence")
    except ValueError:
        try:
            sequence = header_lower.index("peptide")
        except ValueError:
            return -4
    pms = -1
    try:
        pms = header_lower.index("peptide modified sequence") #NOTE: apparently, this changes, so check here.
    except ValueError:
        try:
            pms = header_lower.index("modified sequence")
        except ValueError:
            try:
                pms = header_lower.index("modified peptide")
            except ValueError:
                logging.error("generic file had no peptide modified sequence")
                print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the generic output file."))
                print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                return -4
    #UNPID = -1
    #try:
    #    UNPID = header.index("Protein Accession")
    #except ValueError:
    #    print("\033[91m {}\033[00m".format("No uniprot ID column detected in the Skyline output file."))
    #    print("Proceeding regardless.\n")
    logging.debug("Reading in intensity values")
    intensity_values = [i for i, h in enumerate(header) if "intensit" in h.lower()]  # NOTE: This makes the assumption that the intensity columns are directly after the charge column
    if len(intensity_values) == 0:
        intensity_values = None
    #psm_contributions = defaultdict(int)
    #unmatchedPeps = 0
    #missingPTM = 0
    #totalSeqs = 0
    #unmatchedSequences = []
    #intensity_dict = dict()
    logging.debug("Finished compiling analysis data")
    return sequence, pms, mod_dict, intensity_values, sky_filename, None
#EOF