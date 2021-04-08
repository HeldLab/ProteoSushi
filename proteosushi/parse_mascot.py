"""parse_mascot.py: functions used to parse the mascot file and read it into memory"""

import csv
from collections import defaultdict
import os
from re import match
from time import sleep
try:
    from .ps_utilities import clean_pep_seq, parse_mascot, load_pepdict #most of the parsing actually happens here 
except ImportError:
    from ps_utilities import clean_pep_seq, parse_mascot, load_pepdict


def __create_mod_dict(filename: str, mod_ids: list, var_mod_map: dict, cleave_rule: tuple, PTMs: list) -> dict:
    """Makes the PTM dictionary that connects a peptide sequence with its modifications
    Arguments:
        filename {str} -- the name of the output file
        mod_ids {list} -- a list of mod ids that mascot uses and the user chose
        var_mod_map {dict} -- connects the PTM to the Mascot index

    Returns:
        dict -- a dictionary of sequences with the type and location of modifications
    """
    modDict = {}
    inv_mod_map = {v:k for k, v in var_mod_map.items()}
    with open(filename, 'r') as mascot_output:
        file_reader = csv.reader(mascot_output, quotechar='"')
        header = next(file_reader)
        sequence = header.index("pep_seq")
        var_mods = header.index("pep_var_mod_pos")

        for row in file_reader:
            #Mascot has duplicate rows; check for these
            pep_seq = row[sequence].replace("L","I")
            if len(pep_seq) < 6:  # Skip any peptides smaller than 6
                continue
            #print(row[sequence])
            mods = row[var_mods]
            if not mods:
                continue
            pep_mod_seq = mods.split('.')[1]
            i = len(pep_mod_seq) - 1
            inv_mod_dict = {v:k for k, v in var_mod_map.items()}
            new_seq = row[sequence]
            while i >= 0:
                if pep_mod_seq[i] != "0":  # If there is a mod
                    new_seq = row[sequence][:i+1] + '(' + inv_mod_dict[pep_mod_seq[i]] + ')' + new_seq[i+1:]
                i -= 1
            pep_mod_seq = new_seq
            new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, pep_mod_seq, PTMs, row[sequence])
            if mods:
                n, center, c = [x for x in mods.split(".")]
                for i, aa in enumerate(list(center)): #split mods; note this excludes termini
                    if aa in mod_ids:
                        try:
                            if not tuple((str(aa), i)) in modDict[pep_mod_seq]:
                                modDict[new_mod_seq].append(tuple((inv_mod_map[aa], i-missed_cleave_fix)))
                        except KeyError:
                            modDict[new_mod_seq] = [tuple((inv_mod_map[aa], i-missed_cleave_fix))]
            #if row[sequence] == "AAFTECCQAADK":
            #    input(f"old {row[sequence]}, new {new_mod_seq}")
    return modDict


def compile_data_mascot(search_engine_filepath: str, PTMs: list, cleave_rule: tuple) -> list:
    """Takes the lists and dictionaries needed to parse files

    Arguments:
        search_engine_filepath {str} -- the filepath to the mascot output file
        PTMs {list} -- the list of PTMs that the user has chosen to use
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        dict -- mod_dict
        int -- intensity_start 
        str -- data_filename
        dict -- dictionary of the mods to the index
    """
    in_file = search_engine_filepath
    protease, quant_range, var_mod_map, missed_cleaves = parse_mascot(in_file)

    mods_for_quant = PTMs


    input_filename = "quant_parsed.csv"
    input_file =  open(input_filename, 'r')
    file_reader = csv.reader(input_file)
    header = next(file_reader)

    sequence = header.index('pep_seq') #assert this exists later on
    var_mods = header.index('pep_var_mod_pos')

    mod_dict = __create_mod_dict(input_filename, [var_mod_map[mod] for mod in mods_for_quant], var_mod_map, cleave_rule, PTMs)
    #print("3")
    if quant_range[1] <= quant_range[0]:
        return sequence, var_mods, mod_dict, None, input_filename, var_mod_map

    return sequence, var_mods, mod_dict, quant_range[0], input_filename, var_mod_map
#EOF