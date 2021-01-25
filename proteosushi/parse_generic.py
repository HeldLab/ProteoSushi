#!/usr/bin/env python

"""parse_generic.py: parses and rolls up the results from any search engine"""

from collections import defaultdict
import csv
import os.path
from re import findall, finditer, match
import sys
from time import sleep

from .ps_utilities import clean_pep_seq
from .proteoSushi_constants import cleave_rules
from .ps_utilities import load_pepdict


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

def create_mod_dict(sky_filename: str, user_PTMs: list, cleave_rule: tuple) -> dict:
    """Parses the generic output to retrieve the modification info

    Arguments:
        sky_filename {str} -- the name of the generic output file
        user_PTMs {list} -- a list of PTMs that the user wants to analyze
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
        
        for row in tsv_reader:
            #print(row)
            # This grabs the PTMs in each sequence and builds a list of all of them
            #sequence = row[seq_index].replace("L","I")
            #if sequence == "FACAVVCIQK":
            #    print("here")
            mod_seq = row[mod_index]
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", mod_seq)
            if len(mods) < 1:
                raise ValueError("A sequence in the Peptide Modified Sequence column is missing PTMs")
            
            new_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rule, mod_seq, user_PTMs, row[seq_index])
            #print(new_mod_seq)
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", new_mod_seq)
            breaks = finditer(r"(\w?\[.+?\])|(\w?\(.+?\(?.\)?\))", new_mod_seq)
            cut_sites = []
            for breakp in breaks:
                cut_sites.append(breakp.start())
            #print(cut_sites)
            # This fixes the indices of each mod to connect with the unmodified sequence
            correction = 0
            fixed_indices = []
            num_of_mods = 0
            for site in cut_sites:
                fixed_indices.append(site - correction)  # TODO: Check this!
                correction += len(mods[num_of_mods][0]) - 1
                num_of_mods += 1
            
            i = 0
            while i < len(cut_sites):
                #print(new_mod_seq)
                if not mods[i][0] in user_PTMs:
                    i += 1
                    continue
                #print(new_mod_seq)
                try:
                    mod_dict[new_mod_seq].append(tuple((mods[i][0], fixed_indices[i])))  # TODO: Check This!
                except KeyError:
                    mod_dict[new_mod_seq] = [tuple((mods[i][0], fixed_indices[i]))]
                i += 1
    #print(mod_dict)
    return mod_dict

def __chooseHit(genes_positions: list, mito_genes: list, annotDict: dict) -> list:
    """chooses which of the matched sequences to use. If there is one mito gene, it will be that one.
    If there are more than one non-mito, annotation score decides. If there are more than one mito, annotation score decides.
    Arguments:
        genes_positions {list} -- a list of tuples with gene and position info
        mito_genes {list} -- a list of strings of mito gene names
        annotDict {dict} -- a dictionary connecting genes to annotation scores
    Returns:
        bool -- whether the match(es) is(are) target genes
        *and*
        str -- the gene name of the chosen match
        int -- the start position of the chosen match
        *or*
        list -- a list of tuples with the gene name and start position of each match
    """
    mitoTups = list()
    nonMitoTups = list()
    for tup in genes_positions:
        if tup[0].upper() + '\n' in mito_genes:
            #print("Mito gene prioritized!")
            mitoTups.append(tup)
        else:
            nonMitoTups.append(tup)
    if len(mitoTups) == 1:
        return True, mitoTups[0]
    elif len(mitoTups) > 1:
        return True, __chooseTup(mitoTups, annotDict)
    elif len(nonMitoTups) == 1:
        return False, nonMitoTups[0]
    elif len(nonMitoTups) > 1:  # Only non-mito proteins (and more than 1)
        return False, __chooseTup(nonMitoTups, annotDict)
    assert False, "ERROR: chooseHit had 0 tuples sent in!"
    #return False, None, None

def __chooseTup(tuples: list, annotDict: dict) -> list:
    """chooses which tuple to return from a list
    Arguments:
        tuples {list} -- a list of tuples with gene and start position info
        annotDict {dict} -- a dictionary connecting genes to annotation scores
    Returns:
        str -- the gene name of the chosen match
        int -- the start position of the chosen match
    """
    highestScore = 0
    high2Score = 0
    highest = None
    #print(f"third {annotDict["MUG2"]}")
    for tup in tuples:
        if not tup[2] in annotDict:
            continue
        if int(annotDict[tup[2]]) > highestScore:  # If the current match is has the highest score, set it
            highestScore = int(annotDict[tup[2]])
            highest = tup
        elif int(annotDict[tup[2]]) > high2Score:  # If the current match is 2nd highest score, set it
            high2Score = int(annotDict[tup[2]])
    if highestScore > high2Score:
        return highest
    else:  # If there are tied high scores for annotation
        sharedPeps = list()
        sharedPeps.append(highest)
        for tup in tuples:  # Cycles through the matches and chooses the ones with the highest annotation score
            if not tup[2] in annotDict:
                continue
            elif int(annotDict[tup[2]]) == highestScore and tup[2] != highest[2]:
                sharedPeps.append(tup)
        if len(sharedPeps) == 1:
            return sharedPeps[0]
        return sharedPeps
        ### Here is the paralog code - Not running currently ###
        '''with open("Paralogs_rat_UniprotIDs_genes.csv", 'r') as par:
            paralogDict = dict()
            par.readline()
            for line in par:  # Make the paralog dictionary so we can check whether ties are paralogs
                if line.split(',')[0] in paralogDict:
                    paralogDict[line.split(',')[0]].append(line.strip().split(',')[1])
                else:
                    paralogDict[line.split(',')[0]] = [line.strip().split(',')[1]]
                if line.split(',')[1] in paralogDict:
                    paralogDict[line.split(',')[1]].append(line.strip().split(',')[0])
                else:
                    paralogDict[line.split(',')[1]] = [line.strip().split(',')[0]]
            paralogs = list()
            paralogs.append(highest)
            for tup in tuples:  # Cycles through the matches and chooses the ones with the highest annotation score if they are paralogs
                if not tup[0] in annotDict:
                    continue
                elif int(annotDict[tup[0]]) == highestScore and tup[0] != highest[0]:
                    if tup[0] in paralogDict and highest[0] in paralogDict[tup[0]]:
                        paralogs.append(tup)
            if len(paralogs) > 1:
                return paralogs
            else:
                return None, None'''
'''
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
    sky_filename = search_engine_filepath #__prompt_file()
    mod_dict = create_mod_dict(sky_filename, user_PTMs, cleave_rule)
    #PTMs = user_PTMs#__prompt_PTMs(mod_list)
    #Collecting the PTM analysis list from the user
    # NOTE: is it possible to not use PTMs for this analysis?
    #cleave_rule, maxMissed = __promptCleavenMissed()
    #pep_dict = load_pepdict(cleave_rule, maxMissed)  # TODO: make it find the max missed from file
    #annotDict = __load_annotDict("uniprot_rat_annotScore_ids.csv")
    #annotDict = __load_annotDict("uniprot-proteome_human_annotScore.tsv")
    #gene_results = []
    #dup_check_set = set()
    #threshold = 0.01 #What should this be?
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
                print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the generic output file."))
                print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                return -4
    #UNPID = -1
    #try:
    #    UNPID = header.index("Protein Accession")
    #except ValueError:
    #    print("\033[91m {}\033[00m".format("No uniprot ID column detected in the Skyline output file."))
    #    print("Proceeding regardless.\n")
    intensity_values = [i for i, h in enumerate(header) if "intensit" in h.lower()]  # NOTE: This makes the assumption that the intensity columns are directly after the charge column
    if len(intensity_values) == 0:
        intensity_values = -1
    #psm_contributions = defaultdict(int)
    #unmatchedPeps = 0
    #missingPTM = 0
    #totalSeqs = 0
    #unmatchedSequences = []
    #intensity_dict = dict()
    return sequence, pms, mod_dict, intensity_values, sky_filename, None
#EOF

#/home/rob/Documents/Held_Lab/191010_Rollup/SkylineOutput_GloCys055_EXP0107u_RED.csv