#!/usr/bin/env python

"""parse_Skyline.py: parses and rolls up the results from Skyline"""

from collections import defaultdict
import csv
import os.path
from re import findall, finditer, match
import sys
from time import sleep

from proteoSushi_constants import cleave_rules
from ps_utilities import load_pepdict
'''
def __prompt_file() -> str:
    """Prompts the user for the file name

    Returns:
        str -- the filename
    """
    print("\033[96m {}\033[00m".format("Please provide the filename of the Skyline file"))
    print("\033[96m {}\033[00m".format("For example: \nC:/experiment1/Combined/EXP01.csv, or"))
    sky_filename = input("\033[96m {}\033[00m".format("/home/[user]/Documents/experiment1/EXP01.csv\n"))
    if not os.path.exists(sky_filename):
        print("\033[91m {}\033[00m".format("Invalid path, try again!\n"))
        sleep(.5)
        return __prompt_file()
    return sky_filename
'''
def parse_file(sky_filename: str) -> dict:
    """Parses the generic output to retrieve the modification info
    Arguments:
        sky_filename {str} -- the name of the generic output file
    Returns:
        dict -- a dictionary of sequences with the type and location of modifications
        list -- the list of mods in the file
    """
    mod_dict = {}
    with open(sky_filename, 'r') as sky_output:
        tsv_reader = csv.reader(sky_output, quotechar='"')
        header = next(tsv_reader)
        seq_index = header.index("Peptide Sequence")
        mod_index = -1
        try:
            mod_index = header.index("Peptide Modified Sequence") #NOTE: apparently, this changes, so check here.
        except ValueError:
            try:
                mod_index = header.index("Modified Sequence")
            except ValueError:
                print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the Skyline output file."))
                print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
                sys.exit()
        
        mod_list = []
        for row in tsv_reader:
            # This grabs the PTMs in each sequence and builds a list of all of them
            sequence = row[seq_index].replace("L","I")
            mod_seq = row[mod_index]
            mods = findall(r"(\w?\[.+?\])|(\w?\(.+?\))", mod_seq)
            for mod in mods[0]:
                if mod == '':
                    continue
                stripped_mod = mod#.strip('(').strip(')').strip('[').strip(']')
                if not stripped_mod in mod_list:
                    mod_list.append(stripped_mod)

            breaks = finditer(r"(\w?\[.+?\])|(\w?\(.+?\))", mod_seq)
            cut_sites = []
            for breakp in breaks:
                cut_sites.append(breakp.start())
            # This fixes the indices of each mod to connect with the unmodified sequence
            correction = 0
            fixed_indices = []
            num_of_mods = 0
            for site in cut_sites:
                fixed_indices.append(site - correction - 1)  # TODO: Check this!
                correction += len(mods[num_of_mods])
                num_of_mods += 1
            
            i = 0
            while i < len(cut_sites):
                try:
                    mod_dict[sequence].append(tuple((mods[i], fixed_indices[i])))  # TODO: Check This!
                except KeyError:
                    mod_dict[sequence] = [tuple((mods[i], fixed_indices[i]))]
                i += 1
    return mod_dict, mod_list
'''
def __promptCleavenMissed() -> list:
    """Prompts the user for the cleavage rule and max missed cleavages
    Returns:
        str -- the filename
        int -- number of allowed misses
    """
    cleave_rules = {
        # C terminal proteases
        'trypsin/p': (r'[RK]', 'c'),
        'trypsin!p': (r'[RK](?!P)', 'c'),
        'lys-c': (r'[K]', 'c'),
        # N terminal proteases
        'asp-n': (r'[D]', 'n'),
        'asp-nc': (r'[DC]', 'n'),
        'lys-n': (r'[K]', 'n')
    }
    print("\033[96m {}\033[00m".format("Please provide the cleavage rule."))
    cleaveRule = input("\033[96m {}\033[00m".format("Examples include 'trypsin/p', 'trypsin!p', 'lys-c', 'asp-n', asp-nc', 'lys-n'\n"))
    if not cleaveRule in cleave_rules:
        print("\033[91m {}\033[00m".format("Not an accepted format; please try again.\n"))
        sleep(.5)
        return __promptCleavenMissed()
    maxMissed = __promptMissed()
    return cleaveRule, maxMissed

def __promptMissed() -> int:
    maxMissed = None
    try:
        maxMissed = int(input("\033[96m {}\033[00m".format("Please provide the number of allowed missed cleavages.\n")))
    except Exception:
        print("\033[91m {}\033[00m".format("Not an accepted format; please try again.\n"))
        sleep(.5)
        return __promptMissed()
    if maxMissed < 0:
        print("\033[91m {}\033[00m".format("Max missed cleavages must be a positive number; please try again.\n"))
        sleep(.5)
        return __promptMissed()
    return maxMissed
'''
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


def compile_data(search_engine_filepath: str, user_PTMs: list) -> list:
    """Takes the lists and dictionaries needed to parse files

    Arguments:
        search_engine_filepath {str} -- the filepath for the generic search engine output
        user_PTMs {list} -- the list of PTMs that the user will use for the analysis
    Returns:
        int -- index of unmodified sequence
        int -- index of modified sequence
        dict -- pep_dict
        dict -- mod_dict
        int -- intensity_start
        file -- mito_genes (file)
        str -- sky_filename
        dict -- mascot var mod dict (blank)
    """
    sky_filename = search_engine_filepath #__prompt_file()
    mod_dict, mod_list = parse_file(sky_filename)
    PTMs = user_PTMs#__prompt_PTMs(mod_list)
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
    sequence = header.index("Peptide Sequence")
    try:
        pms = header.index("Peptide Modified Sequence") #NOTE: apparently, this changes, so check here.
    except ValueError:
        try:
            pms = header.index("Modified Sequence")
        except ValueError:
            print("\033[91m {}\033[00m".format("No peptide modified sequence column detected in the Skyline output file."))
            print("\033[91m {}\033[00m".format("Please add or modify header with the name \"Peptide Modified Sequence\"\n"))
            sys.exit()
    #UNPID = -1
    #try:
    #    UNPID = header.index("Protein Accession")
    #except ValueError:
    #    print("\033[91m {}\033[00m".format("No uniprot ID column detected in the Skyline output file."))
    #    print("Proceeding regardless.\n")
    intensity_start = header.index("Precursor Charge") + 1  # NOTE: This makes the assumption that the intensity columns are directly after the charge column
    #psm_contributions = defaultdict(int)
    #unmatchedPeps = 0
    #missingPTM = 0
    #totalSeqs = 0
    #unmatchedSequences = []
    #intensity_dict = dict()
    return sequence, pms, mod_dict, intensity_start, sky_filename, None


    '''
    for row in tsv_reader:
        totalSeqs += 1
        raw_seq = row[sequence]
        pep_mod_seq = row[pms]
        #if pep_mod_seq == "LGYILTC[+57]PSNLGTGLR":
        #    print(f"found {pep_mod_seq}")
        pep_seq = row[sequence].replace("L","I")
        if pep_seq is None:
            break
        genes_positions = pep_dict.get(pep_seq)
        if genes_positions and len(genes_positions) == 1:
            gene, start_pos, unpid = list(genes_positions)[0]
            if not pep_seq in mod_dict:
                print("\033[91m {}\033[00m".format(f"{pep_seq} not in mod_dict!"))
                missingPTM += 1
                continue
            mods = mod_dict[pep_seq]
            for mod in list(set(mods)):
                site = start_pos + mod[1]
                intensity_dict, to_add = __addIntensity(intensity_dict, pep_seq, pep_mod_seq, gene, site, row[intensityStart:], cleave_rules[cleave_rule])
                #gene_true_pos = f"{gene}|{start_pos+mod[1]}"
                if to_add:  # TODO: make it check if this is a target gene: if so, put it in the target gene column
                    if gene.upper() + '\n' in mitoGenes:  # TODO: fetch the uniprot ID from the match
                        gene_results.append([gene, site, "", gene, raw_seq, pep_mod_seq, unpid])
                    else:
                        gene_results.append([gene, site, "", "", raw_seq, pep_mod_seq, unpid])  # NOTE: I changed this from a tuple, so things might be different
        elif genes_positions and len(genes_positions) > 1:
            isTarget, match = __chooseHit(genes_positions, mitoGenes, annotDict)
            if match is None:
                print(row)
                continue
            if not match[0] is None:  # Checks to see that there has been a match
                if isinstance(match, tuple):  # Checks if this is just a single match.
                    if not pep_seq in mod_dict:
                        print("\033[91m {}\033[00m".format(f"{pep_seq} not in mod_dict!"))
                        missingPTM += 1
                        continue
                    mods = mod_dict[pep_seq]
                    for mod in list(set(mods)):
                        site = match[1] + mod[1]
                        intensity_dict, to_add = __addIntensity(intensity_dict, pep_seq, pep_mod_seq, match[0], site, row[intensityStart:], cleave_rules[cleave_rule])
                        #gene_true_pos = f"{match[0]}|{match[1]+mod[1]}"
                        if to_add:  # NOTE: Wait, what does this do if it is not supposed to add the intensity? Isn't it supposed to make its own row
                            if isTarget:
                                gene_results.append([match[0], site, "", match[0], raw_seq, pep_mod_seq, match[2]])
                            else:
                                gene_results.append([match[0], site, "", "", raw_seq, pep_mod_seq, match[2]])
                else:  # There are multiple matches
                    if not pep_seq in mod_dict:
                        print("\033[91m {}\033[00m".format(f"{pep_seq} not in mod_dict!"))
                        missingPTM += 1
                        continue
                    mods = mod_dict[pep_seq]
                    # NOTE: I am starting to think whether handling modifications like this is the right way.
                    # The mods are associated with the peptide sequence, but I wonder whether that necessarily 
                    # Applies to each of the database matches :/
                    for mod in list(set(mods)):  # TODO: this is where I add in the intensity columns
                        additGenes = list()
                        for tup in match:
                            additGenes.append(tup[0])
                        assert len(additGenes) > 1, "The # of matches should be >1, but isn't"
                        #gene_true_pos = f"{match[0][0]}|{match[0][1]+mod[1]}"
                        site = match[0][1] + mod[1]
                        intensity_dict, to_add = __addIntensity(intensity_dict, pep_seq, pep_mod_seq, match[0][0], site, row[intensityStart:], cleave_rules[cleave_rule])
                        if to_add:
                            if isTarget:
                                gene_results.append([match[0][0], site, ' '.join(additGenes[1:]), ' '.join([i[0] for i in match]), raw_seq, pep_mod_seq, ' '.join([i[2] for i in match])])
                            else:
                                gene_results.append([match[0][0], site, ' '.join(additGenes[1:]), "", raw_seq, pep_mod_seq, ' '.join([i[2] for i in match])])
            else:
                unmatchedPeps += 1
                unmatchedSequences.append(tuple([raw_seq, pep_mod_seq]))
        else:
            unmatchedPeps += 1
            unmatchedSequences.append(tuple([raw_seq, pep_mod_seq]))
    skyFile.close()
    print(f"Unmatched Peptides: {unmatchedPeps}\nMissing PTMs: {missingPTM}\nTotal Peptides: {totalSeqs}")
    with open("SkylineRollup.csv", 'w', newline = '') as w1:
        #new_key = f"{build_seq}|{gene.upper()}|{str(site)}"
        out_writer = csv.writer(w1)
        header2 = ["Gene", "Site", "Additional Genes", "Target Genes", "Peptide Sequence", "Peptide Modified Sequence", "Uniprot Accession ID"] + header[intensityStart:]
        out_writer.writerow(header2)
        #for i in sorted(list(set([i for i in gene_results])), key=lambda r: r[0]):
        for i in sorted(gene_results, key=lambda r: r[0]):
            writable_row = i + intensity_dict[f"{i[5]}|{i[0].upper()}|{i[1]}"]
            out_writer.writerow(writable_row)
    with open("unmatchedSequences.csv", 'w', newline = '') as w2:
        out_writer = csv.writer(w2)
        header = ["Sequence", "Peptide Modified Sequence"]
        out_writer.writerow(header)
        for i in unmatchedSequences:
            writable_row = list(i)
            out_writer.writerow(writable_row)
            '''
#EOF

#/home/rob/Documents/Held_Lab/191010_Rollup/SkylineOutput_GloCys055_EXP0107u_RED.csv