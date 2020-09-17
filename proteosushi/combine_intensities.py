#!/usr/bin/env python

"""combine_intensities.py: takes data from the parser and combines the mods and intensities"""

import csv
import os
import pandas as pd
from re import finditer, match
from time import sleep

from download_uniprot_AS import download_AS_file
import parse_mascot, parse_MaxQuant, parse_generic, sparql
from proteoSushi_constants import cleave_rules
#from ruputilities import load_pepdict, parse_mascot, parse_maxquant_summary
import ps_utilities


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
    missed_cleavages = __promptMissed()
    return cleaveRule, missed_cleavages

def __promptMissed() -> int:
    missed_cleavages = None
    try:
        missed_cleavages = int(input("\033[96m {}\033[00m".format("Please provide the number of allowed missed cleavages.\n")))
    except Exception:
        print("\033[91m {}\033[00m".format("Not an accepted format; please try again.\n"))
        sleep(.5)
        return __promptMissed()
    if missed_cleavages < 0:
        print("\033[91m {}\033[00m".format("Max missed cleavages must be a positive number; please try again.\n"))
        sleep(.5)
        return __promptMissed()
    return missed_cleavages

def __chooseHit(genes_positions: list, target_genes: list, annot_dict: dict, use_target: bool) -> list:
    """chooses which of the matched sequences to use. If there is one target gene, it will be that one.
    If there are more than one non-target, annotation score decides. If there are more than one target, annotation score decides.
    Arguments:
        genes_positions {list} -- a list of tuples with gene and position info
        target_genes {list} -- a list of strings of target gene names
        annotDict {dict} -- a dictionary connecting genes to annotation scores
    Returns:
        bool -- whether the match(es) is(are) target genes
        *and*
        str -- the gene name of the chosen match
        int -- the start position of the chosen match
        *or*
        list -- a list of tuples with the gene name and start position of each match
    """
    target_tups = list()
    nontarget_tups = list()
    for tup in genes_positions:
        if use_target and tup[0].upper() + '\n' in target_genes:
            #print("Mito gene prioritized!")
            target_tups.append(tup)
        else:
            nontarget_tups.append(tup)
    if len(target_tups) == 1:
        return True, target_tups[0]
    elif len(target_tups) > 1:
        return True, __chooseTup(target_tups, annot_dict)
    elif len(nontarget_tups) == 1:
        return False, nontarget_tups[0]
    elif len(nontarget_tups) > 1:  # Only non-target proteins (and more than 1)
        return False, __chooseTup(nontarget_tups, annot_dict)
    assert False, "ERROR: chooseHit had 0 tuples sent in!"
    #return False, None, None

def __chooseTup(tuples: list, annot_dict: dict) -> list:
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
        if not tup[2] in annot_dict:
            continue
        if int(annot_dict[tup[2]]) > highestScore:  # If the current match is has the highest score, set it
            highestScore = int(annot_dict[tup[2]])  # NOTE: This is not working in maxquant, it never gets here
            highest = tup
        elif int(annot_dict[tup[2]]) > high2Score:  # If the current match is 2nd highest score, set it
            high2Score = int(annot_dict[tup[2]])
    if highestScore > high2Score:
        return highest
    else:  # If there are tied high scores for annotation
        no_scores = True
        sharedPeps = list()
        sharedPeps.append(highest)
        for tup in tuples:  # Cycles through the matches and chooses the ones with the highest annotation score
            if not tup[2] in annot_dict:
                continue
            elif int(annot_dict[tup[2]]) == highestScore and tup[2] != highest[2]:
                no_scores = False
                sharedPeps.append(tup)
        if not sharedPeps[0] is None and len(sharedPeps) == 1:  # If there is a top scorer, send that one
            return sharedPeps[0]
        elif no_scores:  # If none of the proteins have a score, send all of them
            return tuples
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

def __add_intensity(intensity_dict: dict, pep_seq: str, pep_mod_seq: str, gene: str, site: int, 
                    intensities: list, rule: str) -> list:
    """Adds the intensity of the site to the dictionary, if entry exists, adds to the numbers
    Arguments:
        intensity_dict {dict} -- a dictionary connecting the pep/gene/site to intensity
        pep_seq {str} -- the unmodified sequence
        pep_mod_seq {str} -- the modified sequence of the peptide
        gene {str} -- the gene name that the peptide belongs to
        site {int} -- the position of the mod in the protein
        intensities {list} -- 
        rule {str} -- 
    Returns:
        dict -- the intensity_dict following the update
        bool -- True if a new entry was added, False if not
    """
    key = f"{pep_mod_seq}|{gene.upper()}|{str(site)}"
    #if gene.upper() == "ACTB":
    #    print("ACTB_14")
    if key in intensity_dict:
        old_intensities = intensity_dict[key]
        new_intensities = list()
        new_ints = list()
        i = 0
        # Combines the new values for the peak
        while i < len(intensities):
            old_int = old_intensities[i]
            if not old_int or "#N/A" == old_int or "NaN" == old_int or not old_int:
                old_int = 0
            new_int = intensities[i]
            if not new_int or "#N/A" == new_int or "NaN" == new_int or not new_int:
                new_int = 0
            new_ints.append(new_int)
            #print(f"values {repr(old_int)}, {repr(new_int)}")
            new_intensities.append(float(old_int) + float(new_int))
            i += 1
        # Adds to the number of combined peptides for averaging later (if needed)
        if any(new_ints):
            new_intensities.append(int(old_intensities[i]) + 1)  # NOTE: what does it mean if this is out of range?
        else:
            new_intensities.append(int(old_intensities[i]))
        # Replaces the combined numbers in the dictionary
        intensity_dict[key] = new_intensities
        return intensity_dict, False
    else:
        breaks = finditer(rule[0], pep_mod_seq)
        cut_sites = []
        cut_peptides = []
        for breakp in breaks:
            cut_sites.append(breakp.start())
        if len(cut_sites) < 1:
            new_intensities = []
            for new_int in intensities:
                if "#N/A" == new_int:
                    new_int = 0
                new_intensities.append(new_int)
            # Adds to the number of combined peptides for averaging later (if needed)
            new_intensities.append(1)
            # Replaces the combined numbers in the dictionary
            intensity_dict[key] = new_intensities
            return intensity_dict, True

        # Find sites
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
        # If the skipped portion has a cysteine, it should have its own entry
        i = 0
        peps_with_C = 0
        while i < len(cut_peptides):
            if 'C' in cut_peptides[i]:
                peps_with_C += 1
            i += 1
        if peps_with_C > 1:
            new_intensities = []
            for new_int in intensities:
                if "#N/A" == new_int:
                    new_int = 0
                new_intensities.append(new_int)
            # Adds in the number of peptides combined (1 so far)
            new_intensities.append(1)
            # Puts the peak sums in the dictionary
            intensity_dict[key] = new_intensities
            return intensity_dict, True

        # Check for subsegments in the dictionary and if there is a match, add to the original
        for pep in cut_peptides:  # TODO: I need to check all of the genes that match and make sure that they are the same, even though the order may be different
            new_key = f"{pep}|{gene.upper()}|{str(site)}"
            if 'C' in pep and new_key in intensity_dict:
                old_intensities = intensity_dict[new_key]
                new_intensities = list()
                i = 0
                while i < len(intensities):
                    old_int = old_intensities[i]
                    if "#N/A" == old_int or "NaN" == old_int or not old_int:
                        old_int = 0
                    new_int = intensities[i]
                    if "#N/A" == new_int or "NaN" == new_int or not new_int:
                        new_int = 0
                    #print(repr(old_int))
                    #print(repr(new_int))
                    new_intensities.append(float(old_int) + float(new_int))
                    i += 1
                # Adds in the number of combined peptide peaks (1 so far)
                new_intensities.append(1)
                # Adds a new entry for these peaks
                intensity_dict[new_key] = new_intensities
                return intensity_dict, False
        
        new_intensities = []
        for new_int in intensities:
            if "#N/A" == new_int:
                new_int = 0
            new_intensities.append(new_int)
        # Adds in the number of combined peptide peaks (1 so far)
        new_intensities.append(1)
        # Adds a new entry for these peaks
        intensity_dict[key] = new_intensities
        return intensity_dict, True

def __load_annot_dict(annot_file: str) -> dict:
    """loads the annotation score file to make a dict that connects gene to annotation score

    Arguments:
        annot_file {str} -- the filename for the annotation score file
    Returns:
        dict -- a dictionary connecting genes to annotation scores
    """
    annot_dict = {}
    #print(annot_file)
    if annot_file == "":
        return annot_dict
    with open(annot_file, 'r') as r1:
        tsv_reader = csv.reader(r1, delimiter='\t')
        header = next(tsv_reader)
        uniprot_ID = header.index("Entry")
        annot_score = header.index("Annotation")
        for row in tsv_reader:
            annot_dict[row[uniprot_ID]] = row[annot_score][0]
    return annot_dict


def parse_output(search_engine: str, search_engine_filepath: str) -> list:
    """Parses the search engine output and grabs info to return to the gui

    Arguments:
        search_engine {str} -- "maxquant", "mascot", or "generic"
        search_engine_filepath {str} -- the filepath for search engine output
    Returns:
        int -- max missed cleavages allowed
        str -- protease used to cleave the proteins
    """
    rollup_file = search_engine
    if rollup_file == "generic":
        missed_cleavages = -1
        enzyme = ""
        PTMs = parse_generic.get_PTMs(search_engine_filepath)
    elif rollup_file == "maxquant":
        MQ_dir = search_engine_filepath
        sum_file = os.path.join(MQ_dir, "summary.txt")
        enzyme, missed_cleavages = ps_utilities.parse_maxquant_summary(sum_file)
        mod_dict, PTMs = parse_MaxQuant.parse_evidence(os.path.join(MQ_dir, "evidence.txt"))  # This grabs the PTMs from the evidence file
    elif rollup_file == "mascot":
        enzyme, quant_range, var_mod_map, missed_cleavages = ps_utilities.parse_mascot(search_engine_filepath)
        PTMs = list(var_mod_map.keys())
    else:
        assert False, "Not a valid file from search"

    return missed_cleavages, enzyme, PTMs


def rollup(search_engine: str, search_engine_filepath: str, use_target_list: bool, 
           target_list_filepath: str, max_missed_cleavages: int, protease: str, 
           fdr_threshold: float, use_quant: bool, user_PTMs: list, 
           proteome_fasta_filepath: str, intensity_method: str, add_annotation: bool, 
           species_id: str):
    """starts proteoSushi rollup when called by run_proteoSushi

    Arguments:
        search_engine {str} -- "maxquant", "mascot", or "generic"
        search_engine_filepath {str} -- the filepath for search engine output
        use_target_list {bool} -- whether to use the target list to prioritize matches
        target_list_filepath {str} -- the filepath for the target gene list
        max_missed_cleavages {int} -- the maximum allowed missed cleavages in a peptide
        protease {str} -- the protease used previously to cleave the proteins in the sample
        fdr_threshold {float} -- threshold used for the pep_expect and PEP columns in mascot/maxquant
        use_quant {bool} -- whether to sum/average the quant values
        user_PTMs {list} -- a list of PTMs that will be used in analysis
        proteome_fasta_filepath {str} -- the filepath of the proteome fasta file
        intensity_method {str} -- whether to "sum" or "average" the peaks for combined peptides
        add_annotation {bool} -- whether to query uniprot and add the annotation onto the results
        species_id {str} -- the species ID (e.g. 9606) to get the annotation score dictionary from Uniprot
    """
    threshold = fdr_threshold
    use_target = use_target_list
    use_intensities = use_quant
    if species_id != "":
        annot_filename = download_AS_file(species_id)
    else:
        annot_filename = ""
    rollup_file = search_engine

    if rollup_file == "generic":
        output_filename = "generic_rollup.csv"
        sequence_index, modified_sequence_index, mod_dict, intensity_start, data_filename, \
            var_mod_dict = parse_generic.compile_data(search_engine_filepath, user_PTMs)
        data_file = open(data_filename, 'r')
        tsv_reader = csv.reader(data_file, quotechar='"')
    elif rollup_file == "maxquant":
        output_filename = "maxquant_rollup.csv"
        sequence_index, modified_sequence_index, mod_dict, intensity_start, data_filename, \
            var_mod_dict = parse_MaxQuant.compile_data(search_engine_filepath, user_PTMs)
        data_file = open(data_filename, 'r')
        tsv_reader = csv.reader(data_file, delimiter='\t', quotechar='"')
    elif rollup_file == "mascot":
        output_filename = "mascot_rollup.csv"
        sequence_index, modified_sequence_index, mod_dict, intensity_start, data_filename, \
            var_mod_dict = parse_mascot.compile_data(search_engine_filepath, user_PTMs)
        data_file = open(data_filename, 'r')
        tsv_reader = csv.reader(data_file, quotechar='"')
    else:
        assert False, "Not a valid file from search"

    enzyme = protease
    missed_cleavages = max_missed_cleavages
    pep_dict = ps_utilities.load_pepdict(proteome_fasta_filepath, enzyme, missed_cleavages)
    annotDict = __load_annot_dict(annot_filename)

    if use_target:
        target_genes = open(target_list_filepath, 'r').readlines()
    else:
        target_genes = []

    header = next(tsv_reader)
    # Use the PEP column if the user chose Maxquant output
    if rollup_file == "maxquant":
        false_disc = header.index("PEP")
    intensity_dict = dict()
    gene_results = list()
    sparql_input = list()
    unmatched_peps = 0
    missing_PTM = 0
    total_seqs = 0
    over_threshold = 0
    unmatched_sequences = []

    for row in tsv_reader:
        # If a maxquant file is used and the false disc rate of this peptide is not below the threshold, skip it
        if rollup_file == "maxquant" and not threshold is None and float(row[false_disc]) > threshold:
            over_threshold += 1
            continue
        total_seqs += 1
        raw_seq = row[sequence_index]
        if raw_seq == "FACAVVCIQK":
            print("start")
        pep_mod_seq = row[modified_sequence_index]
        pep_seq = row[sequence_index].replace("L","I")
        if pep_seq is None:
            print("HALT!")  # A handled exception would be preferable here
            break
        genes_positions = pep_dict.get(pep_seq)
        if use_intensities and rollup_file == "generic":
            intensities = [e for i, e in enumerate(row) if i in intensity_start]
            if intensities == '' or intensities[0] == '':
                continue
            intensity_header = [e for i, e in enumerate(header) if i in intensity_start]
        elif rollup_file == "maxquant":
            pep_mod_seq = pep_mod_seq.strip('_')

            if use_intensities:
                #intensities = row[intensity_start]
                intensities = [e for i, e in enumerate(row) if i in intensity_start]
                if intensities == '' or intensities[0] == '':
                    continue
                #intensities = [intensities]
                intensity_header = [e for i, e in enumerate(header) if i in intensity_start]
        elif rollup_file == "mascot":
            if pep_mod_seq == "":
                continue
            new_seq = pep_seq
            pep_mod_seq = pep_mod_seq.split('.')[1]
            i = len(pep_mod_seq) - 1
            inv_mod_dict = {v:k for k, v in var_mod_dict.items()}
            while i >= 0:
                if pep_mod_seq[i] != "0":  # If there is a mod
                    new_seq = new_seq[:i+1] + '(' + inv_mod_dict[pep_mod_seq[i]] + ')' + new_seq[i+1:]
                i -= 1
            pep_mod_seq = new_seq

            if use_intensities:
                if len(row) < len(header):  # If the row is cut short (of the intensity cells), skip to the next one
                    continue
                while intensity_start < len(header) and '/' in row[intensity_start]:
                    intensity_start += 2  # TODO: See if you can send this as an error to the GUI
                assert intensity_start < len(header), "Mascot intensity values may be missing, please check"

                intensities = row[intensity_start+1::2]
                if intensities[0] == "---":  # If there are no intensity values in this line, go to the next one
                    continue
                intensity_header = row[intensity_start::2]  # Ideally, this would happen outside of this loop

        if genes_positions and len(genes_positions) == 1:
            gene, start_pos, unpid = list(genes_positions)[0]
            if not pep_mod_seq in mod_dict:
                print("\033[91m {}\033[00m".format(f"{pep_seq} not in modDict!"))
                missing_PTM += 1
                continue
            mods = mod_dict[pep_mod_seq]
            if not (any(mods) and set([m[0] for m in mods]) & set(user_PTMs)): 
                continue
            for mod in list(set(mods)):
                if not mod[0] in user_PTMs:
                    continue
                site = start_pos + mod[1]
                #if gene.upper() == "ACTL6A":
                #    print("ACTL6A")
                if use_intensities:
                    intensity_dict, to_add = __add_intensity(intensity_dict, pep_seq, pep_mod_seq, 
                                                             gene, site, intensities, cleave_rules[enzyme])
                else:
                    key = f"{pep_mod_seq}|{gene.upper()}|{str(site)}"
                    if key in intensity_dict:
                        to_add = False
                    else:
                        to_add = True
                        intensity_dict[key] = 0
                if to_add:
                    if use_target and gene.upper() + '\n' in target_genes:  # TODO: fetch the uniprot ID from the match
                        gene_results.append([gene, site, "", gene, raw_seq, pep_mod_seq, unpid])
                        if len(unpid) >= 5:
                            sparql_input.append(tuple((unpid, site, gene)))
                    else:
                        gene_results.append([gene, site, "", "", raw_seq, pep_mod_seq, unpid])  # NOTE: I changed this from a tuple, so things might be different
                        if len(unpid) >= 5:
                            sparql_input.append(tuple((unpid, site, gene)))
        elif genes_positions and len(genes_positions) > 1:
            isTarget, match = __chooseHit(genes_positions, target_genes, annotDict, use_target)
            # If there was > 1 target genes, non-target genes, or a combination, AND none was chosen
            if match is None:  
                continue
            if not match[0] is None:  # Checks to see that there has been a match
                if isinstance(match, tuple):  # Checks if this is just a single match.
                    gene = match[0]
                    if not pep_mod_seq in mod_dict:
                        print("\033[91m {}\033[00m".format(f"{pep_seq} not in modDict!"))
                        missing_PTM += 1
                        continue
                    mods = mod_dict[pep_mod_seq]
                    for mod in list(set(mods)):
                        if not mod[0] in user_PTMs:
                            continue
                        site = match[1] + mod[1]
                        #if gene.upper() == "ACTL6A":
                        #    print("ACTL6A")
                        to_add = None
                        if use_intensities:
                            intensity_dict, to_add = __add_intensity(intensity_dict, pep_seq, pep_mod_seq, match[0], site, intensities, cleave_rules[enzyme])
                        else:
                            key = f"{pep_mod_seq}|{gene.upper()}|{str(site)}"
                            if key in intensity_dict:
                                to_add = False
                            else:
                                to_add = True
                                intensity_dict[key] = 0
                        if to_add:  # NOTE: Wait, what does this do if it is not supposed to add the intensity? Isn't it supposed to make its own row
                            if isTarget:
                                gene_results.append([gene, site, "", gene, raw_seq, pep_mod_seq, match[2]])
                                if len(match[2]) >= 5:
                                    sparql_input.append(tuple((match[2], site, gene)))
                            else:
                                gene_results.append([gene, site, "", "", raw_seq, pep_mod_seq, match[2]])
                                if len(match[2]) >= 5:
                                    sparql_input.append(tuple((match[2], site, gene)))
                else:  # There are multiple matches
                    if not pep_mod_seq in mod_dict:
                        print("\033[91m {}\033[00m".format(f"{pep_seq} not in modDict!"))
                        missing_PTM += 1
                        continue
                    mods = mod_dict[pep_mod_seq]
                    for mod in list(set(mods)):
                        if not mod[0] in user_PTMs:
                            continue
                        additGenes = list()
                        for tup in match:
                            additGenes.append(tup[0])
                        assert len(additGenes) > 1, "The # of matches should be >1, but isn't"
                        additGenes = list(set(additGenes))
                        #gene_true_pos = f"{match[0][0]}|{match[0][1]+mod[1]}"
                        site = match[0][1] + mod[1]
                        if use_intensities:
                            intensity_dict, to_add = __add_intensity(intensity_dict, pep_seq, pep_mod_seq, match[0][0], site, intensities, cleave_rules[enzyme])
                        else:
                            if any(f"{pep_mod_seq}|{gene.upper()}|{str(site)}" in intensity_dict for gene in additGenes):
                                to_add = False
                            else:
                                to_add = True
                                for gene in additGenes:
                                    intensity_dict[f"{pep_mod_seq}|{gene.upper()}|{str(site)}"] = 0
                        if to_add:
                            additGenes.remove(match[0][0])
                            if isTarget:
                                gene_results.append([match[0][0], site, ' '.join(additGenes), ' '.join([i[0] for i in match]), raw_seq, pep_mod_seq, ' '.join([i[2] for i in match])])
                                if len(match[0][2]) >= 5:
                                    sparql_input.append(tuple((match[0][2], site, match[0][0])))
                            else:
                                gene_results.append([match[0][0], site, ' '.join(additGenes), "", raw_seq, pep_mod_seq, ' '.join([i[2] for i in match])])
                                if len(match[0][2]) >= 5:
                                    sparql_input.append(tuple((match[0][2], site, match[0][0])))
            else:
                unmatched_peps += 1
                unmatched_sequences.append(tuple([raw_seq, pep_mod_seq]))
        else:
            unmatched_peps += 1
            unmatched_sequences.append(tuple([raw_seq, pep_mod_seq]))
    data_file.close()
    # Prints the stats from the rollup
    print(f"Unmatched Peptides: {unmatched_peps}\nMissing PTMs: {missing_PTM}\nTotal Peptides: {total_seqs}")
    # TODO: Eventually make it so this file isn't made in the first place
    #os.remove("pepdict.pepdict")

    ########################
    #Just to annotate EGFR #
    ########################
    if add_annotation and False:
        sparql_dict = dict()
        #with open("EGFR_annotation_uniprot_unfixed.csv", 'r') as annotation_file:
        '''
        with open("sparql-complete.csv", 'r') as annotation_file:
            annotations_full = pd.DataFrame(columns=annotation_file.readline().split(','))
            #print("here")
            for line in annotation_file:
                if line[:5] == "entry":
                    continue
                else:
                    for line_split in csv.reader([line], delimiter=',', quotechar='"'):
                        annotations_full = annotations_full.append(pd.Series(line_split, index=annotations_full.columns), ignore_index=True)
        '''
        annotations_full = pd.read_csv("sparql-complete.csv")
        #annotations_full = pd.read_json("sparql-complete.srj")
        annotations_full = annotations_full[annotations_full.entry != "entry"]
        sparql_output, sparql_dict = sparql.process_sparql_output(annotations_full, sparql_dict)

    # If the user chose, it combines the annotation onto the rollup results (eventually)
    if add_annotation and True:
        print("Querying Uniprot for Annotations!")
        batch = 5
        i = 0
        results_annotated = 0
        sparql_output_list = list()
        sparql_input = sorted(list(set(sparql_input)), key=lambda x: x[2])
        sparql_dict = dict()
        # Separates the input into batches then sends those batches
        while i + batch <= len(sparql_input):
            # Makes the request and sends it to uniprot
            batch_output = sparql.sparql_request(sparql_input[i:i+batch])
            # If after all attempts to get annotations for this batch has failed, this is reported and the next batch will be sent
            if batch_output is None or (not isinstance(batch_output, str) and batch_output.empty):
                print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+batch+1} not annotated!"))
                i += batch
                continue
            # This processes and combines the annotations to 1 per site
            sparql_output, sparql_dict = sparql.process_sparql_output(batch_output, sparql_dict)
            if not sparql_output:
                print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+batch+1} not annotated!"))
                i += batch
                continue
            sparql_output_list.append(sparql_output)
            i += batch
            results_annotated += batch
            print("\033[96m {}\033[00m" .format(f"{round(float(results_annotated)/len(sparql_input)*100, 2)}% of results annotated"))
        batch_output = sparql.sparql_request(sparql_input[i:])
        if batch_output is None or (not isinstance(batch_output, str) and batch_output.empty):
            print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+len(sparql_input[i:])+1} not annotated!"))
            pass
        else:
            sparql_output, sparql_dict = sparql.process_sparql_output(batch_output, sparql_dict)
            if not sparql_output:
                print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+batch+1} not annotated!"))
                pass
            else:
                sparql_output_list.append(sparql_output)
                results_annotated += len(sparql_input[i:])
        print("\033[96m {}\033[00m" .format(f"{round(float(results_annotated)/len(sparql_input)*100, 2)}% of results annotated"))
        # Writes the annotations to a separate file in case the user wants to view those in a more vertical way
        with open("sparql_annotations.csv", 'w') as spql_annot:
            out_writer = csv.writer(spql_annot)
            for line in sparql_output_list:
                out_writer.writerow(line)
    print("Writing the rollup output file")

    # Prints out the completed rollup with annotations from Uniprot (if requested)
    with open(output_filename, 'w', newline = '') as w1:
        out_writer = csv.writer(w1)
        header2 = ["Gene", "Site", "Shared Genes", "Target Genes", "Peptide Sequence", 
            "Peptide Modified Sequence", "Uniprot Accession ID"]
        if use_intensities:
            header2 += intensity_header
        if add_annotation:
            annotations_length = max(len(v) for k, v in sparql_dict.items())
            if annotations_length > 0:
                header2 += ["entry", "position", "lengthOfSequence"]
                annotations_length -= 3
            while annotations_length > 0:
                header2 += ["begin", "end", "regionOfInterest", "catalyicActivity", "location", 
                            "ec", "rhea", "type", "comment"]
                annotations_length -= 12
        out_writer.writerow(header2)
        # This builds the rollup output file depending on what the user chose
        for i in sorted(gene_results, key=lambda r: r[0]):
            if rollup_file == "maxquant":
                if not any(ptm[:2].lower() in i[5] for ptm in user_PTMs):
                    continue
            elif not any(ptm in i[5] for ptm in user_PTMs):  # If none of the chosen ptms are in the rollup line
                continue
            # Start by adding in the base data
            writable_row = i
            if use_intensities:  # Add to that the intensity data if requested
                if intensity_method == "sum":  # Reports the sum of each peak
                    writable_row += intensity_dict[f"{i[5]}|{i[0].upper()}|{i[1]}"][:-1]
                elif intensity_method == "average":  # Calculates the average for each peak and reports
                    N = intensity_dict[f"{i[5]}|{i[0].upper()}|{i[1]}"][-1]
                    intensities = intensity_dict[f"{i[5]}|{i[0].upper()}|{i[1]}"][:-1]
                    writable_row += [float(x)/N for x in intensities]
            if add_annotation:
                try:
                    writable_row += sparql_dict[i[6] + '|' + str(i[1])]
                except KeyError:
                    #print(i[6] + '|' + str(i[1]) + " not in dict")
                    pass
            out_writer.writerow(writable_row)
        #print(sparql_dict)

    # Puts all of the unmatched sequences into a new file
    with open(f"unmatched_sequences_{rollup_file}", 'w', newline = '') as w2:
        out_writer = csv.writer(w2)
        header = ["Sequence", "Peptide Modified Sequence"]
        out_writer.writerow(header)
        for i in unmatched_sequences:
            writable_row = list(i)
            out_writer.writerow(writable_row)
#EOF