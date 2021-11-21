#!/usr/bin/env python

"""combine_intensities.py: takes data from the parser and combines the mods and intensities"""

import csv
import logging
import os
from re import findall

try:
    from .download_uniprot_AS import download_AS_file
    from .parse_mascot import compile_data_mascot
    from .parse_MaxQuant import compile_data_maxquant, compile_localization_data_maxquant, get_MQ_PTMs
    from .parse_generic import compile_data_generic, compile_localization_data_generic, get_PTMs
    from .sparql import process_sparql_output, sparql_request
    from .proteoSushi_constants import cleave_rules, annotation_type_dict, secondary_annotations
    from .ps_utilities import clean_pep_seq, parse_mascot, load_pepdict, parse_maxquant_summary
except ImportError:
    from download_uniprot_AS import download_AS_file
    from parse_mascot import compile_data_mascot
    from parse_MaxQuant import compile_data_maxquant, compile_localization_data_maxquant, get_MQ_PTMs
    from parse_generic import compile_data_generic, compile_localization_data_generic, get_PTMs
    from sparql import process_sparql_output, sparql_request
    from proteoSushi_constants import cleave_rules, annotation_type_dict, secondary_annotations
    from ps_utilities import clean_pep_seq, parse_mascot, load_pepdict, parse_maxquant_summary


def __chooseHit(genes_positions: list, target_genes: list, annot_dict: dict, use_target: bool) -> list:
    """Chooses which of the matched sequences to use. If there is one target gene, it will be that one.
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



def consolidate_sequence(new_pep_mod_seq: str, new_user_PTMs: list) -> str:
    """Changes the modified peptide sequence to only have the relevant PTMs for indexing

    Arguments:
        new_pep_mod_seq {str} -- the modified peptide sequence
        new_user_PTMs {list} -- the PTMs selected by the user
    Returns:
        str -- the modified peptide sequence with irrelevant PTMs removed (for indexing)
    """
    try:
        mods = findall(r"(\[.+?\])|(\(.+?\(?.\)?\))", new_pep_mod_seq).remove('')
    except ValueError:
        mods = findall(r"(\[.+?\])|(\(.+?\(?.\)?\))", new_pep_mod_seq)
    # If there is a PTM in the sequence that isn't in the user list
    #if all([mod[0] in new_user_PTMs for mod in mods]):
    #    return new_pep_mod_seq

    for mod in mods:
        if not mod[1][1:-1] in new_user_PTMs:
            new_pep_mod_seq = new_pep_mod_seq.replace(mod[1], "")
    return new_pep_mod_seq

def __add_intensity(intensity_dict: dict, new_pep_mod_seq: str, genes: list, site: int, 
                    intensities: list) -> list:
    """Adds the intensity of the site to the dictionary, if entry exists, adds to the numbers

    Arguments:
        intensity_dict {dict} -- a dictionary connecting the pep/gene/site to intensity
        new_pep_mod_seq {str} -- the consolidated, modified sequence of the peptide
        genes {list} -- a list of gene names associated with this peptide
        site {int} -- the position of the mod in the protein
        intensities {list} -- 
    Returns:
        dict -- the intensity_dict following the update
        bool -- True if a new entry was added, False if not
    """
    to_adds = []
    for gene in genes:
        key = f"{new_pep_mod_seq}|{gene.upper()}|{str(site)}"
        #if gene.upper() == "SQRDL":
        #    print("SQRDL_379")
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
                new_intensities.append(float(old_int) + float(new_int))
                i += 1

            # Adds to the number of combined peptides for averaging later (if needed)
            if any(new_ints):
                new_intensities.append(int(old_intensities[i]) + 1)  # NOTE: what does it mean if this is out of range?
            else:
                new_intensities.append(int(old_intensities[i]))

            # Replaces the combined numbers in the dictionary
            intensity_dict[key] = new_intensities
            to_adds.append(False)

        else:
            new_intensities = []
            for new_int in intensities:
                if "#N/A" == new_int:
                    new_int = 0
                new_intensities.append(new_int)

            # Adds in the number of combined peptide peaks (1 so far)
            new_intensities.append(1)
            # Adds a new entry for these peaks
            intensity_dict[key] = new_intensities
            to_adds.append(True)

    return intensity_dict, any(to_adds)

def __load_annot_dict(annot_file: str) -> dict:
    """loads the annotation score file to make a dict that connects gene to annotation score

    Arguments:
        annot_file {str} -- the filename for the annotation score file
    Returns:
        dict -- a dictionary connecting genes to annotation scores
    """
    annot_dict = {}
    if annot_file == "":
        return annot_dict
    if annot_file == "ERROR: Invalid identifier":
        return None
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
    if search_engine == "generic":
        missed_cleavages = -1
        protease = ""
        PTMs = get_PTMs(search_engine_filepath)
        if PTMs == -3:
            return -3, None, None  # "A sequence in the Peptide Modified Sequence column is missing PTMs"
        if PTMs == -4:
            return -4, None, None
    elif search_engine == "maxquant":
        MQ_dir = search_engine_filepath
        sum_file = os.path.join(MQ_dir, "summary.txt")
        protease, missed_cleavages = parse_maxquant_summary(sum_file)
        PTMs = get_MQ_PTMs(os.path.join(MQ_dir, "evidence.txt"))  # This grabs the PTMs from the evidence file
    elif search_engine == "mascot":
        protease, quant_range, var_mod_map, missed_cleavages = parse_mascot(search_engine_filepath)
        if protease == -5:
            return -5, None, None
        PTMs = list(var_mod_map.keys())
    else:  # It shouldn't be able to get here
        assert False, "Not a valid file from search"

    return missed_cleavages, protease, PTMs


def __compress_annotations(annotation_list: list) -> list:
    """Takes a list of annotations and makes a single combined group

    Arguments:
        annotation_list {list} -- a list of annotations comprising 1 or more groups
    Returns:
        list -- the compressed list of annotations
    """
    #entry,position,lengthOfSequence,location,ec,rhea,type,comment,begin,end,regionOfInterest
    #entry,position,lengthOfSequence,begin,end,regionOfInterest,location,ec,rhea,type,comment(,begin,end,...)
    # NOTE: you may need to update these numbers if you add or delete a column

    begin_index = 3
    end_index = 4
    type_index = 9
    comment_index = 10
    length_uniprot_annotations = 8  # The number of 
    # Grabs the entry, position, and lengthofsequence columns
    new_annotations = annotation_list[:begin_index]
    # Compresses the begin and end columns into a range column
    new_annotations.append(f"{annotation_list[begin_index]}-{annotation_list[end_index]}")
    # Adds the rest of the columns (up to type)
    new_annotations += annotation_list[end_index+1:type_index]
    # Adds a blank spot for each annotation (as listed in proteosushi_constants.py)
    new_annotations += [""]*len(annotation_type_dict)

    try:  # Attempts to put the comment in the appropriate column using the type as reference
        new_annotations[annotation_type_dict[annotation_list[type_index]] + length_uniprot_annotations] = annotation_list[comment_index]
    except KeyError:
        if annotation_list[type_index] == "nan":
            pass
        else:
            new_annotations[annotation_type_dict["Other"] + length_uniprot_annotations] = annotation_list[comment_index]

    # NOTE: you will need to update these numbers if you add or remove a column from query
    new_begin_index = begin_index - 3
    new_end_index = end_index - 3
    range_index = 3
    region_index = 4
    location_index = 5
    ec_index = 6
    rhea_index = 7
    new_type_index = 8
    new_comment_index = 9
    secondary_structure = ""
    i = 1  # Basically which group it is on
    while (i * length_uniprot_annotations) + 3 < len(annotation_list):  # Essentially, we are moving through the annotation list and compressing it to new annotations
        if not f"{annotation_list[new_begin_index+(i*length_uniprot_annotations) + 3]}-{annotation_list[new_end_index+(i*length_uniprot_annotations)+3]}" in new_annotations[range_index].split(','):
            new_annotations[range_index] += f",{annotation_list[new_begin_index+(i*length_uniprot_annotations) + 3]}-{annotation_list[new_end_index+(i*length_uniprot_annotations)+3]}"
        if not f"{annotation_list[region_index-2+(i*length_uniprot_annotations) + 3]}" in new_annotations[region_index].split(','):
            new_annotations[region_index] += f",{annotation_list[region_index-2+(i*length_uniprot_annotations) + 3]}"
        if not f"{annotation_list[location_index-2+(i*length_uniprot_annotations)+3]}" in new_annotations[location_index].split(','):
            new_annotations[location_index] += f",{annotation_list[location_index-2+(i*length_uniprot_annotations)+3]}"
        if not f"{annotation_list[ec_index-2+(i*length_uniprot_annotations)+3]}" in new_annotations[ec_index].split(','):
            new_annotations[ec_index] += f",{annotation_list[ec_index-2+(i*length_uniprot_annotations)+3]}"
        if not f"{annotation_list[rhea_index-2+(i*length_uniprot_annotations)+3]}" in new_annotations[rhea_index].split(','):
            new_annotations[rhea_index] += f",{annotation_list[rhea_index-2+(i*length_uniprot_annotations)+3]}"

        try:  # Attempts to put the comment in the appropriate column using the type as reference
            current_comment = annotation_list[new_comment_index-2+(i*length_uniprot_annotations)+3]
            current_type = annotation_list[new_type_index-2+(i*length_uniprot_annotations)+3]

            if (not current_comment in new_annotations[annotation_type_dict[current_type] + length_uniprot_annotations] and
                not current_type in new_annotations[annotation_type_dict[current_type] + length_uniprot_annotations]):
                if new_annotations[annotation_type_dict[current_type]
                                + length_uniprot_annotations]:
                    if current_comment != "nan":
                        new_annotations[annotation_type_dict[current_type]
                                        + length_uniprot_annotations] += f",{current_comment}"
                else:
                    if current_comment == "nan":
                        new_annotations[annotation_type_dict[current_type]
                                        + length_uniprot_annotations] += f"{current_type.replace('_Annotation', '')}"
                    else:
                        new_annotations[annotation_type_dict[current_type]
                                        + length_uniprot_annotations] += f"{current_comment}"
        except KeyError:
            if current_type == "nan":
                pass
            else:
                # If the current comment is not in the Other cell
                if current_type in secondary_annotations:
                    secondary_structure = current_type.replace("_Annotation", "")
                elif (not current_comment in new_annotations[annotation_type_dict["Other"] + length_uniprot_annotations] and
                    not current_type in new_annotations[annotation_type_dict["Other"] + length_uniprot_annotations]):
                    if new_annotations[annotation_type_dict["Other"] + length_uniprot_annotations]:
                        if current_comment != "nan":
                            new_annotations[annotation_type_dict["Other"] 
                                            + length_uniprot_annotations] = f",{current_type.replace('_Annotation', '')}: {current_comment}"
                    else:
                        if current_comment == "nan":
                            new_annotations[annotation_type_dict["Other"] 
                                            + length_uniprot_annotations] = f"{current_type.replace('_Annotation', '')}"
                        else:
                            new_annotations[annotation_type_dict["Other"] 
                                            + length_uniprot_annotations] = f"{current_type.replace('_Annotation', '')}: {current_comment}"
        i += 1
    return [s.replace("nan", "") for s in new_annotations[:8]] + [secondary_structure] + [s.replace("nan", "") for s in new_annotations[8:]]


def batch_write(batch_results: list, search_engine: str, user_PTMs: list, use_quant: bool,
                intensity_method: str, intensity_dict: dict, add_annotation: bool, sparql_input: list,
                use_target_list: bool):
    """Write a batch to the output to save memory (RAM)
    Arguments:
        batch_results {list} -- a portion of the results to be written
        search_engine {str} -- the string of the search engine ("maxquant", "mascot", "generic")
        user_PTMs {list} -- a list of the user chosen PTMs
        use_quant {bool} -- whether the user chose to use quant
        intensity_method {str} -- whether to sum or average the quant
        intensity_dict {dict} -- a dictionary that connects site info to quant
        add_annotation {bool} -- whether to provide uniprot annotations to the results
        sparql_input {list} -- data needed to get uniprot annotations
        use_target_list {bool} -- whether the user chose to use a gene list to affect matches
    Returns:
        list -- row to write
    """
    writable_rows = list()
    if add_annotation:
        sparql_dict = batch_annotate(sparql_input)
        if isinstance(sparql_dict, int) and sparql_dict == 502:
            return 502
        if isinstance(sparql_dict, int) and sparql_dict == 4:
            return 4
        
    # This builds the rollup output file depending on what the user chose
    for i in sorted(batch_results, key=lambda r: r[0]):
        gene = i[0].upper()
        pos = i[1]
        target_index = 4
        pep_mod_seq = i[6]
        uniprot_id = i[8]

        if search_engine == "maxquant":
            if not any(ptm[:2].lower() in pep_mod_seq for ptm in user_PTMs):
                continue
        elif not any(ptm in pep_mod_seq for ptm in user_PTMs):  # If none of the chosen ptms are in the rollup line
            continue
        # Start by adding in the base data
        writable_row = i
        if use_quant:  # Add to that the intensity data if requested
            if intensity_method == "sum":  # Reports the sum of each peak
                writable_row += intensity_dict[f"{pep_mod_seq}|{gene}|{pos}"][:-1]
            elif intensity_method == "average":  # Calculates the average for each peak and reports
                N = intensity_dict[f"{pep_mod_seq}|{gene}|{pos}"][-1]
                intensities = intensity_dict[f"{pep_mod_seq}|{gene}|{pos}"][:-1]
                writable_row += [float(x)/N for x in intensities]
        if not use_target_list:
            writable_row = writable_row[:target_index] + writable_row[target_index+1:]
        if add_annotation:
            try:
                compressed_annotations = __compress_annotations(sparql_dict[uniprot_id.split(' ')[0] + '|' + str(pos)])
                writable_row += compressed_annotations[2:]
                logging.debug(compressed_annotations)
                #print("\033[96m {}\033[00m" .format(f"{round(float(results_annotated)/len(sparql_input)*100, 2)}% of results annotated"))
            except KeyError:
                pass
        writable_rows.append(writable_row)
    return writable_rows

def batch_annotate(sparql_input: list) -> dict:
    """Annotates the provided rollup results returning a dictionary
    Arguments:
        sparql_input {list} -- a segment of rollup results
    Returns:
        dict -- dictionary connecting rollup to annotation
    """
    # If the user chose, it combines the annotation onto the rollup results (eventually)
    #print("\033[95m {}\033[00m".format("\nQuerying Uniprot for Annotations!"))
    batch = 50
    i = 0
    results_annotated = 0
    sparql_input = sorted(list(set(sparql_input)), key=lambda x: x[2])
    sparql_dict = dict()
    # Separates the input into batches then sends those batches
    while i + batch <= len(sparql_input):
        # Makes the request and sends it to uniprot
        #print(f"\nInput range: {i}-{i+batch},\nInput list = {sparql_input[i:i+batch]}")
        logging.debug(f"Input range: {i}-{i+batch},\nInput list = {sparql_input[i:i+batch]}")
        batch_output = sparql_request(sparql_input[i:i+batch])

        #If there is a 502 error, send that to the GUI to display
        if isinstance(batch_output, int) and batch_output == 502:
            return 502

        # If after all attempts to get annotations for this batch has failed, this is reported and the next batch will be sent
        if batch_output is None or (not isinstance(batch_output, str) and batch_output.empty):
            #print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+batch+1} not annotated!"))
            i += batch
            continue
        logging.debug(f"Columns are:\n{batch_output.columns}\nAnnotations are:\n{batch_output}")
        # This processes and combines the annotations to 1 per site
        #print(f"\nIn batch_annotate uniprot annotation is {type(batch_output)}")
        logging.debug(f"In batch_annotate uniprot annotation is {type(batch_output)}")
        sparql_output, sparql_dict = process_sparql_output(batch_output, sparql_dict)
        #print(sparql_output)
        if isinstance(sparql_output, int) and sparql_output == 4:
            return 4
        if not sparql_output:
            #print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+batch+1} not annotated!"))
            i += batch
            continue
        i += batch
        results_annotated += batch
        #print("\nBatch processed successfully!")
        logging.debug("Batch processed successfully!")
    batch_output = sparql_request(sparql_input[i:])

    #If there is a 502 error, send that to the GUI to display
    if isinstance(batch_output, int) and batch_output == 502:
        return 502

    if batch_output is None or (not isinstance(batch_output, str) and batch_output.empty):
        #print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+len(sparql_input[i:])+1} not annotated!"))
        pass
    else:
        sparql_output, sparql_dict = process_sparql_output(batch_output, sparql_dict)
        if isinstance(sparql_output, int) and sparql_output == 4:
            return 4
        if not sparql_output:
            #print("\033[91m {}\033[00m".format(f"Lines {i+2} to {i+batch+1} not annotated!"))
            pass
        else:
            results_annotated += len(sparql_input[i:])
    return sparql_dict
        

def rollup(search_engine: str, search_engine_filepath: str, use_target_list: bool, 
           target_list_filepath: str, max_missed_cleavages: int, protease: str, 
           fdr_threshold: float, use_quant: bool, user_PTMs: list, 
           proteome_fasta_filepath: str, intensity_method: str, add_annotation: bool, 
           species_id: str, output_filename: str, localization_threshold: float) -> int:
    """starts proteoSushi rollup when called by runsparql_input
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
        localization {float} -- the threshold for maxquant localization probability
    Returns:
        int -- possible error flag
    """
    print("Preparing data for rollup...\n")
    logging.info(f"Search Engine is {search_engine}")
    logging.info(f"Localization threshold is {str(localization_threshold)}")
    if search_engine == "generic":
        if not localization_threshold is None:
            sequence_index, modified_sequence_index, localization_indices, mod_dict, intensity_start, \
                data_filename = compile_localization_data_generic(search_engine_filepath, 
                                                                  user_PTMs, 
                                                                  cleave_rules[protease],
                                                                  localization_threshold)
        else:
            sequence_index, modified_sequence_index, mod_dict, intensity_start, data_filename, \
                var_mod_dict = compile_data_generic(search_engine_filepath, user_PTMs, cleave_rules[protease])
        logging.debug("Reading in data file")
        data_file = open(data_filename, 'r')
        tsv_reader = csv.reader(data_file, quotechar='"')
        if intensity_start is None and use_quant:
            return 2
    elif search_engine == "maxquant":
        if not localization_threshold is None:
            sequence_index, modified_sequence_index, localization_indices, mod_dict, intensity_start, \
                data_filename = compile_localization_data_maxquant(search_engine_filepath, 
                                                                   user_PTMs, 
                                                                   cleave_rules[protease],
                                                                   localization_threshold)#
        else:
            sequence_index, modified_sequence_index, mod_dict, intensity_start, data_filename, \
                var_mod_dict = compile_data_maxquant(search_engine_filepath, user_PTMs, cleave_rules[protease])
        data_file = open(data_filename, 'r')
        tsv_reader = csv.reader(data_file, delimiter='\t', quotechar='"')
    elif search_engine == "mascot":
        sequence_index, modified_sequence_index, mod_dict, intensity_start, data_filename, \
            var_mod_dict = compile_data_mascot(search_engine_filepath, user_PTMs, cleave_rules[protease])
        data_file = open(data_filename, 'r')
        tsv_reader = csv.reader(data_file, quotechar='"')
        if intensity_start is None and use_quant:
            return 2
    else:
        assert False, "Not a valid file from search"


    if species_id != "":
        logging.debug("Downloading AS file")
        annot_filename = download_AS_file(species_id)
    else:
        annot_filename = ""

    logging.debug("Loading in pep_dict")
    pep_dict = load_pepdict(proteome_fasta_filepath, protease, max_missed_cleavages)
    logging.debug("Loading in annot_dict")
    annotDict = __load_annot_dict(annot_filename)

    if use_target_list:
        target_genes = open(target_list_filepath, 'r').readlines()
    else:
        target_genes = []

    header = next(tsv_reader)
    # Use the PEP column if the user chose Maxquant output
    if search_engine == "maxquant":
        false_disc = header.index("PEP")
    intensity_dict = dict()
    gene_results = list()
    sparql_input = list()
    unmatched_peps = 0
    missing_PTM = 0
    total_seqs = 0
    over_threshold = 0
    missing_intensities = 0
    total_sites = 0
    batch_size = 100
    unmatched_sequences = []

    print("Beginning PTM site rollup")
    logging.debug("Beginning PTM site rollup")

    # Prints out the completed rollup with annotations from Uniprot (if requested)
    with open(output_filename, 'w', newline = '') as w1:
        out_writer = csv.writer(w1)
        header2 = ["Gene", "Site", "Protein_Name", "Shared_Genes"]
        if use_target_list:
            header2 += ["Target_Genes"] 
        header2 += ["Peptide_Sequence", "Peptide_Modified_Sequence", "Annotation_Score", 
                    "Uniprot_Accession_ID"]
        if use_quant:
            if search_engine == "mascot":
                header2 += [header[intensity_start] + f" ({intensity_method})"]
            else:
                intensity_header = [e for i, e in enumerate(header) if i in intensity_start]
                header2 += [ih + f" ({intensity_method})" for ih in intensity_header]
        if add_annotation:
            header2 += ["Length_Of_Sequence", "Range_of_Interest", "Region_of_Interest", 
                        "Subcellular_Location", "Enzyme_Class", "rhea", "Secondary_Structure", 
                        "Active_Site_Annotation", "Alternative_Sequence_Annotation", 
                        "Chain_Annotation", "Compositional_Bias_Annotation", 
                        "Disulfide_Bond_Annotation", "Domain_Extent_Annotation", 
                        "Lipidation_Annotation", "Metal_Binding_Annotation", 
                        "Modified_Residue_Annotation", "Motif_Annotation", 
                        "Mutagenesis_Annotation", "Natural_Variant_Annotation", 
                        "NP_Binding_Annotation", "Other", "Region_Annotation", "Repeat_Annotation", 
                        "Topological_Domain_Annotation", "Zinc_Finger_Annotation"]
        out_writer.writerow(header2)
        logging.debug("Header written to output file")

        for row in tsv_reader:
            if len(gene_results) >= batch_size:
                writable_rows = batch_write(gene_results, search_engine, user_PTMs, use_quant, 
                                            intensity_method, intensity_dict, add_annotation, 
                                            sparql_input, use_target_list)
                if isinstance(writable_rows, int) and writable_rows == 502:
                    return 502
                if isinstance(writable_rows, int) and writable_rows == 4:
                    return 4
                for writable_row in writable_rows:
                    out_writer.writerow(writable_row)
                gene_results = list()
                sparql_input = list()
                print("\033[92m {}\033[00m".format(f"\n{total_sites} sites rolled-up and written"), end='')
            # If a maxquant file is used and the false disc rate of this peptide is not below the threshold, skip it
            if search_engine == "maxquant" and not fdr_threshold is None and float(row[false_disc]) > fdr_threshold:
                over_threshold += 1
                continue
            total_seqs += 1
            raw_seq = row[sequence_index]
            if len(raw_seq) < 6:
                continue
            #print(raw_seq)
            pep_mod_seq = row[modified_sequence_index]
            pep_seq = row[sequence_index].replace("L","I")
            if pep_seq is None:
                logging.error(f"{search_engine} file is missing pep_seq")
                print("HALT!")  # A handled exception would be preferable here
                break
            genes_positions = pep_dict.get(pep_seq)
            if use_quant and search_engine == "generic":
                intensities = [e for i, e in enumerate(row) if i in intensity_start]
                # If there is no intensity, skip that site
                if intensities == '' or intensities[0] == '':
                    missing_intensities += 1
                    continue
                intensity_header = [e for i, e in enumerate(header) if i in intensity_start]
            elif search_engine == "maxquant":
                pep_mod_seq = pep_mod_seq.strip('_')

                if use_quant:
                    intensities = [e for i, e in enumerate(row) if i in intensity_start]
                    # If there is no intensity, skip that site
                    if intensities == '' or intensities[0] == '':
                        missing_intensities += 1
                        continue
                    intensity_header = [e for i, e in enumerate(header) if i in intensity_start]
            elif search_engine == "mascot":
                if pep_mod_seq == "":
                    missing_PTM += 1
                    continue
                new_seq = raw_seq
                pep_mod_seq = pep_mod_seq.split('.')[1]
                i = len(pep_mod_seq) - 1
                inv_mod_dict = {v:k for k, v in var_mod_dict.items()}
                while i >= 0:
                    if pep_mod_seq[i] != "0":  # If there is a mod
                        new_seq = new_seq[:i+1] + '(' + inv_mod_dict[pep_mod_seq[i]] + ')' + new_seq[i+1:]
                    i -= 1
                pep_mod_seq = new_seq

                if use_quant:
                    if len(row) < len(header):  # If the row is cut short (of the intensity cells), skip to the next one
                        missing_intensities += 1
                        continue
                    while intensity_start < len(header) and '/' in row[intensity_start]:
                        intensity_start += 2  # 
                    # NOTE: This is now an error to the GUI and this line of code shouldn't run.
                    assert intensity_start < len(header), "Mascot intensity values may be missing, please check"

                    intensities = row[intensity_start+1::2]
                    # If there is no intensity, skip that site
                    if intensities[0] == "---":  # If there are no intensity values in this line, go to the next one
                        missing_intensities += 1
                        continue
                    intensity_header = row[intensity_start::2]  # Ideally, this would happen outside of this loop

            new_user_PTMs = user_PTMs
            if search_engine == "maxquant":
                new_user_PTMs = [ptm.lower()[:2] for ptm in user_PTMs]
            
            new_pep_mod_seq, new_pep_seq, missed_cleave_fix = clean_pep_seq(cleave_rules[protease], pep_mod_seq, new_user_PTMs, raw_seq)
            #input("old " + pep_mod_seq + ", new " + new_pep_mod_seq)
            if genes_positions and len(genes_positions) == 1:
                gene, start_pos, unpid, protein_name = list(genes_positions)[0]
                if not new_pep_mod_seq in mod_dict:
                    #print("\033[91m {}\033[00m".format(f"\n{raw_seq} does not have the PTM(s) selected!"), end='')
                    missing_PTM += 1
                    continue
                mods = mod_dict[new_pep_mod_seq]
                #if new_pep_seq == "ACLNPASPIVK":
                #    input(mods)
                # Skip any peptide sequences without a user-chosen PTM
                if not (any(mods) and set([m[0] for m in mods]) & set(user_PTMs)): 
                    missing_PTM += 1
                    continue
                for mod in list(set(mods)):
                    if not mod[0] in user_PTMs:
                        missing_PTM += 1
                        continue
                    site = start_pos + mod[1] + missed_cleave_fix + 1  # The last +1 is to change from 0-indexing to 1-indexing (like humans use)
                    #input("1: " + unpid + ' ' + str(site))
                    if use_quant:  # If the user chose to combine/average intensities
                        intensity_dict, to_add = __add_intensity(intensity_dict, 
                                                                new_pep_mod_seq, 
                                                                [gene], site, intensities)
                    else:
                        key = f"{new_pep_mod_seq}|{gene.upper()}|{str(site)}"
                        if key in intensity_dict:
                            to_add = False
                        else:
                            to_add = True
                            intensity_dict[key] = 0
                    if to_add:
                        if use_target_list and gene.upper() + '\n' in target_genes:  # TODO: fetch the uniprot ID from the match
                            gene_results.append([
                                gene, 
                                site, 
                                protein_name,
                                "", 
                                gene, 
                                new_pep_seq, 
                                new_pep_mod_seq, 
                                annotDict[unpid] if unpid in annotDict else "", 
                                unpid
                                ])
                            total_sites += 1
                            if len(unpid) >= 5:
                                sparql_input.append(tuple((unpid, site, gene)))
                        else:
                            gene_results.append([
                                gene, 
                                site, 
                                protein_name,
                                "", 
                                "", 
                                new_pep_seq, 
                                new_pep_mod_seq, 
                                annotDict[unpid] if unpid in annotDict else "", 
                                unpid
                                ])  # NOTE: I changed this from a tuple, so things might be different
                            total_sites += 1
                            if len(unpid) >= 5:
                                sparql_input.append(tuple((unpid, site, gene)))
            elif genes_positions and len(genes_positions) > 1:
                isTarget, match = __chooseHit(genes_positions, target_genes, annotDict, use_target_list)
                # If there was > 1 target genes, non-target genes, or a combination, AND none was chosen
                if match is None:
                    unmatched_peps += 1
                    continue
                if not match[0] is None:  # Checks to see that there has been a match
                    if isinstance(match, tuple):  # Checks if this is just a single match.
                        gene = match[0]
                        if not new_pep_mod_seq in mod_dict:
                            #print("\033[91m {}\033[00m".format(f"\n{raw_seq} does not have the PTM(s) selected!"), end='')
                            missing_PTM += 1
                            continue
                        mods = mod_dict[new_pep_mod_seq]
                        for mod in list(set(mods)):
                            if not mod[0] in user_PTMs:
                                continue
                            site = match[1] + mod[1] + missed_cleave_fix + 1  # The last +1 is to change from 0-indexing to 1-indexing (like humans use)
                            #input("2: " + match[2] + ' ' + str(site))
                            to_add = None
                            if use_quant:
                                intensity_dict, to_add = __add_intensity(intensity_dict, 
                                                                        new_pep_mod_seq, 
                                                                        [match[0]], 
                                                                        site, intensities)
                            else:
                                key = f"{new_pep_mod_seq}|{gene.upper()}|{str(site)}"
                                if key in intensity_dict:
                                    to_add = False
                                else:
                                    to_add = True
                                    intensity_dict[key] = 0
                            if to_add:
                                if isTarget:
                                    gene_results.append([
                                        gene, 
                                        site, 
                                        match[3],  # Protein name
                                        "", 
                                        gene, 
                                        new_pep_seq, 
                                        new_pep_mod_seq, 
                                        annotDict[match[2]] if match[2] in annotDict else "", 
                                        match[2]  # UNPID
                                        ])
                                    total_sites += 1
                                    if len(match[2]) >= 5:
                                        sparql_input.append(tuple((match[2], site, gene)))
                                else:
                                    gene_results.append([
                                        gene, 
                                        site, 
                                        match[3],
                                        "", 
                                        "", 
                                        new_pep_seq, 
                                        new_pep_mod_seq, 
                                        annotDict[match[2]] if match[2] in annotDict else "", 
                                        match[2]
                                        ])
                                    total_sites += 1
                                    if len(match[2]) >= 5:
                                        sparql_input.append(tuple((match[2], site, gene)))
                    else:  # There are multiple matches
                        if not new_pep_mod_seq in mod_dict:
                            #print("\033[91m {}\033[00m".format(f"\n{raw_seq} does not have the PTM(s) selected!"), end='')
                            missing_PTM += 1
                            continue
                        mods = mod_dict[new_pep_mod_seq]
                        for mod in list(set(mods)):
                            if not mod[0] in user_PTMs:
                                continue
                            additGenes = list()
                            for tup in match:
                                additGenes.append(tup[0])
                            assert len(additGenes) > 1, "The # of matches should be >1, but isn't"
                            additGenes = list(set(additGenes))
                            site = match[0][1] + mod[1] + missed_cleave_fix + 1  # The last +1 is to change from 0-indexing to 1-indexing (like humans use)
                            #input("3: " + match[0][2] + ' ' + str(site))
                            if use_quant:
                                intensity_dict, to_add = __add_intensity(intensity_dict,
                                                                        new_pep_mod_seq, 
                                                                        additGenes, 
                                                                        site, intensities)
                            else:
                                if any(f"{new_pep_mod_seq}|{gene.upper()}|{str(site)}" in intensity_dict for gene in additGenes):
                                    to_add = False
                                else:
                                    to_add = True
                                    for gene in additGenes:
                                        intensity_dict[f"{new_pep_mod_seq}|{gene.upper()}|{str(site)}"] = 0
                            if to_add:
                                for addit_gene in additGenes:
                                    gene_list = [x for x in additGenes if x != addit_gene]
                                    current_match = [x for x in match if x[0] == addit_gene][0]
                                    if isTarget:
                                        gene_results.append([
                                            addit_gene,
                                            site,
                                            current_match[3],
                                            ' '.join(gene_list),
                                            ' '.join([i[0] for i in match]),
                                            new_pep_seq,
                                            new_pep_mod_seq,
                                            annotDict[current_match[2]] if current_match[2] in annotDict else "",
                                            ' '.join([i[2] for i in match])
                                            ])
                                        total_sites += 1
                                        if len(current_match[2]) >= 5:
                                            sparql_input.append(tuple((current_match[2], site, current_match[0])))
                                    else:
                                        gene_results.append([
                                            addit_gene, 
                                            site, 
                                            current_match[3], 
                                            ' '.join(gene_list), 
                                            "", 
                                            new_pep_seq, 
                                            new_pep_mod_seq, 
                                            annotDict[current_match[2]] if current_match[2] in annotDict else "", 
                                            ' '.join([i[2] for i in match])
                                            ])
                                        total_sites += 1
                                        if len(current_match[2]) >= 5:
                                            sparql_input.append(tuple((current_match[2], site, current_match[0])))
                else:  # NOTE: it generally shouldn't get here
                    unmatched_peps += 1
                    unmatched_sequences.append(tuple([raw_seq, pep_mod_seq]))
            else:  # There was no match for the peptide in the pepdict
                unmatched_peps += 1
                unmatched_sequences.append(tuple([raw_seq, pep_mod_seq]))
        # Adds in the last few results
        if len(gene_results) >= 1:
            writable_rows = batch_write(gene_results, search_engine, user_PTMs, use_quant, 
                                        intensity_method, intensity_dict, add_annotation, 
                                        sparql_input, use_target_list)
            for writable_row in writable_rows:
                out_writer.writerow(writable_row)
            print("\033[92m {}\033[00m".format(f"\n{total_sites} sites rolled-up and written"), end='')

    data_file.close()
    # Prints the stats from the rollup
    if unmatched_peps > 0:
        print("\033[93m {}\033[00m".format(f"\n\nUnmatched Peptides: {unmatched_peps}"), end='')
    if missing_PTM > 0:
        print("\033[93m {}\033[00m".format(f"\nMissing Selected PTMs: {missing_PTM}"), end='')
    if over_threshold > 0:
        print("\033[93m {}\033[00m".format(f"\nOver the FDR Threshold: {over_threshold}"), end='')
    if missing_intensities > 0:
        print("\033[93m {}\033[00m".format(f"\nMissing Intensity Value(s): {missing_intensities}"), end='')
    print("\033[93m {}\033[00m".format(f"\nTotal Peptides: {total_seqs}"))

    #print("\033[95m {}\033[00m".format("\nWriting the rollup output file"))
    
    # Puts all of the unmatched sequences into a new file
    with open(f"unmatched_sequences_{search_engine}.csv", 'w', newline = '') as w2:
        out_writer = csv.writer(w2)
        header = ["Sequence", "Peptide Modified Sequence"]
        out_writer.writerow(header)
        for i in unmatched_sequences:
            writable_row = list(i)
            out_writer.writerow(writable_row)
    print("\nUnmatched sequences written\n", end='')
    return 0
#EOF