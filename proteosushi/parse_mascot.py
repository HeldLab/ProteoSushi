
import csv
from collections import defaultdict
import os
from re import match
from time import sleep

from ps_utilities import parse_mascot, load_pepdict #most of the parsing actually happens here 


def __promptFile() -> str:
    """Prompts the user for the file name
    Returns:
        str -- the filename
    """
    print("\033[96m {}\033[00m".format("Please provide the filename of the Mascot file"))
    print("\033[96m {}\033[00m".format("For example: \nC:/experiment1/Combined/F011000.csv, or"))
    filename = input("\033[96m {}\033[00m".format("/home/[user]/Documents/experiment1/F011000.csv\n"))
    if not os.path.exists(filename):
        print("\033[91m {}\033[00m".format("Invalid path, try again!\n"))
        sleep(.5)
        return __promptFile()
    return filename

def __parseFile(filename: str, mod_ids: list, var_mod_map: dict) -> dict:
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
                
            mods = row[var_mods]
            if mods:
                n, center, c = [x for x in mods.split(".")]
                for i, aa in enumerate(list(center)): #split mods; note this excludes termini
                    if aa in mod_ids:
                        try:
                            if not tuple((str(aa), i)) in modDict[pep_seq]:
                                modDict[pep_seq].append(tuple((inv_mod_map[aa], i)))
                        except KeyError:
                            modDict[pep_seq] = [tuple((inv_mod_map[aa], i))]
    return modDict

def __promptPTMs(PTMs: list) -> list:
    """prompts the user for which PTMs should be used
    Arguments:
        PTMs {list} -- a list of PTM names
    Returns:
        list -- a list of PTMs specified by the user
    """
    print("\033[96m {}\033[00m".format("Possible PTMs for analysis include:"))
    print("\033[93m {}\033[00m".format("\n ".join([f"[{i}] {ptm}" for i, ptm in enumerate(PTMs)])))
    print("\033[96m {}\033[00m".format("Enter the number for each of the PTMs of interest separated by a comma."))
    ptmIndices = input("\033[96m {}\033[00m".format("For example: 0,1,2\n"))
    if match(r"{^0-9|,}", ptmIndices):  # r"[0-9](?:,[0-9])*"
        print("\033[91m {}\033[00m".format("Invalid input, try again!\n"))
        sleep(.5)
        return __promptPTMs(PTMs)  # NOTE: the end of this recursion is getting it right
    try:
        modPTMs = [PTMs[int(x.strip())] for x in ptmIndices.strip().split(',')]
    except (ValueError, IndexError):
        print("\033[91m {}\033[00m".format("Invalid input, try again!\n"))
        sleep(.5)
        return __promptPTMs(PTMs)
    return modPTMs


def compile_data(search_engine_filepath: str, PTMs: list) -> list:
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
    ###Set input file; perhaps change to an argument and put in function

    #in_file = "MascotExportForRollup-UniHumanRefFASTA-20190731_4Rob.csv"
    in_file = search_engine_filepath #__promptFile()
    enzyme, quant_range, var_mod_map, missed_cleaves = parse_mascot(in_file)

    ###Change this to a variable for later, perhaps let the user choose from var_mod_map
    ###See code below
    #mod_for_quant = "IYn_CleavedPC_TMT6 (C)".lower()
    '''
    if not mod_for_quant: #make this a choice for later
        print(*var_mod_map.items(), sep = '\n')
        mod_choice = input("Choose which modification you want to use for quant (use number): ")
        mod_for_quant = [key for key, value in var_mod_map.items() if value == mod_choice][0]
        print(mod_for_quant)
        input()
    '''
    mods_for_quant = PTMs #__promptPTMs(list(var_mod_map.keys()))

    #assert mod_for_quant in var_mod_map, "Make sure mod is in original file"
    #mod_for_quant_id = var_mod_map.get(mod_for_quant)

    #pep_dict = load_pepdict(enzyme, missed_cleaves)

    #lociDict = dict()

    #print("Quant range", quant_range)#Just see if it found the correct quant columns

    #gene_results = dict()
    #dup_check_set = set()
    #threshold = 0.01 #Make this a mod, too

    input_filename = "quant_parsed.csv"
    input_file =  open(input_filename, 'r')
    file_reader = csv.reader(input_file)
    header = next(file_reader)

    sequence = header.index('pep_seq') #assert this exists later on
    var_mods = header.index('pep_var_mod_pos')
    #score = header.index('pep_expect')
    print(sequence, var_mods)
    #failedSeqs = 0
    #totalSeqs = 0

    #psm_contributions = defaultdict(int) #count number of psms that end up being rolled up

    mod_dict = __parseFile(input_filename, [var_mod_map[mod] for mod in mods_for_quant], var_mod_map)

    return sequence, var_mods, mod_dict, quant_range[0], input_filename, var_mod_map

'''
    for row in file_reader:
        start, stop = quant_range
        dup_check = '.'.join(row[start:stop+1]) #Mascot has duplicate rows; check for these
        totalSeqs += 1
        if dup_check not in dup_check_set and float(row[score]) < threshold:
            dup_check_set.add(dup_check)
            
            pep_seq = row[sequence].replace("L","I")
            genes_positions = pep_dict.get(pep_seq)
            
            if genes_positions and len(genes_positions) == 1: #Gene Specific; maybe don't do this?
                mods = row[var_mods]
                gene, start_pos = list(genes_positions)[0]
                if mods: #Now only deal with mods...
                    #Consider changing if you care about unmodified sequences
                    n, center, c = [x for x in mods.split(".")]
                    for i, aa in enumerate(list(center)): #split mods; note this excludes termini
                        if aa == mod_for_quant_id:
                            gene_true_pos = f"{gene}|{start_pos+i}" #find actual position based on offset
                            qvals = row[quant_range[0]+1:quant_range[1]+2:2]
                            try:
                                #Skip the whole row if any negative values are found
                                #Rationale is that negative values are low and near limit of detection
                                #Can definitely handle this another way if you choose
                                q_test = [float(x) for x in qvals]
                                if any([x < 0 for x in q_test]):
                                    failedSeqs += 1
                                    continue
                            except ValueError:
                                failedSeqs += 1
                                continue

                            #Sum psms
                            tmp_list = gene_results.get(gene_true_pos,
                                                        [0.00 for x in range((quant_range[1]+1-quant_range[0])//2)])
                            qvals = [float(x) if not any([y in x for y in ('-', '#')]) else 0.0 for x in qvals]
                            combined = [x+y for x, y in zip(tmp_list, qvals)]
                            gene_results[gene_true_pos] = combined
                            psm_contributions[gene_true_pos] += 1
                            lociDict[gene_true_pos] = pep_seq
                else:
                    failedSeqs += 1
    print(f"Failed Peptides: {failedSeqs}\nTotal Peptides: {totalSeqs}")

    with open('quantout.csv','w', newline = '') as w1:
        out_writer = csv.writer(w1)

        header = ["Genes","PSMs"]
        start, stop = quant_range
        header.extend(["Sample Quant" for x in range((stop + 1 - start) // 2)])
        header.extend(["sequence"])
        out_writer.writerow(header)

        for i in gene_results:
            writable_row = [i]
            writable_row.extend([psm_contributions[i]])
            writable_row.extend(gene_results[i])
            writable_row.extend([lociDict[i]])
            out_writer.writerow(writable_row)
'''
#EOF