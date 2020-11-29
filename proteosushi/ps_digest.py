"""ps_digest.py: digests in silico the uniprot proteome fasta file"""

__author__ = "Arshag Mooradian"
__email__ = "mooradian@wustl.edu"

import os
from collections import defaultdict
import re
import pickle

from .proteoSushi_constants import cleave_rules

# NOTE: ONLY HANDLES SINGLE ORGANISM. Use ProteoClade for multi-organism/PDX files.

def fasta_producer(proteome_fasta_filepath: str, seq_list: list) -> list:
    """Converts the proteome fasta into a list of tuples with gene, organism,
    sequence, and uniprot ID plainly listed
    
    Arguments:
        proteome_fasta_filepath {str} -- filepath of the fasta proteome file
        seq_list {list} -- List in which (gene, organism, sequence) are stored
        as a tuple.
    Returns:
        list -- the sequence list following modification
    """

    def fasta_entry_check(header, sequence):
        unpid = "N/A"
        organism = None
        protein_name = ""

        # TODO: make it so this is no longer dependent on whether an organism is listed
        # Check organism
        if 'OX=' in header:
            organism = int(header.split('OX=')[1].split()[0])
        elif "OS=" in header:
            organism = header.split('OS=')[1].split(" GN=")[0]

        # Check gene; if no GN, use Uniprot Id
        if 'GN=' in header:
            gene = header.split('GN=')[1].split()[0]
        #elif header.count('|') >= 2:
            # normally would be ==2, but there are some genes with | in them
            #gene = header.split('|')[1]
        else:
            gene = "N/A"

        if header.count('|') >= 2:
            unpid = header.split('|')[1]
        
        if "OS=" in header:
            protein_name = ' '.join(header.split(' ')[1:]).split(" OS=")[0]

        seq_list.append((gene, organism, sequence, unpid, protein_name))

    # open and read in the proteome fasta
    if os.path.exists(proteome_fasta_filepath):
        with open(proteome_fasta_filepath) as input_file:
            print(f"Reading from file: {proteome_fasta_filepath}")
            header = None
            seq = ''
            for line in input_file:
                if line.startswith('>') and seq:
                    fasta_entry_check(header, seq)
                    header = line.strip()
                    seq = ''
                elif line.startswith('>'):  # first entry of file
                    header = line.strip()
                else:
                    seq += line.strip()
            fasta_entry_check(header, seq)  # For last entry
    return seq_list
                


def cleave_rule_determination(rule):
    '''Handle whether the user chooses a built in digest rule or supplies their
    own.

    Arguments:
        rule {str|tuple} --
        if string: a built in rule, ex: "trypsin/p"
        if tuple: a custom rule ("regexcutsites","terminus") ex: ("[RK]","c")

    Returns:
        rules_to_use {tuple} -- Tuple of strings containing
        (regex cutsites, terminus)
    '''
    if rule not in cleave_rules:
        assert type(rule) in (tuple, list), \
            'Cleave rule must be a tuple of strings ("regexrule","terminus"), or be built in.'
        cutsites, direction = rule
        assert direction.lower() in ('c', 'n'), \
            'Second argument of cleave rule needs to be a valid protein terminus.'
        rules_to_use = rule  # assume the user has put in a reg_ex string
    else:
        rules_to_use = cleave_rules.get(rule)
    return rules_to_use


def digest(sequence: str, min_length: int, max_length: int,
           missed_cleavages: int, m_cleave: bool, li_swap: bool,
           rule_to_use: tuple, reverse: bool) -> set:
    '''Take a protein sequence and cut it into pieces, then hash for database
    insertion.

    Arguments:
        sequence {str} -- Protein sequence in all capital letters to chop up.
        min_length {int} -- Minimum amino acid count of peptides to keep for
        database.
        max_length {int} -- Maximum amino acid count of peptides to keep for
        database.
        missed_cleavages {int} -- Number of times a protease is allowed to miss
        a specific site.
        m_cleave {bool} -- Whether or not protein N-terminal methionines are
        removed.
        li_swap {bool} -- Whether peptides will be stored with all leucines
        converted to isoleucines.
        rule_to_use {tuple} -- Tuple of strings with tuple[0] being
        a regex expression for amino acid specificity and tuple[1] as either
        'n' or 'c' cut direction.
        reverse {bool} -- Whether protein sequence is to be reversed before
        storage
    Returns
        cut_set {set} -- Set of integers (hashed peptides) that result from the
        cut rules used.
    '''
    #AMVSEFLKQAWFIENEEQEYVQTVKSSKGGPGSAVSPYPTFNPSSDVAALHKAIMVKGVDEATIIDILTKRNNAQRQQIKAAYLQETGKPLDETLKKALTGHLEEVVLALLKTPAQFDADELRAAMKGLGTDEDTLIEILASRTNKEIRDINRVYREELKRDLAKDITSDTSGDFRNALLSLAKGDRSEDFGVNEDLADSDARALYEAGERRKGTDVNVFNTILTTRSYPQLRRVFQKYTKYSKHDMNKVLDLELKGDIEKCLTAIVKCATSKPAFFAEKLHQAMKGVGTRHKALIRIMVSRSEIDMNDIKAFYQKMYGISLCQAILDETKGDYEKILVALCGGN
    #if "ILVALCGGN" in sequence:
    #    print("ILVALCGGN")
    site_specificity, cut_terminus = rule_to_use
    # Set parameters for a decoy database; only used for ProteoClade really
    if reverse:
        if m_cleave and sequence[-1] == 'M':
            sequence = sequence[:-1]
    else:
        if m_cleave and sequence[0] == 'M':
            sequence = sequence[1:]
    cut_sites = []
    cut_peptides = []
    sites_matched = re.finditer(site_specificity, sequence)
    for site in sites_matched:
        cut_sites.append(site.start())
    # Find sites
    last_site = 0
    if cut_terminus.lower() == 'c':
        for i in cut_sites:
            cut_peptides.append(sequence[last_site:i+1])
            last_site = i + 1
        cut_peptides.append(sequence[last_site:])
    elif cut_terminus.lower() == 'n':
        for i in cut_sites:
            cut_peptides.append(sequence[last_site:i])
            last_site = i
        cut_peptides.append(sequence[last_site:])
    cut_and_missed = list(cut_peptides)  # duplicate to add to for iteration
    #begin_sites = [y - len(x) + 2 for x, y in zip(cut_peptides, cut_sites)]  # Wouldn't it be simpler to add a 1 first and go from there?
    begin_sites = [0] + [y + 2 for x, y in zip(cut_peptides, cut_sites)]
    #if begin_sites:
    #    begin_sites.append(begin_sites[-1] + 1)
    cut_and_missed = [(seq, pos)
                      for seq, pos in zip(cut_and_missed, begin_sites)]
    missed_counter = 1
    while missed_counter <= missed_cleavages:
        missed_peptides = [
            (''.join(cut_peptides[i:i+1+missed_counter]), begin_sites[i])
            for i in range(0, len(cut_peptides) - missed_counter)
            ]  # subtract missed counter here to not duplicate c-terminal peps
        cut_and_missed += missed_peptides
        missed_counter += 1
    if m_cleave is True:
        cut_and_missed = [(x, y + 1) for x, y in cut_and_missed]
    # Filter peptides by length and swap Leucine to Isoleucine, since Mass Spec
    # Can't tell the difference
    if li_swap:
        cut_set = set([
            (i[0].replace("L", "I"), i[1])
            for i in cut_and_missed
            if min_length <= len(i[0]) <= max_length  # Chaining multiple comparators is allowed in python... technically.
            ])
    else:
        cut_set = set([
            (i[0], i[1])
            for i in cut_and_missed
            if min_length <= len(i[0]) <= max_length
            ])
    return cut_set
#EOF