"""proteoSushi_constants.py: provides the constants used in different parts of the program"""

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

annotation_type_dict = {
    "Active_Site_Annotation": 0,
    "Alternative_Sequence_Annotation": 1,
    "Chain_Annotation": 2,
    "Compositional_Bias_Annotation": 3,
    "Disulfide_Bond_Annotation": 4,
    "Domain_Extent_Annotation": 5,
    "Lipidation_Annotation": 6,
    "Metal_Binding_Annotation": 7,
    "Modified_Residue_Annotation": 8,
    "Motif_Annotation": 9,
    "Mutagenesis_Annotation": 10,
    "Natural_Variant_Annotation": 11,
    "NP_Binding_Annotation": 12,
    "Other": 13,
    "Region_Annotation": 14,
    "Repeat_Annotation": 15,
    "Topological_Domain_Annotation": 16,
    "Zinc_Finger_Annotation": 17
}

secondary_annotations = [
    "Beta_Strand_Annotation",
    "Helix_Annotation",
    "Turn_Annotation"
]