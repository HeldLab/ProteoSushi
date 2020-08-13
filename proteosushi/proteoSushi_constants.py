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