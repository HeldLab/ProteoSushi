"""run_proteoSushi.py: The starting point to run proteoSushi"""

from proteosushi_gui import run_gui
'''
def run_proteoSushi(search_engine: str, search_engine_filepath: str, use_target_list: bool, 
                    target_list_filepath: str, max_missed_cleavages: int, protease: str, fdr_threshold: float):
    """starts proteoSushi rollup when called by the GUI

    Arguments:
        search_engine {str} -- "maxquant", "mascot", or "generic"
        search_engine_filepath {str} -- the filepath for search engine output
        use_target_list {bool} -- whether to use the target list to prioritize matches
        target_list_filepath {str} -- the filepath for the target gene list
        max_missed_cleavages {int} -- the maximum allowed missed cleavages in a peptide
        protease {str} -- the protease used previously to cleave the proteins in the sample
        fdr_threshold {float} -- threshold used for the pep_expect and PEP columns in mascot/maxquant
    """
    pass
'''
if __name__ == "__main__":
    run_gui()
#EOF