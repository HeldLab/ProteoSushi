"""run_proteoSushi.py: The starting point to run proteoSushi"""
try:
    from proteosushi.proteosushi_gui import run_gui
    from proteosushi.combine_intensities import rollup
except ImportError:  # Allows the program to be run not as a module
    from proteosushi_gui import run_gui
    from combine_intensities import rollup


if __name__ == "__main__":
    run_gui()
    '''
    rollup(
        search_engine="maxquant",
        search_engine_filepath="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/examples/txt",
        use_target_list=False,
        target_list_filepath="",
        max_missed_cleavages=3,
        protease="trypsin/p",
        fdr_threshold=.01,
        use_quant=True,
        user_PTMs=["Carbamidomethyl (C)"],
        proteome_fasta_filepath="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/examples/fastas/Uni-Hum-Ref-20141022.fasta",
        intensity_method="sum",
        add_annotation=False,
        species_id="9606",
        output_filename="/home/rob/Documents/Held_Lab/ProteoSushi/localization_test.csv",
        localization_threshold=.5
    )
    '''
#EOF