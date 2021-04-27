"""run_proteoSushi.py: The starting point to run proteoSushi"""
import logging
import os

try:
    from proteosushi.proteosushi_gui import run_gui
    from proteosushi.combine_intensities import rollup
except ImportError:  # Allows the program to be run not as a module
    from proteosushi_gui import run_gui
    from combine_intensities import rollup


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, filename=os.path.join("logs","ps.log"), filemode='w')
    run_gui()
    '''
    rollup(
        search_engine="generic",
        search_engine_filepath="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/examples/EGFR_Skyline_data.csv",
        use_target_list=False,
        target_list_filepath="",
        max_missed_cleavages=2,
        protease="trypsin/p",
        fdr_threshold=None,
        use_quant=False,
        user_PTMs=["C[+57]", "C[+IAC/NEM]"],
        proteome_fasta_filepath="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/examples/fastas/Uni-Hum-Ref-20141022.fasta",
        intensity_method="sum",
        add_annotation=True,
        species_id="9606",
        output_filename="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/output/skyline_output_annotated.csv",
        localization_threshold=None
    )
    '''
#EOF