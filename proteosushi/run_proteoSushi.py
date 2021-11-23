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
    try:
        os.mkdir("logs")
    except FileExistsError:
        pass
    logging.basicConfig(level=logging.DEBUG, filename=os.path.join("logs","ps.log"), filemode='w')
    #run_gui(multithread=True)
    #'''
    rollup(
        search_engine="generic",
        #search_engine_filepath="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/examples/EGFR_Skyline_data.csv",
        search_engine_filepath="/home/rob/Downloads/Proteosushi_test_deamid_seqOnly.csv",
        use_target_list=False,
        target_list_filepath="",
        max_missed_cleavages=2,
        protease="trypsin/p",
        fdr_threshold=None,
        use_quant=False,
        #user_PTMs=["C[+57]", "C[+IAC/NEM]"],
        user_PTMs=["N[+0.984]"],
        #proteome_fasta_filepath="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/examples/fastas/Uni-Hum-Ref-20141022.fasta",
        proteome_fasta_filepath="/home/rob/Downloads/9606_HomoSapiens_150721.fasta",
        #intensity_method="sum",
        intensity_method="",
        add_annotation=True,
        species_id="9606",
        output_filename="/home/rob/Documents/Held_Lab/ProteoSushi/proteosushi/output/ps_test_211120.csv",
        localization_threshold=None
    )
    #'''
#EOF