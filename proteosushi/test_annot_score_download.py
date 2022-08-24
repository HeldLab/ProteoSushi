"""test_annot_score_download.py: runs the download_uniprot_AS.py script independent of ProteoSushi"""

try:
    from .download_uniprot_AS import download_AS_file
except ImportError:
    from download_uniprot_AS import download_AS_file

print(download_AS_file(species="8296"))
#EOF