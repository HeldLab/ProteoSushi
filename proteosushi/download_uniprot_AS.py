"""download_uniprot_AS.py: downloads a file from uniprot with the annotation score for proteins of 
a species"""

import os
import urllib.parse
import urllib.request

def download_AS_file(species: str) -> str:
    """downloads the uniprotKB file for a species with the annotation score column

    Arguments:
        species {str} -- the species number (e.g. 9606)
    Returns:
        str -- the filepath where the file was downloaded to
    """
    '''
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': 'ACC+ID',
    'format': 'tab',
    "organism": species,
    "columns": "id,genes,annotation_score",
    "compress": "no",
    "query": '',
    "contact": "rseymour@wustl.edu"  # NOTE: not sure if this one works
    }
    '''
    url = "https://www.uniprot.org/uniprot/?query=organism:" + species + "&columns=id,genes,annotation_score&format=tab"
    #url = "https://www.uniprot.org/uniprot/?query=organism:" + species + "+AND+columns=id,genes,annotation_score&format=tab&compress=no&email=rseymour@wustl.edu"
    #headers = {"user-agent": "rseymour@wustl.edu"}
    #data = urllib.parse.urlencode(params)
    #data = data.encode('utf-8')
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        response = f.read()
        annot_score_filename = species + "_annot_score.tsv"
        with open(annot_score_filename, 'w') as as_file:
            decoded_response = response.decode("utf-8")
            if decoded_response[:5] == "Entry":
                as_file.write(decoded_response)
            else:
                return "ERROR: Invalid identifier"
        #print(response.decode('utf-8'))
    return os.path.join(os.getcwd(), annot_score_filename)

if __name__ == "__main__":
    print(download_AS_file("8296"))
#EOF