"""download_uniprot_AS.py: downloads a file from uniprot with the annotation score for proteins of 
a species"""

import os
import urllib.parse
import urllib.request

def download_AS_file(species: str, attempts_left = 5) -> str:
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
    if attempts_left <= 0:
        print("\033[91m {}\033[00m".format("Unable to retrieve Uniprot Annotation Score! Proceeding regardless..."))
        return ""
    url = "https://www.uniprot.org/uniprot/?query=organism:" + species + "&columns=id,genes,annotation_score&format=tab"
    #url = "https://www.uniprot.org/uniprot/?query=organism:" + species + "+AND+columns=id,genes,annotation_score&format=tab&compress=no&email=rseymour@wustl.edu"
    #headers = {"user-agent": "rseymour@wustl.edu"}
    #data = urllib.parse.urlencode(params)
    #data = data.encode('utf-8')
    try:
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
    except urllib.error.URLError:
        return download_AS_file(species, attempts_left-1)
    return os.path.join(os.getcwd(), annot_score_filename)

if __name__ == "__main__":
    print(download_AS_file("8296"))
#EOF