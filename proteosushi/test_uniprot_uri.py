import requests

req = requests.get("https://www.uniprot.org/uniprot/?query=organism:9606&columns=id,genes,annotation_score&format=tab")