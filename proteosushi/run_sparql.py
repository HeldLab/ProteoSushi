"""run_sparql.py: runs the sparql script independent of ProteoSushi"""

from sparql_script import sparql_request

from time import sleep

with open("sparqllist.tsv", 'r') as sparqllist:
    sparql_input = sparqllist.readlines()
    batch = 1
    i = 0
    while i + batch <= len(sparql_input):
        print(sparql_input[i:i+batch])
        sparql_output = sparql_request(sparql_input[i:i+batch])
        i += batch
        print(f"{float(i)/len(sparql_input)*100}% of results annotated")
        sleep(1)
        print(sparql_output.text)
#EOF