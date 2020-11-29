"""run_sparql.py: runs the sparql script independent of ProteoSushi"""

from .sparql import sparql_request

from time import sleep

'''with open("sparqllist.tsv", 'r') as sparqllist:
    sparql_input = sparqllist.readlines()
    batch = 1
    i = 0
    while i + batch <= len(sparql_input):
        print(sparql_input[i:i+batch])
        sparql_output = sparql_request(sparql_input[i:i+batch])
        i += batch
        print(f"{float(i)/len(sparql_input)*100}% of results annotated")
        sleep(1)
        print(sparql_output.text)'''

sparql_input = ["(uniprotkb:O00468 317)",
    "(uniprotkb:O00468 948)",
    "(uniprotkb:O00468 844)",
    "(uniprotkb:O00468 103)",
    "(uniprotkb:O00468 296)",
    "(uniprotkb:O00468 1113)",
    "(uniprotkb:O00468 750)",
    "(uniprotkb:O00468 285)",
    "(uniprotkb:O00468 440)",
    "(uniprotkb:O00468 847)",
    "(uniprotkb:O00468 399)",
    "(uniprotkb:P23526 195)",
    "(uniprotkb:O43865 211)",
    "(uniprotkb:O43865 172)",
    "(uniprotkb:O43865 166)",
    "(uniprotkb:O43865 151)",
    "(uniprotkb:O43865 373)",
    "(uniprotkb:Q09666 1833)",
    "(uniprotkb:P02765 219)",
    "(uniprotkb:P02765 132)",
    "(uniprotkb:P02765 230)",
    "(uniprotkb:Q9Y4K1 503)",
    "(uniprotkb:Q12904 161)",
    "(uniprotkb:Q13155 306)",
    "(uniprotkb:O00170 208)",
    "(uniprotkb:O00170 90)",
    "(uniprotkb:P54819 92)",
    "(uniprotkb:P14550 134)",
    "(uniprotkb:P14550 46)",
    "(uniprotkb:P15121 187)",
    "(uniprotkb:P15121 200)",
    "(uniprotkb:P15121 299)",
    "(uniprotkb:P15121 304)",
    "(uniprotkb:Q04828 188)",
    "(uniprotkb:Q04828 206)",
    "(uniprotkb:Q04828 242)",
    "(uniprotkb:Q04828 193)",
    "(uniprotkb:P42330 242)",
    "(uniprotkb:O43488 132)",
    "(uniprotkb:O43488 156)",
    "(uniprotkb:O95154 186)",
    "(uniprotkb:P31749 296)",
    "(uniprotkb:P02768 500)",
    "(uniprotkb:P02768 289)",
    "(uniprotkb:P02768 501)",
    "(uniprotkb:Q13740 313)",
    "(uniprotkb:Q13740 392)",
    "(uniprotkb:Q13740 354)",
    "(uniprotkb:Q13740 270)",
    "(uniprotkb:Q8IZ83 350)"]

sparql_output = sparql_request(sparql_input)
print(sparql_output)
#EOF