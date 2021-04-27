"""sparql.py: takes the results from proteoSushi and annotates it using the sparql API to access uniprot"""

__author__ = "Rob Seymour, Arshag Mooradian"
__email__ = "rseymour@wustl.edu"

from io import StringIO
import json
import logging

import csv
import pandas as pd
from re import findall
import requests
import sys  # Remove later
from time import sleep
import urllib


# Global constants
ENDPOINT = "https://sparql.uniprot.org"

def sparql_json_to_df(sparql_json_dict: dict):
    """Converts the messy json format into a dataframe

    Arguments:
        sparql_json_dict {dict} -- a json object retrieved through requests.get query
    Returns:
        sparql_df {pd.DataFrame} -- dataframe containing the results
    """
    
    cols = sparql_json_dict["head"]["vars"]

    temp_row_list = [] # list holding rows
    for row in sparql_json_dict["results"]["bindings"]:
        item = []
        for c in cols:
            item.append(row.get(c, {}).get("value"))
        temp_row_list.append(item)

    sparql_df = pd.DataFrame(temp_row_list, columns = cols)
    return sparql_df


def sparql_request(unpid_site_list: list):
    """Sends the request to uniprot to get data on each of the sites
    Arguments:
        unpidSiteList {list} -- a list of tuples with uniprot ID and site of PTM
    Returns:
        pandas.DataFrame -- the complete annotated dataframe to be applied to the rolled-up PTM sites
    """
    unpid_site_list_str = ""
    if not unpid_site_list:
        return None
    for tup in unpid_site_list:
        if len(tup) >= 2 and tup[1] == '0':
            tup[1] = '1'
        #if tup[0] == "P60709" or tup[0] == "P63261":
        #    print("Multiple")
        if len(tup[0]) >= 5:
            unpid_site_list_str += f"(uniprotkb:{tup[0]} {tup[1]})\n"
    # TODO: Delete this later
    #with open("sparqllist.tsv", 'w') as sparqllist:
    #    sparqllist.write(unpid_site_list_str)
    #print(unpid_site_list_str)



    prefix="""PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX ec: <http://purl.uniprot.org/enzyme/>"""

    #TODO: Add up:Similarity_Annotation, 
    new_entry_length_pos_query = prefix + """
SELECT
    ?entry
    ?position ### New the position of the C that was given in the values
    (STRLEN(?iupac) AS ?lengthOfSequence)
    ?begin
    ?end
    (SUBSTR(?iupac, ?begin, ?end - ?begin + 1) AS ?regionOfInterest)
FROM <http://sparql.uniprot.org/uniprot>
WHERE {
    VALUES (?entry ?position) {""" + unpid_site_list_str + """}
    ?entry up:sequence ?sequence .
    ?annotation up:range ?range .
    ?range faldo:begin
        [ faldo:position ?begin ; faldo:reference ?sequence ] ;
            faldo:end
        [ faldo:position ?end ; faldo:reference ?sequence ] .
    FILTER (?begin <= ?position && ?position <= ?end)
    # get the IUPAC AAs associated with the identifier
    ?sequence rdf:value ?iupac .
    # get the AA subsequence of the annotation as a new variable
} ORDER BY DESC(?entry) ASC(?position)"""

    query_stripped_entry_pos = prefix + """
SELECT
    ?entry
    ?position ### New the position of the C that was given in the values
    ?lengthOfSequence
    ?begin
    ?end
    ?regionOfInterest ### subsequence of the annotation that had the C
FROM <http://sparql.uniprot.org/uniprot>
WHERE {
     VALUES (?entry ?position) {""" + unpid_site_list_str + """}
     ?entry up:sequence ?sequence .
     ?annotation up:range ?range1, ?range2;
                a ?type .
     ?range1 faldo:begin
         [ faldo:position ?begin ; faldo:reference ?sequence ] .
     ?range2 faldo:end
         [ faldo:position ?end ; faldo:reference ?sequence ] .
     FILTER (?begin <= ?position && ?position <= ?end)
     # get the IUPAC AAs associated with the identifier
     ?sequence rdf:value ?iupac .
} ORDER BY DESC(?entry) ASC(?position)"""



    query_entry_subcellular = prefix + """
SELECT
    ?entry 
    ?location
FROM <http://sparql.uniprot.org/uniprot>
WHERE {
    VALUES (?entry ?position) {""" + unpid_site_list_str + """
    }
    ?entry up:annotation ?subAnnotation .
    ?subAnnotation a up:Subcellular_Location_Annotation .
    ?subAnnotation up:locatedIn/up:cellularComponent ?location .
} ORDER BY DESC(?entry)"""




    query_ec_rhea_type = prefix + """
SELECT
    ?entry
    ?position
    ?ec
    ?rhea
FROM <http://sparql.uniprot.org/uniprot>
WHERE {
    VALUES (?entry ?position) {""" + unpid_site_list_str + """
    }
    ?entry up:sequence ?sequence .
    ?entry up:annotation ?caAnnotation .
    ?caAnnotation up:catalyticActivity ?ca .

    OPTIONAL {
        ?ca up:enzymeClass ?ec .
    }

    OPTIONAL {
        ?ca up:catalyzedReaction ?rhea .
    }

#Order desc by entry and ascending by position afterwards
} ORDER BY DESC(?entry) ASC(?position)"""  # NOTE: remove the ordering to help reduce time for query

    query_comment = prefix + """
SELECT
    ?entry
    ?position
    ?type
    ?comment
FROM <http://sparql.uniprot.org/uniprot>
WHERE {
    VALUES (?entry ?position) {""" + unpid_site_list_str + """}
    ?entry up:sequence ?sequence .
    ?entry up:annotation ?annotation .
    ?annotation up:range ?range ;
               a ?type .
    ?range faldo:begin
        [ faldo:position ?begin ; faldo:reference ?sequence ] ;
            faldo:end
        [ faldo:position ?end ; faldo:reference ?sequence ] .
    FILTER (?begin <= ?position && ?position <= ?end)
	OPTIONAL {
        ?annotation rdfs:comment ?comment .
    }
#Order desc by entry and ascending by position afterwards
} ORDER BY DESC(?entry) ASC(?position)"""

    # Grabs the annotation by parts to maximize the amount we receive from the uniprot server
    region_annot = request_annot(new_entry_length_pos_query)
    region_annot_stripped = request_annot(query_stripped_entry_pos)
    if len(region_annot.index) < len(region_annot_stripped.index):  # This should check to see if there are rows missing in the full query
        region_annot = region_annot.append(region_annot_stripped)
        region_annot.drop_duplicates(subset=["entry", " position", " begin", " end"], keep="first", inplace=True)
    subcell_annot = request_annot(query_entry_subcellular)
    extras_annot = request_annot(query_ec_rhea_type)
    comment_annot = request_annot(query_comment)
    # Creates blank dataframes if uniprot did not return that info
    if region_annot is None or list(region_annot.columns) != ["entry", " position", " lengthOfSequence", " begin", " end", " regionOfInterest"]:
        region_annot = pd.DataFrame(columns=["entry", " position", " lengthOfSequence", " begin", " end", " regionOfInterest"])
    else:
        region_annot.dropna(how="all", inplace=True)
    if subcell_annot is None or list(subcell_annot.columns) != ["entry", " location"]:
        subcell_annot = pd.DataFrame(columns=["entry", " location"])
    else:
        subcell_annot.dropna(how="any", inplace=True)
    if extras_annot is None or list(extras_annot.columns) != ["entry", " position", " ec", " rhea"]:
       extras_annot = pd.DataFrame(columns=["entry", " position", " ec", " rhea"])
    else:
        extras_annot.dropna(how="all", inplace=True)  # TODO: I will likely need to change this back to all later
    
    if comment_annot is None or list(comment_annot.columns) != ["entry", " position", " type", " comment"]:
        comment_annot = pd.DataFrame(columns=["entry", " position", " type", " comment"])
    else:
        comment_annot.dropna(how="all", inplace=True)

    # Full outer joins the annotations to preserve all info possible
    full_annot = region_annot
    del(region_annot)
    try:
        full_annot = full_annot.merge(subcell_annot, how="outer", on="entry")
        del(subcell_annot)
        full_annot = full_annot.merge(extras_annot, how="outer", on=["entry", " position"])
        del(extras_annot)
        full_annot = full_annot.merge(comment_annot, how="outer", on=["entry", " position"])
        del(comment_annot)
    except KeyError:  # The proper columns are missing in one of the query results; usually means uniprot isn't working correctly
        #print(full_annot)
        if isinstance(full_annot, pd.DataFrame) and "502 Proxy Error" in full_annot.iloc[1][0]:
            logging.exception("502 Error")
            print("502 ERROR: try Uniprot annotations later")
            return 502
        logging.exception("Uniprot data formatting error: See sparql_request function in sparql.py")
        print("Uniprot data formatting error: please notify ProteoSushi administrator")
        print(list(full_annot.columns))
        sys.exit()
    #print(full_annot.columns)
    return full_annot

def request_annot(query: str, attempts_left=10):
    """Sends the query to uniprot and checks the output

    Arguments:
        query {str} -- the query with the sites to annotate
        attempts_left {int} -- the number of requests that can be sent out
    Returns:
        pandas.dataframe -- the returned annotation
    """
    headers = {"user-agent": "rseymour@wustl.edu"}
    time_to_sleep = 5
    try:
        r = requests.post(ENDPOINT, data = {"format": "csv", "query": query}, headers = headers)
        if "<!DOCTYPE html SYSTEM \"about:legacy-compat\">" in r.text:  # This just means it returned a 500 error
            if attempts_left > 0:
                sleep(time_to_sleep)
                print("\033[93m {}\033[00m".format(f"\nAttempts used: {11 - attempts_left}"), end='')
                return request_annot(query, attempts_left - 1)
            else:
                return None  # TODO: remove this once Uniprot fixes their stuff
                raise pd.errors.ParserError("Ran out of Uniprot accession attempts")
        csv_file = StringIO(r.text)
        sparql_from_csv_df = pd.read_csv(csv_file)
    except (pd.errors.ParserError, requests.exceptions.ConnectionError):
        if attempts_left > 0:
            sleep(time_to_sleep)
            print("\033[93m {}\033[00m".format(f"\nAttempts used: {11 - attempts_left}"), end='')
            return request_annot(query, attempts_left - 1)
        else:
            return None  # TODO: remove this once Uniprot fixes their stuff
            raise pd.errors.ParserError("Ran out of Uniprot accession attempts")
    return sparql_from_csv_df

def process_sparql_output(output_df, sparql_dict: dict) -> list:
    """Process the output from uniprot to make it consistent with the rollup output

    Arguments:
        output_df {pandas.dataframe} -- the string output from uniprot, may have multiple
        sparql_dict {dict} -- the growing dictionary of uniprot annotations
    Returns:
        list -- a list of annotations with the data properly separated
    """
    #print(f"\nIn process_sparql_output beginning uniprot annotation is {type(output_df)}")
    logging.debug(f"\nIn process_sparql_output beginning uniprot annotation is {type(output_df)}")

    def __parse_sparql_line(line: str) -> str:
        """changes the complicated uniprot annotation output into a cleaner line

        Arguments:
            line {str} -- a line from the output
        Returns:
            str -- a cleaned up version of that same line
        """
        #entry,lengthOfSequence,catalyicActivity,location,ec,rhea,type,comment,position,begin,end,regionOfInterest
        new_str = list()
        for chunk in line[1]:
            if isinstance(chunk, str) and "http://purl.uniprot.org/uniprot/" in chunk:
                new_str.append(chunk.replace("http://purl.uniprot.org/uniprot/", ""))
            elif isinstance(chunk, str) and "http://purl.uniprot.org/core/" in chunk:
                new_str.append(chunk.replace("http://purl.uniprot.org/core/", ""))
            elif isinstance(chunk, str) and "http://purl.uniprot.org/locations/" in chunk:
                new_chunk = chunk.replace("http://purl.uniprot.org/locations/", "")
                try:
                    new_str.append(subcellular_location_dict[new_chunk.zfill(4)])  # Converts the number location to descriptor
                except KeyError:
                    new_str.append("")
            elif isinstance(chunk, str) and "http://purl.uniprot.org/enzyme/" in chunk:
                new_chunk = chunk.replace("http://purl.uniprot.org/enzyme/", "")
                new_str.append(ec_id_to_description_dict[new_chunk])
            elif isinstance(chunk, str) and "^^<http://www.w3.org/2001/XMLSchema#int>" in chunk:
                new_str.append(findall(r"(\d+?)\^", chunk)[0])
            else:
                new_str.append(str(chunk))
        return '\t'.join(new_str)

    output_list = list()
    comments_dict = dict()
    #sparql_dict = dict()  # This gets used by the main program to connect these annotations to the rest of the data.
    if len(output_df) == 1 or (not isinstance(output_df, str) and output_df.empty):
        print("\nUniprot Error: Failed to get annotations", end='')
        return output_list, sparql_dict
    
    try:
        position_index = output_df.columns.get_loc(" position")
    except (AttributeError, KeyError):
        logging.exception("Uniprot annotations missing position column or wrong format")
        print("\nUniprot Error: Missing the Position column", end='')
        logging.debug(f"In process_sparql_output error uniprot annotation is {type(output_df)}")
        #print(f"\nIn process_sparql_output error uniprot annotation is {type(output_df)}")
        return output_list, sparql_dict

    # get the subcellular location from UniProt
    def retrieve_uniprot_subcellular_location(attempts_left=5):
        '''Retrieve subcellular location information from UniProt directly'''
        if attempts_left <= 0:
            return 4
        try:
            r = requests.get('https://www.uniprot.org/locations/?format=tab')
        except requests.exceptions.ConnectionError:
            return retrieve_uniprot_subcellular_location(attempts_left-1)
        subcellular_location_df = pd.read_csv(StringIO(r.text), sep = '\t')  # TODO: potential ParserError
        try:
            subcellular_location_df["location_id_reformatted"] = (
                                                                subcellular_location_df["Subcellular location ID"]
                                                                .str.split("-")
                                                                .str[1]
                                                                )
        except KeyError:
            subcellular_location_df["location_id_reformatted"] = ("")
        return subcellular_location_df
    
    subcellular_location_df = retrieve_uniprot_subcellular_location()
    if isinstance(subcellular_location_df, int) and subcellular_location_df == 4:
        return 4, None  # This is an error related to uniprots's 
    subcellular_location_dict = dict(subcellular_location_df[["location_id_reformatted", "Alias"]].values.tolist())

    # Fetches a file from expasy to deal with the 'ec' column
    try:
        r = urllib.request.urlopen("ftp://ftp.expasy.org/databases/enzyme/enzyme.dat", "enzyme.dat")
    except urllib.error.URLError:
        sleep(10)
        r = urllib.request.urlopen("ftp://ftp.expasy.org/databases/enzyme/enzyme.dat", "enzyme.dat")
    enzyme_dat_list = r.read().decode("utf-8").split('\n')

    # Turn enzyme.dat into a dictionary for mapping
    ec_id_to_description_dict = dict()
    for line_number, line in enumerate(enzyme_dat_list):
        if line.startswith("ID "):
            current_id = line[5:]
            current_description = enzyme_dat_list[line_number + 1][5:]
            ec_id_to_description_dict[current_id] = current_description

    #entry,position,lengthOfSequence,begin,end,regionOfInterest,location,ec,rhea,type,comment
    #header = output_lines.colnames()
    #print(output_lines.columns)
    
    #i = 1
    repeat_index = 3
    #while i < len(output_lines):
    for output_line in output_df.iterrows():  # TODO: optimize by vectorizing
        output_line = __parse_sparql_line(output_line)
        # This should only have 1 iteration and is used to properly split the tsv type line
        for output_split in csv.reader([output_line], delimiter='\t', quotechar='"'):
            key = ','.join([*output_split[:repeat_index]])
            try:
                comments_dict[key] += output_split[repeat_index:]
            except KeyError:
                comments_dict[key] = output_split[repeat_index:]
        #i += 1
    
    for key in comments_dict:
        output_list.append(key + ',' + ','.join(comments_dict[key]))
        sparql_dict[ key.split(',')[0] + '|' + key.split(',')[1] ] = key.split(',') + comments_dict[key] #unpid, pos
    #print(f"sparql_dict has {sparql_dict}")
    return output_list, sparql_dict
#EOF