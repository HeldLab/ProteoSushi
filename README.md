# ProteoSushi

ProteoSushi is designed to be simple to use and easy to access. Prior to using ProteoSushi, there are some common, minor requirements so it will work correctly. The machine running ProteoSushi must have stable internet access and at least 4GB of RAM (note here explaining that 8GB or higher is recommended, but truly depends on how the user uses the machine. More is better up to a certain point). In addition, the computer being used must have Python version 3.8 or higher installed (note here explaining that earlier versions of python may work, but nothing as early as 2.7 will work). Install from https://www.python.org/downloads/.

In order to run ProteoSushi, there are some required files in specific formats:

- The CSV file output from Mascot

  - The file must have the header lines with the information from the search, such as the protease used with the sample and the maximum number of missed cleavages

- Or the output folder from MaxQuant

  - This folder must have the summary.txt and evidence.txt files. Other files from the output are not used.

- Or the output from any other search engine

  - This file must have a column for peptide sequence named “peptide sequence”

  - Also a column for peptide modified sequence (with PTMs included) named “peptide modified sequence”

  - If you elect to use the quantitation values in the analysis, there must be a column for this as well

  - Add these columns in if needed

- A FASTA Uniprot proteome file

  - Ideally this is the same proteome file used in the peptide search

- The following files are optional:

  - A list of gene names in a TXT file, one gene name per line

    - This file can be used to prioritize the provided genes whenever there are multiple matches once ProteoSushi performs a search

Once you have the required files, there are two ways to install ProteoSushi. The (likely) easier way is by opening the terminal (or command prompt in Windows) and running the command

`pip install proteosushi`

While in the terminal, run python with the command

`python`

(or start python in Windows) and run the next two commands to import the module and run the GUI:

`from run_proteosushi import run_proteosushi`

`run_proteosushi()`

Alternatively, download the files directly from [GitHub](https://github.com/heldlab/proteosushi). 
Once unpacked, start the terminal in the proteosushi folder and run the command
`python run_proteosushi.py`
to start the GUI.

At this point, the GUI should pop up and look like this:
![](empty_gui.png)

First, choose the search engine used and select the output to use with the window that pops up. For a given search engine, they will need to be one of the following:

- Mascot

    - Choose the annotated Mascot output file, should be a CSV file

- MaxQuant

    - Choose the MaxQuant output folder with the evidence.txt and summary.txt files inside

- Generic

    - Choose the output from any search engine, however, it must have a peptide sequence and peptide modified sequence columns (along with a quantitation column labeled “Intensity” if you choose that option)

ProteoSushi will then parse the file to autofill some information (Max missed cleavages, protease) and add the option to choose the PTM(s) to use in the analysis

![](search_engine_gui.png)

Choose the PTM(s) that will be used in the analysis. The options available are dynamic and will change based on the file provided. ( note: If there are many different PTMs in the file, you may need to move the ProteoSushi window horizontally as they will all be in the same line)

Next, choose the FASTA Uniprot proteome to use in ProteoSushi.

Following this is the **Options** section with some settings that can be added or ignored based on your analysis.

First, choose whether to use a prioritized gene list. If so, choose the file to be used. These genes will be used if there is a tie of annotation score between multiple matches for a PTM site of a peptide. If one of the PTM sites is part of a gene from this list, it will be chosen.

Second, choose whether to use the quantitation values. You will need to specify whether to sum or average values that will be combined. The columns used for quantitation must have "intensity" in the header.

If not already filled in, specify the number of maximum allowed missed cleavages for a given peptide (usually about 3). A higher number will cause the analysis to take longer, but can allow for more matches (possibly). It is usually best to stay consistent with what was chosen for the search engine originally.

If not already filled in, specify the protease used in the sample digestion step. Possible proteases are specified if you hover the cursor over the text here. These include:

-trypsin/p

-trypsin!p

-lys-c

-asp-n

-asp-nc

-lys-n

Specify the threshold for FDR, if using Mascot or Maxquant. This value can be left blank if you do not want to specify a threshold

Once all of the necessary options are included, click on the “Rollup!” button to start the analysis

Results will be returned as a CSV spreadsheet with one of the following filenames depending on the search engine used:

- “MQ_Rollup.csv” for MaxQuant

- “Mascot_Rollup.csv” for Mascot

- “Other_Rollup.csv” for any other search engine

The resulting file will include information for each modified residue of interest including:

![Image of sample results](sample_results.png)

- Peptide

- Uniprot ID

- Gene Name

- Secondary structure at residue

- Whether the residue is involved in a binding site

- Cellular location of the protein
