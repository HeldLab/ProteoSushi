ProteoSushi GUI should be working again.

# ProteoSushi

ProteoSushi transforms peptide-centric, PTM-enriched peptide data into condensed and annotated PTM site output that is easy to understand and analyze.

ProteoSushi is designed to be simple to use and easy to access. 

## Requirements

Prior to using ProteoSushi, there are some common, minor requirements so it will work correctly. 
The machine running ProteoSushi must have: 

- stable internet access 

- at least 8GB of RAM, 12GB for Windows 10

  - 16GB or higher is recommended, but ultimately depends on the user's memory usage (having other programs open will decrease the memory available for ProteoSushi). More memory allows for more flexibility in the number and scope of programs running at once.
  
  - With only 8GB of RAM, the user will likely need to free up memory by closing other open programs

- Python version 3.8 or higher installed 

  - earlier versions of python may work, but nothing as early as 2.7 will work
  
## Installation

Download python from [here](https://www.python.org/downloads/) and install if needed.

*Those who have python 2.7 installed along with python 3.x will need to replace the python and pip commands with python3 and pip3 respectively*

There are two ways to install ProteoSushi: 

1. Through **PIP**

This is (likely) the easier way to install.

**In Windows:** 

Open the command prompt by typing `cmd` into the search bar (probably at the bottom of the screeen) and clicking on "Command Prompt" when it pops up. Run the command:

`py -m pip install proteosushi`

While in the command prompt, run ProteoSushi with the command:

`py -m proteosushi`

*Tip: If the py command doesn't work, try replacing py with python*

**In MacOS/Linux:**

Open the terminal in MacOS by either searching for it in spotlight or manually finding it in the application list.

In linux, the terminal is among one of the installed applications. If you use linux, you most likely already know how to use the terminal.

Once the terminal is open, run the command:

`python -m pip install proteosushi`

While the terminal is open, run ProteoSushi with the command:

`python -m proteosushi`

2. Through **Github**

**NOTE:** If you have installed ProteoSushi via pip, it will always use that version first. Even if you download the files and run `python run_proteoSushi.py`, it will STILL use your installed version through pip and not use the files you downloaded from GitHub. In this case, you will need to first uninstall the pip version of ProteoSushi with the command `pip uninstall proteosushi` to be able to use the downloaded files.

Download the files directly from [GitHub](https://github.com/HeldLab/ProteoSushi). Click on the Code button and then on Download ZIP.

Once unpacked (unzipped), we will need to use the terminal (command prompt in Windows) to run ProteoSushi.

Before running ProteoSushi, you will need to install some dependencies. In the terminal/command prompt, run the following commands:
```
pip install PyQt5
pip install requests
pip install pandas
```

If that doesn't install correctly, use the following commands:
```
python -m pip install PyQt5
python -m pip install requests
python -m pip install pandas
```

In Windows, run these commands instead:
```
py -m pip install PyQt5
py -m pip install requests
py -m pip install pandas
```

Now we will run ProteoSushi itself.

For **MacOS**, use Finder to navigate to the 'proteosushi' folder in the downloaded files. Right click (or command click) on the 'proteosushi' folder -> Services -> New Terminal at Folder to open a new terminal in the 'proteosushi' folder. Run the command 

`python run_proteoSushi.py`

to start the GUI.

Alternatively, in **MacOS** or **Linux**, you can open the terminal and navigate within the terminal to the 'proteosushi' folder in the downloaded files. Use the cd command to change the current folder as in 

`cd Downloads/ProteoSushi-master/proteosushi` 

and use the ls command to list the contents of the current folder as in

`ls`

, then use the command

`python run_proteoSushi.py`

to run ProteoSushi.

Finally, in **Windows**, first try to open the "run_proteosushi.py" file in the "proteosushi" folder by double-clicking on it. 
In case it doesn't start, open the command prompt by clicking on the search bar in the toolbar at the bottom (usually) of the screen. Type in 'cmd' and click on Command Prompt when it pops up. Once Command Prompt pops up, use the command cd to change to the 'proteosushi' folder in the downloaded files, as in

`cd Downloads/ProteoSushi-master/proteosushi`

and the dir command to list the contents of the current folder, as in

`dir`

, then use the command

`py run_proteoSushi.py`

to run ProteoSushi.

*Again, if the command isn't working, try replacing py with python*

## Files Needed

In order to run ProteoSushi, there are some required files in specific formats. Example files are included in the 'examples' folder in the downloaded files:

- The CSV file output from **Mascot**

  - The file must have the header lines with the information from the search. In order to make sure that the file will be processed correctly, be sure to specify the following settings: 

    - The protease used with the sample 

    - The maximum number of missed cleavages (usually 1 or 2)

    - The variable modifications present in the sample(s)

- Or the txt output folder from **MaxQuant**

  - This folder must have the *summary.txt* and *evidence.txt* files. Other files from the output are not used.

  - *NOTE: It is recommended that you use the newest version of MaxQuant*

- Or the output from any other search engine

  - This file must have a column for peptide sequence named “peptide sequence”

  - Also a column for peptide modified sequence (with PTMs included in parenthesis () or brackets []) named “peptide modified sequence”

  - If you elect to use the quantitation values in the analysis, there must be a column with numbers (presumably positive) with "intensity" or "intensities" in the name

  - Add these columns in if needed

- A FASTA Uniprot Proteome file

  - Ideally this is the same proteome file used in the peptide search

- The following files are optional:

  - A **list of gene names** in a TXT file, one gene name per line

    - This file can be used to prioritize the provided genes whenever there are multiple matches once ProteoSushi performs a search. For example, this can be used to highlight mitochondrial genes if the sample was from the mitochondria

## Using ProteoSushi

Run ProteoSushi and the GUI should pop up and look like this:

![Blank GUI](empty_gui.png)

First, choose the **search engine** used and select the output to use with the window that pops up. For a given search engine, they will need to be one of the following:

- MaxQuant

    - Choose the MaxQuant output folder with the evidence.txt and summary.txt files inside (any extra files will not be used)

    - In the example files from [GitHub](https://github.com/HeldLab/ProteoSushi/tree/master/proteosushi/examples), the maxquant folder is located in `Proteosushi-master/proteosushi/examples` and named "txt".

- Mascot

    - Choose the annotated Mascot output file, should be a CSV file

    - The example file from [GitHub](https://github.com/HeldLab/ProteoSushi/tree/master/proteosushi/examples) is located in `Proteosushi-master/proteosushi/examples` and named "MascotEGFR.csv"

- Generic

    - Choose the output from any search engine, however, it must have a peptide sequence and peptide modified sequence columns (along with a quantitation column labeled “Intensity” if you choose that option)

    - The example file from [GitHub](https://github.com/HeldLab/ProteoSushi/tree/master/proteosushi/examples) is located in `Proteosushi-master/proteosushi/examples` and named "EGFR_Skyline_data.csv"

ProteoSushi will then parse the file to autofill some information (Species ID, Max missed cleavages, protease) and add the option to choose the PTM(s) to use in the analysis

![GUI with Search Engine Selected](search_engine_gui.png)

**Choose the PTM(s)** that will be used in the analysis. The options available are dynamic and will change based on the file provided. (Note: If there are many different PTMs in the file, you may need to move the ProteoSushi window horizontally as they will all be in the same line)

For the example files, the PTMs that appear are different
- For Maxquant, choose the "Carbamidomethyl (C)" PTM 
- For Mascot, choose the "carbamidomethyl (c)" PTM
- For Generic, choose the "C[+57]" PTM

Next, choose the FASTA Uniprot proteome to use in ProteoSushi.

The example files are based on the "Uni-Hum-Ref-20141022.fasta" file located in "Proteosushi-master/proteosushi/examples/fastas", so this is the file that should be used.

If not already filled in, specify the number of **maximum allowed missed cleavages** for a given peptide (usually about 3). 
A higher number will cause the analysis to take longer, but can allow for more matches (possibly). 
It is usually best to stay consistent with what was chosen for the search engine originally.
The Generic search engine option will not fill in the missed cleavages, so for the example, use 3.

If not correctly filled in, specify the **protease** used in the sample digestion step. 
Possible proteases are specified if you hover the cursor over the text here. 
These include:

- trypsin/p

    - Cleaves after Lysine (K) or Arginine (R)

- trypsin!p

    - Cleaves after Lysine (K) or Arginine (R), but not before Proline (P)

- lys-c

    - Cleaves after Lysine (K)

- asp-n

    - Cleaves before Aspartate (D)

- asp-ne

    - Cleaves before Aspartate (D) and Glutamate (E)

- lys-n

    - Cleaves before Lysine (K)

The examples should use trypsin/p

After that, choose the ProteoSushi output file name and location. Click the "Output name and location" button and a window will pop up to choose the location (folder/directory) and type in the name. The output will be saved as a CSV file that can be easily opened in any spreadsheet program.

Following this is the **Options** section with some settings that can be added or ignored based on your analysis. If 

First, if using MaxQuant, a **localization threshold** can be set as a number between 0 and 1. 
In the maxquant *evidence.txt* file, the localization probability is in the columns "*PTM* Probabilities". 
The localization probability indicates the likelihood that the PTM referenced in the header is at the indicated site. 
Any PTM sites with a localization probability below the provided threshold will not be rolled up. 
If this is not provided, ProteoSushi will use the localization determination by MaxQuant.

Second, choose whether to use a prioritized **gene list**. If so, choose the file to be used. 
These genes will be used if there is a tie of annotation score between multiple matches for a PTM site of a peptide. 
If one of the PTM sites is part of a gene from this list, it will be chosen.
This may be helpful if the user is specifically interested in proteins from a specific organelle or pathway.

Third, choose whether to use the **quantitation values**. 
You will need to specify whether to sum or average values that will be combined. 
**The columns used for quantitation must have "intensity" or "intensities" in the header.**

Fourth, choose whether ProteoSushi should **annotated** the rolled-up sites using data from Uniprot

Fifth, if using Mascot or Maxquant, you can optionally specify a **False Discovery Rate (FDR) threshold** between 0 and 1.
In Maxquant, this is taken from the "PEP" column and in Mascot, this is taken from the "pep_expect" column.
More specifically, the PEP column in maxquant means Posterior Error Probability and is the likelihood that the peptide is wrongly assigned when using a target decoy database.
The pep_expect column in Mascot is the [Expectation value of the protein match (PMF only)](https://www.matrixscience.com/help/csv_headers.html).

Once all of the necessary options are included, click on the **Rollup!** button to start the analysis.

After the analysis finishes, you can compare your results to the example files provided in the **[output folder](https://github.com/HeldLab/ProteoSushi/tree/master/proteosushi/output)**.

### Running Examples

A quick overview to run the examples in ProteoSushi.

Example files must be downloaded from [GitHub](https://github.com/HeldLab/ProteoSushi/tree/master/proteosushi/examples).

Install and run ProteoSushi following the instructions listed in **[Installation](https://github.com/HeldLab/ProteoSushi#installation)** and **[Using ProteoSushi](https://github.com/HeldLab/ProteoSushi#using-proteosushi)**. 

Select any of the 3 options for search engine and click on the button that appears on the right.

  - For Maxquant, choose the txt folder listed in the examples folder

  - For Mascot, choose the MascotEGFR.csv file

  - For Generic, choose the EGFR_Skyline_data.csv file

At this point, the "PTMs for Current Analysis" section will pop up and show checkboxes where you can choose any or all of the PTMs as listed. I personally recommend clicking on the carbamidomethyl or c[+57] PTM

Click on the "Uniprot Proteome FASTA" button and choose the "Uni-Hum-Ref-20141022.fasta" file.

The next 3 lines should autofill unless you chose the Generic option earlier. In that case, put 3 in the "Max Missed Cleavages" box.

Choose the output location and filename. Click on the "Output Name and Location" button to pop open the window and type in the filename once it is in the folder you want. Click save.

Anything in the **Options** section can be ignored unless you want to use one of the specific functions. If so, please refer to the appropriate section in **Using ProteoSushi**

Click the "Rollup!" button once you are ready.

Updates as ProteoSushi is running are displayed in the Terminal/Command Prompt.

Once you have the results, you can compare to the example files in the [output folder](https://github.com/HeldLab/ProteoSushi/tree/master/proteosushi/output)

## Results

Results will be returned as a CSV spreadsheet with the name and location based on what the user chose earlier.

![Image of sample results](sample_results.png)

The resulting file will include information for each modified residue of interest including:

- **Gene Name**

- **Site**

  - This is the site of the PTM referenced in this line

- **Protein_Name**

- **Shared_Genes**

  - The other gene names that also match up with this site. They also have a separate entry in the output.

- Target_Genes

  - If the user provided a list of genes to prioritize matching and at least one of those genes matched for this site, it is displayed here. *Note that the annotation score rule is overidden if a target gene is present.*

- **Peptide_Sequence**

  - *Note that each cleaved peptide displayed must have a minimum length of 6 and a maximum length of 55. If the peptide is too short, the shortest cleaved version 6 amino acids or longer will be shown.*

- **Peptide_Modified_Sequence**

  - This is the peptide sequence with the PTM(s) inserted in parentheses `()` or brackets `[]` after the amino acid it is attached to

- **Annotation_Score**

  - The score provided by Uniprot for this Uniprot entry out of 5, with 5 being the highest score (or highest quality annotation), and 1 being the lowest score (or lowest quality annotation)

- **Uniprot_Accession_ID**

  - The unique identifier by Uniprot for a specific protein entry. Information about the protein entry can be accessed by with the URL `www.uniprot.org/uniprot/[Uniprot_Accession_ID]`

- Intensit(y|ies) ... (sum|average)

  - Column(s) showing either the summed or averaged quantification value for the site, if provided and chosen by the user

- Length_Of_Sequence

  - The length of the protein

- Range_of_Interest

  - The number range(s) of relevant peptide sequence(s) for this PTM site

- Region_of_Interest

  - The actual sequence(s) of the same relevant peptide sequence(s) in the last column

- Subcellular_Location

  - The location(s) within the cell where this protein is typically found. For example, mitochondria, nucleus, endoplasmic reticulum, etc.

- Enzyme_Class

  - The type of enzyme that this protein functions as (if applicable)

- rhea

  - Hyperlinks to the [RHEA database](https://www.rhea-db.org)

- Secondary_Structure

  - The protein secondary structure that this PTM site is a part of, if applicable. For example, (alpha) Helix, Beta_Strand, or (beta) Turn.

- Active_Site_Annotation

  - Annotation if the PTM site is an active site
  
- Alternative_Sequence_Annotation

  - Annotation related to protein isoforms

- Chain_Annotation
  
  - Annotation for the chain in the mature protein after processing at the PTM site

- Compositional_Bias_Annotation

  - Annotation indicating overrepresentation of certain amino acids
  
- Disulfide_Bond_Annotation

  - Annotation if the PTM site part of a disulfide bond

- Domain_Extent_Annotation

  - Description of the domain

- Lipidation_Annotation

  - Annotation if the residue is known to be lipidated
  
- Metal_Binding_Annotation

  - Type of metal binding at the PTM site
  
- Modified_Residue_Annotation

  - Known type of PTM at site 

- Motif_Annotation

  - Description of a short, conserved sequence motif of biological significance

- Mutagenesis_Annotation

  - Mutations and known effect on protein function at the site

- Natural_Variant_Annotation

  - Known natural variant at the PTM site (if there is one)

- NP_Binding_Annotation

  - If the site is a known nucleotide binding site

- Other

  - If there is an annotation different than the other _Annotation columns

- Region_Annotation
  
  - Annotation related to the ‘Region_of_Interest’ column

- Repeat_Annotation

  - The types of repeated sequence motifs or repeated domains

- Topological_Domain_Annotation

  - Orientation in the plasma membrane (cytosolic or extracellular)

- Zinc_Finger_Annotation

  - If the modified site is within a zinc finger domain


More detailed information on any of the above annotations is available on the help section of the [UniProt website](https://www.uniprot.org/help).
