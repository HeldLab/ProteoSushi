"""combine_output_annotations.py: since annotation is still imperfect, I combine
the patchy output into a file that has as much as possible
"""

import os
import pandas as pd

main_rollup_file = None
annotations_dict = dict()
for f in os.listdir(os.getcwd()):
    if "generic_rollup" in f:
        print(f"Combining {f}")
        rollup_output = pd.read_csv(open(f, 'rU'), engine="python")
        if main_rollup_file is None:
            main_rollup_file = rollup_output
            '''
            for line in main_rollup_file.iterrows():
                try:
                    if line["entry"]:
                        annotations_dict[f"{line["Site"]}|{line["Peptide Modified Sequence"]}"] = line["entry",]
                except KeyError:
                    continue
                '''
        else:
            main_rollup_file.update(rollup_output, overwrite=False)
            '''
            for line in rollup_output.iterrows():
                try:
                    annotations_dict[f"{line["Site"]}|{line["Peptide Modified Sequence"]}"]
                except KeyError:
                    if line["entry"]:
                        annotations_dict[f"{line["Site"]}|{line["Peptide Modified Sequence"]}"] = line["entry",]
                    continue
                '''

main_rollup_file.to_csv("generic_combined_rollup.csv")
#EOF