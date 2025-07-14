Peptide.py is a Python script that filters out proteins with medium or low FDR confidence, as well as peptides with an odd number of dehydro of Cysteine, due to their low probability. It only keeps peptides that have "Unambiguous" for "PSM Ambiguity" and "High" for "Confidence (by Search Engine): Sequest HT". Peptide.py then merges the peptide isoforms into peptides without considering the modification of peptide isoforms and sums up their abundances. Finally, Pepitde.py produces the most abundant peptide among all exclusive peptides of each protein, which is considered the candidate mature peptide for that protein. The exclusive peptide must not be labeled as "No Quan Values", "Not Reliable", or "Shared" in the "Quan Info" column.

Usage: Make sure the input file "Search_file.xlsx" is in this directory, and execute the following command in the Linux terminal: "./Peptide.py". The output folder is "Peptide"

The output folder "Peptide" consists of the following files/folders:
1. File "Master_proteins.txt": the original list of master proteins
2. Folder "Peptides": each master protein and its corresponding list of filtered peptides
3. Folder "Peptides_sorting": merges the isoforms of filtered peptides into peptides without considering the modification of peptide isoforms, and sums up their abundances, and then sorts them by abundance from high to low.
4. Folder "Peptides_fasta": outputs the amino acid sequence of each file in "Peptides_sorting"
5. Folder "Mature_Peptides": the list of candidate mature peptides

