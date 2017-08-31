# PSCodeHPEV
The model implements the mechanism of co-transcriptional folding in the context of an assemblysome (Shakeel et al, 2017) in which nascent RNAs fold and contact coat-protein as the nucleotides are incorporated by the polymerase. The presented code is an example that was used with the HPeV1-8 circulating strain variants with the publication Leonov et al, 2017.

To execute the code used for the analysis of the circulating HPeV1-8 strains using PS-factor 4 and PS-factor 1 (with and without affinity for protein during cotranscriptional folding), check the docstring at the top in the Python module for dependencies, and execute as:

python determine_PSs_cotranslationally.py HPEV.fasta
