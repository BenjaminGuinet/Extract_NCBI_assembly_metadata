# Extract_NCBI_assembly_metadata

This script will download metadata from a table with assembly accession numbers.

It can also download the assemblies.

# Usage example :

’’’python3 Download_NCBI_assemblies.py -t Accession_table.txt -d False -o /crex/Test_folder/test’’’

* Will only create a table called ***Genome_metadata.csv*** directed to the directory ***/crex/Test_folder/test***

’’’python3 Download_NCBI_assemblies.py -t Accession_table.txt -d True -o /crex/Test_folder/test’’’

* Will create a table called ***Genome_metadata.csv*** directed to the directory ***/crex/Test_folder/test***
* Will download all the assembly fasta files directed to the directory ***/crex/Test_folder/test***
