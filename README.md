# Extract_NCBI_assembly_metadata

This script will download metadata from a table with assembly accession numbers.

This table should look like (file provided) : 

```
Accessions
GCA_023650665.1
GCA_022132215.1
GCA_007725185.1
```

It can also download the assemblies if you set the option ```-d True```

Prior running the script, you will need ```NCBI-datasets```

In UPPMAX it can be loaded with : ```module load NCBI-datasets/15.29.0```

# Usage example :

```python3 Download_NCBI_assemblies_and_metadata.py -t Accession_table.txt -d False -o /crex/Test_folder/test```

* Will only create a table called ***Genome_metadata.csv*** directed to the directory ***/crex/Test_folder/test***

```python3 Download_NCBI_assemblies_metadata.py -t Accession_table.txt -d True -o /crex/Test_folder/test```

* Will create a table called ***Genome_metadata.csv*** directed to the directory ***/crex/Test_folder/test***
* Will download all the assembly fasta files directed to the directory ***/crex/Test_folder/test***


With the 3 accessions, the Genome_detadata.csv should look like that :
```
        Accessions            Path    geographic location                host isolation source  strain     BioSample Taxonomy                      Organism Supressed    sample type
0  GCA_023650665.1  Not downloaded      Germany:magdeburg        Homo_sapiens       wound swab  319078  SAMN28202791     1648  Erysipelothrix rhusiopathiae        NO        Missing
1  GCA_022132215.1  Not downloaded        Usa:orangebeach  Tursiops_truncatus    dolphin brain  10DISL  SAMN23594042     1648  Erysipelothrix rhusiopathiae        NO  tissue sample
2  GCA_007725185.1  Not downloaded  China:sichuanprovince                 Pig          missing      ZJ  SAMN12347781     1648  Erysipelothrix rhusiopathiae        NO        Missing
```
This code also works for Genbank accessions : 

```Download_Genbank_assemblies_metadata.py -t Accession_table.txt  -o .```
