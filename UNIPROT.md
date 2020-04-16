# Using the UniRef dataset

A different way to make a reduced dataset is to use someone that others have already made!

[UniProt](https://www.uniprot.org/) has the terrific [UniRef](https://www.uniprot.org/uniref/) clusters that are automatically updated.

They release several versions of those data (these descriptions are from the above website)

- [uniref 100](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/)
   - UniRef100 combines identical sequences and sub-fragments with 11 or more residues from any organism into a single UniRef entry.
- [uniref 90](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/)
   - UniRef90 is built by clustering UniRef100 sequences such that each cluster is composed of sequences that have at least 90% sequence identity to, and 80% overlap with, the longest sequence (a.k.a. seed sequence).
- [uniref 50](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/)
   - UniRef50 is built by clustering UniRef90 seed sequences that have at least 50% sequence identity to, and 80% overlap with, the longest sequence in the cluster.

The data is in a very standard format, with a header line, like this example

```
>UniRef50_A0A345MT53 Adenylation kDNA ligase-like protein n=1 Tax=CrAssphage sp. TaxID=2212563 RepID=A0A345MT53_9CAUD
```

the fields `n` (number of proteins in the group), `Tax` (taxonomy string), `TaxID` (taxonomy ID), and `RepID` are present in every entry.

Here is the `Virus` entry with the most number of proteins:

```
>UniRef50_P0C0U1 Protein PB1-F2 n=32816 Tax=Influenza A virus TaxID=11320 RepID=PB1F2_I34A1
```

With that information, it is relatively trivial to split these into groups based on their taxonomic profiles. Basically, we need to get the taxonomy from NCBI, and then
trace up through that taxonomy to our desired level, then we can split the fasta file into smaller files, one per taxonomy.

For example, if we split at the `superkingdom` level, we end up with:

- root.fasta
- Archaea.fasta
- Bacteria.fasta
- Eukaryota.fasta
- Viruses.fasta

In this case, `root.fasta` contains those sequences that are conserved across at least two domains, or can not be assigned a taxonomy for some reason. Perhaps not unsurprisingly, a lot of Virus sequences end up in root, because they are present both in viruses and in another domain (e.g. phage genes are common in Bacteria, of course).

We have made a snakefile to prepare these files for you, at any rank you would like. (For a complete list of ranks, see the `nodes.dmp` file in the NCBI taxonomy). However, be aware that you will end up with more and more sequences in the `root.fasta` as you go further down the taxonomy.

We typically start by splitting at the `superkingdom` level (to get the files noted above), and then resplit those at different phylogenetic levels. 

# Installation

To run the snakefile you will need 
- [snakemake](https://snakemake.readthedocs.io/)
- [cURL](https://curl.haxx.se/) [this is probably already installed on your computer]
- [md5sum](https://en.wikipedia.org/wiki/Md5sum) [this is probably already installed on your computer]
- [tar](https://www.gnu.org/software/tar/) [this is probably already installed on your computer]
- The [snakefile](uniref.snakefile)
- The [config file](uniref_config.yaml)

# Creating the databases


You need to edit the config file to set the appropriate paths. This config file is in [YAML](https://yaml.org/), and only has a few parameters:

```
---
Paths:
	'NCBI taxonomy directory': 'ncbi'
	'Uniprot directory' : 'uniref'
	'Output directory': 'Results'
	'Log directory': 'log'

'Taxonomic rank': 'superkingdom'
```

You can set the directories to whatever you would like, that is where the files are written to. 

The `Taxonomic rank` option is which of the ranks you would like to include. The defaul, `superkingdom` is a good place to start so that you don't end up with too many proteins in `root`.

To run the command, just run it like any other `snakefile`:

```
snakemake --configfile uniref_config.yaml --snakefile uniref.snakefile --cores 4
```

You don't need that many cores as there are not many options to parallelize things. More than one helps as the `uniref50` can be downloading in parallel with the NCBI taxonomy being downloaded and parsed.

We run this on our cluster with the command like this:

```
snakemake --configfile uniref_config.yaml --snakefile uniref.snakefile \ 
--cluster 'qsub -cwd -o sge_out -e sge_err -V -q important' --cores 10 --local-cores 4 --latency-wait 60 
```
















