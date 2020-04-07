# Some Bacterial Proteins

[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-Green)](https://edwards.sdsu.edu/research)
[![DOI](https://www.zenodo.org/badge/252589739.svg)](https://www.zenodo.org/badge/latestdoi/252589739)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


We have created a database of a reduced set of bacterial (and perhaps archaeal) proteins that can be used freely by anyone.

Our overall goal is to make a small database that has exemplar bacterial proteins that can be used to screen for some features. It is not expected to be a complete set of bacterial proteins, nor is it expected to contain _every_ protein. However, it should contain _most_ bacterial proteins, depending on how you define that. Essentially we want to screen contigs against common bacterial proteins, and are less interested in the single hypothetical proteins that contribute a lot of proteins.

The database uses the upper case protein sequence's [md5sum](https://en.wikipedia.org/wiki/Md5sum) as an identifier for the sequence. There are many advantages to this: in particular it is computable based on the protein sequence, so you can easily test if you have seen the sequence before. There is a minor issue using the md5sum, which is that there is a chance of false positives - you may think you have seen a sequence before but you have not. However, in testing that did not occur with this protein database.

We make a separate file with the IDs and their definition lines from PATRIC, called [id.map](id.map.gz) that you can download separately. You should probably store that in a [SQLite](https://www.sqlite.org/) or similar database for rapid access, but we haven't provided that yet. Let us know if you need it. 


# Release 2020/03/28 Statistics

You can download the individual fasta files:
- [proteins.70.cdhit.gz](https://edwards.sdsu.edu/data/BacterialProteins/latest/id.map.gz) [1.9 G] Fasta file of the proteins dereplicated at 70%.
- [id.map](https://edwards.sdsu.edu/data/BacterialProteins/latest/id.map.gz) [2.7G] This file contains the id and the sequence definition line
- [proteins.70.cdhit.clstr.gz](https://edwards.sdsu.edu/data/BacterialProteins/latest/proteins.70.cdhit.clstr.gz) [760M] The cd-hit cluster definitions including all proteins and the clusters to which they belong.

## Reduced data sets

You can reduce this set further by extracting clusters with at least _n_ members.

e.g. all clusters of 6 or more proteins (note the parameter `-m` is the minimum that must be exceeded):

```
python3 scripts/extract_clusters.py -f proteins.70.cdhit -c proteins.70.cdhit.clstr -m 5 -v >  proteins.70.cdhit.5members.faa
```

We have not provided these files, but you can easily recreate them from the files available above.

File | Number of proteins
--- | ---
All proteins in complete genomes | 
After dereplicating at 100% identity | 31,177,581
After dereplicating at 70% identity  | 10,262,128
After removing singletons (clusters with 1 member) | 3,228,900
After removing clusters with 3 or fewer members | 1,352,145
After removing clusters with 5 or fewer members | 846,123




# Creating the database

The steps to make the database are:

1. Download the genome summary from PATRIC
```
curl -Lo genome_summary.tsv ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_summary
```

2. Identify complete genomes from that list (this eliminates draft genomes with lots of frameshifts and redundancy)
```
perl -F"\t" -lane 'print $F[0] if ($F[4] eq "Complete")' genome_summary.tsv > complete_ids.txt
```

3. Create a script called `download_faa.sh`

```
ID=$(head -n $SGE_TASK_ID complete_ids.txt | tail -n 1)
curl -Lo faa/${ID}.PATRIC.faa ftp://ftp.patricbrc.org/genomes/${ID}/${ID}.PATRIC.faa
```

4. Submit this to a cluster and whack PATRIC to download fasta protein files of all the complete bacteria. You should probably talk to the folks there and let them know they may experience a DOS attack. 

```
mkdir faa
chmod +x ./download_faa.sh
qsub -cwd -o sge_out -e sge_err -t 1-27173:1 ./download_faa.sh
```

5. Convert protein IDs to md5 and write an id map file. By definition this step also dereplicates at 100%. We also ignore proteins with stop codons, and (by default) proteins < 100 amino acids

a. first make `combine_all.sh` that has the python command:
```
python3 scripts/protein_md5.py -d faa -i id.map -o proteins.faa
```

b. now submit this to the cluster with a `-hold_jid` for the download jobs (since we are using SGE):

```
qsub -cwd -o comb.sge.out -e comb.sge.err -hold_jid 20693 ./combine_all.sh
```

6. Combine these with [cd-hit](http://www.cd-hit.org/) at 70% identity

```
cd-hit -d 0 -M 0 -i proteins.faa -o proteins.70.cdhit -c 0.7 -T 0
```

(_Note_: I also tried to do this with `mmseqs2` but it crashed)

