# BacterialProteins
A database of a reduced set of bacterial (and perhaps archaeal) proteins

This database houses a subset of bacterial proteins, and can be used freely by anyone.

The overall goal is to make a small database that has exemplar bacterial proteins that can be used to screen for some features. It is not expected to be a complete set of bacterial proteins, nor is it expected to contain _every_ protein. However, it should contain _most_ bacterial proteins.

The steps to make the database are:

1. Download the genome summary from PATRIC
curl -Lo genome_summary.tsv ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_summary

2. Identify complete genomes from that list (this eliminates draft genomes with lots of frameshifts and redundancy)
perl -F"\t" -lane 'print $F[0] if ($F[4] eq "Complete")' genome_summary.tsv > complete_ids.txt

3. Create a script called `download_faa.sh`

```
ID=$(head -n $SGE_TASK_ID complete_ids.txt | tail -n 1)
curl -Lo faa/${ID}.PATRIC.faa ftp://ftp.patricbrc.org/genomes/${ID}/${ID}.PATRIC.faa
```

4. Submit this to the cluster and whack PATRIC. You should probably talk to the folks there and let them know they may experience a DOS attack. 

```
mkdir faa
chmod +x ./download_faa.sh
qsub -cwd -o sge_out -e sge_err -t 1-27173:1 ./download_faa.sh
```

5. Convert protein IDs to md5 and write an id map file. By definition this step also dereplicates at 100%. We also ignore proteins with stop codons, and (by default) proteins < 100 amino acids

a. first make `combine_all.sh` that has the python command:
```
python3 ~/EdwardsLab/proteins/protein_md5.py -d faa -i id.map -o proteins.faa
```

b. now submit this to the cluster with a `-hold_jid` for the download jobs:
```
qsub -cwd -o comb.sge.out -e comb.sge.err -hold_jid 20693 ./combine_all.sh
```

6. Combine these with [cd-hit](http://www.cd-hit.org/) at 70% identity

```
cd-hit -d 0 -M 0 -i proteins.faa -o proteins.70.cdhit -c 0.7 -T 0
```

(_Note_: I also tried to do this with `mmseqs2` but it crashed)

7. Extract clusters with a set number of members

e.g. all clusters of 6 or more proteins:
```
python3 extract_clusters.py -f proteins.70.cdhit -c proteins.70.cdhit.clstr -m 5 -v >  proteins.70.cdhit.5members.faa
```

Here are some counts:


File | Number of proteins
--- | ---
After dereplicating at 100% identity | 31,177,581
After dereplicating at 70% identity  | 10,262,128
After removing singletons (clusters with 1 member) | 3,228,900
After removing clusters with 3 or fewer members | 1,352,145
After removing clusters with 5 or fewer members | 846,123
