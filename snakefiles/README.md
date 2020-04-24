# Make a UniRef database

These snakefiles make a UniRef database for you. 

There are several essential steps:

1. Download the [NCBI taxonomy](ftp.ncbi.nlm.nih.gov/pub/taxonomy/) files to the location defined in the [config](uniref_config.yaml) file under Paths : NCBI taxonomy directory
2. Parse the NCBI taxonomy to create a file of the ranks of every tax ID. This assigns the taxids to the rank defined in the [config](uniref_config.yaml) file under `Taxonomic rank`.
3. Download the data from UniProt. We use the [UniRef50](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/) database that is compressed at 50% identity. This is stored in the location defined in the [config](uniref_config.yaml) file under Paths : Uniprot directory.
4. Count the entries in the UniRef50 file
5. Create several temporary files to parallize the parsing of the UniRef. See the notes below.
6. Assign the sequences from each of those temporary files into a file appropriate for the taxonomic rank you have defined.
7. Combine all the output files into individual files, one per taxonomic rank
8. Remove all temporary files.

*Important Notes*:
- The splitting speeds things up, but you can only parallelize things so far since we have to read the disk to get the data. I recommend somewhere between 8-16 cores.
- The splitting is controlled by the number of `cores` as defined by the `--cores` flag.
   - If you are running on a multicore machine, that defines the number of cores you want to use.
   - If you are running on a cluster, that defines the number of nodes you want to use.
   - In either case, I do not recommend setting it much about 16. You will not realize much of a speed up!

# Installation

You will need:
- gunzip
- curl
- awk
- md5sum (yes, we check the file downloads md5sums are correct!)
- tar
- snakemake

All of these are available on most standard Linux installs (except perhaps `snakemake`, which you can install with `pip install --user snakemake`.

# Creating the databases

Run the snakemake command pointing it at the appropriate config file.

On a multicore machine:

```
snakemake --configfile ~/GitHubs/ReducedBacterialProteinsDatabase/snakefiles/uniref_config.yaml --snakefile ~/GitHubs/ReducedBacterialProteinsDatabase/snakefiles/uniref.snakefile --cores 6
```

or on a cluster:

```
snakemake --configfile ~/GitHubs/ReducedBacterialProteinsDatabase/snakefiles/uniref_config.yaml --snakefile ~/GitHubs/ReducedBacterialProteinsDatabase/snakefiles/uniref.snakefile --cluster 'qsub -cwd -o sge_out -e sge_err -V -q important' --cores 12 --local-cores 6 --latency-wait 60
```

(As usual on a cluster, increase `--latency-wait` to ensure NFS completes the file transfers). This will run on 12 nodes of your cluster.
