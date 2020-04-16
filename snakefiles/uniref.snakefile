"""
Download uniref50 and NCBI's taxonomy, and then separate the files
in uniref50 based on their taxonomic rank.

Rob Edwards, April 2020

"""

import os
import sys
import re
from datetime import datetime

parent = {}
taxnames = {1 : "root"}
ranks    = {}

def get_all_nodes(ncbidir, rank, outputfile, logfile):
    """
    Read the NCBI taxonomy nodes file, and find the nodes
    that we want to collect.

    ncbidir should be the path to where the NCBI taxonomy files
    are stored. We will only use "names.dmp" and "nodes.dmp"


    We also save the information about each node and its parents
    :param ncbidir: The directory that has the NCBI taxonomy information.
    :param rank: the rank that we want to collect the information at
    :param output file to write: the ranks output file to write (ranks.tsv)
    :param logfile: where to write the log to
    """

    e = False
    for f in ["nodes.dmp", "names.dmp"]:
        if not os.path.exists(os.path.join(ncbidir, f)):
            sys.stderr.write(f'FATAL: Can not find {os.path.join(ncbidir, f)}\n')
            e = True
    if e:
        sys.stderr.write(f"Did you unpack the NCBI taxonomy into {ncbidir}?\n")
        sys.exit(-1)



    with open(logfile, 'a') as out:
        out.write(f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} Reading taxonomic nodes\n')

    # read the nodes and store the parent trace
    with open(os.path.join(ncbidir, "nodes.dmp"), 'r') as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            if p[2] == rank:
                taxnames[p[2]] = p[2]
            parent[p[0]] = p[1]

    with open(logfile, 'a') as out:
        out.write(f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} Reading taxonomic names\n')

    # read the names
    with open(os.path.join(ncbidir, "names.dmp"), 'r') as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            if p[3] == 'scientific name' and p[0] in taxnames: 
                taxnames[p[0]] = p[1]

    # assign the ranks for everything else
    with open(logfile, 'a') as out:
        out.write(f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} Assigning ranks to everything else\n')

    for t in parent:
        adding = {t}
        tid = t
        while tid not in ranks and parent[tid] not in taxnames:
            adding.add(t)
            tid = parent[tid]
        if tid in ranks:
            for x in adding:
                ranks[x] = ranks[tid]
        else:
            for x in adding:
                ranks[x] = parent[tid]

    with open(logfile, 'a') as out:
        out.write(f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} Writing the output file {outputfile}\n')

    with open(outputfile, 'w') as out:
        for x in ranks:
            out.write(f"{x}\t{ranks[x]}\n")

def stream_fasta(fastafile):
    """
    Stream a fasta file, one read at a time. Saves memory!

    :param fastafile: The fasta file to stream
    :type fastafile: str
    :return:A single read
    :rtype:str, str
    """

    try:
        if fastafile.endswith('.gz'):
            f = gzip.open(fastafile, 'rt')
        else:
            f = open(fastafile, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fastafile)

    posn = 0
    while f:
        # first line should start with >
        idline = f.readline()
        if not idline:
            break
        if not idline.startswith('>'):
            sys.exit("Do not have a fasta file at: {}".format(idline))
        idline = idline.strip().replace('>', '', 1)
        posn = f.tell()
        line = f.readline()
        seq = ""
        while not line.startswith('>'):
            seq += line.strip()
            posn = f.tell()
            line = f.readline()
            if not line:
                break
        f.seek(posn)
        yield idline, seq



def rewrite_output(unireffa, rankf, rank, logfile):
    """
    Rewrite the fasta file into separate files
    :param unireffa: the uniref fasta file
    :param rankf: the ranks file
    :param rank: the desired rank
    :param logfile: where to write the log
    """
    
    logout = open(logfile, 'a')
    logout.write(f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} Starting to rewrite the fasta file\n')

    if not ranks:
        # we got here some other way or ranks.tsv already exists
        with open(rankf, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                ranks[p[0]] = p[1]


    s = re.compile('TaxID=(\d+)')

    logout.write(f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} Rewriting {unireffa}\n')

    fhs = {}
    for seqid, seq in stream_fasta(unireffa):
        m = s.search(seqid)
        rnk = "root"
        if m:
            tid = m.groups()[0]
            if tid in ranks:
                rnk = ranks[tid]
            else:
                logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} ERROR: No rank for {tid}\n")
        else:
            logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} ERROR: No taxonomy in {seqid}\n")

        if rnk not in fhs:
            fhs[rnk] = open(os.path.join(outdir, rnk + ".fasta"), 'w')
        fhs[rnk].write(f">{seqid}\n{seq}\n")

    for fh in fhs:
        fhs[fh].close()

    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Rewriting complete")


if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


NCBIDIR   = config['Paths']['NCBI taxonomy directory']
UNIDIR    = config['Paths']['Uniprot directory']
OUTPUTDIR = config['Paths']['Output directory']
LOGDIR    = config['Paths']['Log directory']
RANK      = config['Taxonomic rank']


rule all:
    input:
        os.path.join(OUTPUTDIR, RANK)


rule write_ranks_file:
    input:
        os.path.join(NCBIDIR, "nodes.dmp"),
        os.path.join(NCBIDIR, "names.dmp")
    output:
        r = os.path.join(OUTPUTDIR, "ranks.tsv")
    log:
        l = os.path.join(LOGDIR, "parsing_taxonomy.log")
    run:
        get_all_nodes(NCBIDIR, RANK, output.r, log.l)

rule get_ncbi_taxonomy:
    output:
        temp(os.path.join(NCBIDIR, "citations.dmp")),
        temp(os.path.join(NCBIDIR, "delnodes.dmp")),
        temp(os.path.join(NCBIDIR, "division.dmp")),
        temp(os.path.join(NCBIDIR, "gencode.dmp")),
        temp(os.path.join(NCBIDIR, "merged.dmp")),
        temp(os.path.join(NCBIDIR, "gc.prt")),
        temp(os.path.join(NCBIDIR, "readme.txt")),
        os.path.join(NCBIDIR, "names.dmp"),
        os.path.join(NCBIDIR, "nodes.dmp")
    params:
        ncbidir = NCBIDIR
    shell:
        """
        cd {params.ncbidir} &&
        curl -LO ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz &&
        curl -LO ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5 &&
        md5sum -c taxdump.tar.gz.md5 &&
        tar xf taxdump.tar.gz
        """

rule get_uniref:
    output:
        os.path.join(UNIDIR, "uniref50.fasta.gz")
    shell:
        "curl -Lo {output} ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"

rule assign_to_taxonomy:
    input:
        u = os.path.join(UNIDIR, "uniref50.fasta.gz"),
        r = os.path.join(OUTPUTDIR, "ranks.tsv")
    output:
        directory(os.path.join(OUTPUTDIR, RANK)),
        touch("Complete")
    log:
        l = os.path.join(LOGDIR, "rewriting_uniref.log")
    run:
        rewrite_output(input.u, input.r, RANK, log.l)


