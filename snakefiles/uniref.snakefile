"""
Download uniref50 and NCBI's taxonomy, and then separate the files
in uniref50 based on their taxonomic rank.

Rob Edwards, April 2020

"""

import os
import sys
import re
import gzip
from datetime import datetime

parent = {}
taxnames = {'1' : "root", '131567' : "CellularOrganisms"}
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

    
    logout = open(logfile, 'w')

    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Reading taxonomic nodes\n")

    # read the nodes and store the parent trace
    with open(os.path.join(ncbidir, "nodes.dmp"), 'r') as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            if p[2] == rank:
                taxnames[p[0]] = p[2]
            parent[p[0]] = p[1]
    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Reading taxonomic names\n")
    logout.write(f"Taxnames are {taxnames}\n")

    # read the names
    with open(os.path.join(ncbidir, "names.dmp"), 'r') as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            if p[3] == 'scientific name' and p[0] in taxnames and p[0] != '131567': 
                taxnames[p[0]] = p[1]

    # assign the ranks for everything else
    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Assigning ranks to everything else\n")
    logout.write(f"Taxnames are {taxnames}\n")

    for t in parent:
        adding = {t}
        tid = t
        while parent[tid] != '1' and parent[tid] != '131567' and tid not in ranks:
            adding.add(tid)
            tid = parent[tid]

        if tid == '2':
            logout.write(f"Found {tid} with parent {parent[tid]}\n")
            if tid in taxnames:
                logout.write(f"\t{tid} is in taxnames as {taxnames[tid]}\n")
        
        if tid in ranks:
            logout.write(f"Already had {ranks[tid]} for {tid}\n")
            for x in adding:
                ranks[x] = ranks[tid]
        elif tid in taxnames:
            for x in adding:
                ranks[x] = taxnames[tid]
        else:
            logout.write(f"No taxonomy for {tid} (type: {type(tid)}. Returned root\n")
            for x in adding:
                ranks[x] = "root"


    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Writing the output file {outputfile}\n")

    with open(outputfile, 'w') as out:
        for x in ranks:
            out.write(f"{x}\t{ranks[x]}\n")
    
    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Parsing taxonomy complete {outputfile}\n")
    logout.close()

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


def rewrite_output(unireffa, rankf, rank, outdir, logfile):
    """
    Rewrite the fasta file into separate files
    :param unireffa: the uniref fasta file
    :param rankf: the ranks file
    :param rank: the desired rank
    :param outdir: the output directory
    :param logfile: where to write the log
    """
    
    logout = open(logfile, 'a')
    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Starting to rewrite the fasta file\n")

    if not ranks:
        # we got here some other way or ranks.tsv already exists
        with open(rankf, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                ranks[p[0]] = p[1]


    sys.stderr.write(f"Outdir is a {type(outdir)} with values {outdir}\n")
    sys.stderr.write(f"Uniref is a {type(unireffa)} with values {unireffa}\n")

    s = re.compile('TaxID=(\d+)')
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    logout.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')} Rewriting {unireffa}\n")

    fhs = {"unknown" : open(os.path.join(outdir, "Unknown.fasta"), 'w')}
    for seqid, seq in stream_fasta(unireffa):
        # deal with "unknown" sequences
        if 'Tax=unknown TaxID= ' in seqid:
            fhs["unknown"].write(f">{seqid}\n{seq}\n")
            continue

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
REMOVETMP = config['Remove Temporary Files']
NUMFILES  = config['Number of Temporary Files'] + 1

# NUMFILES = min(workflow.cores, workflow.nodes) + 1

rule all:
    input:
        expand(os.path.join(OUTPUTDIR, RANK, "{filenumber}"), filenumber=range(1, NUMFILES)),
        "Complete"

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
    shell:
        """
        cd {NCBIDIR} &&
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


rule count_uniref:
    """
    We make this its own rule so we only run it once!
    """
    input:
        u = os.path.join(UNIDIR, "uniref50.fasta.gz")
    output:
        # ue = os.path.join(UNIDIR, "tmp", "uniref50.entries")
        ue = temp(os.path.join(UNIDIR, "tmp", "uniref50.entries")) if REMOVETMP else os.path.join(UNIDIR, "tmp", "uniref50.entries")
    shell:
        "gunzip -c {input.u} | grep -c \> > {output.ue}"

rule split_uniref_to_files:
    """ 
    By using {filenumber} as a wildcard, we will start this FILENUMBER
    times, but it allows us to remove the outputs when we are done
    """
    input:
        u = os.path.join(UNIDIR, "uniref50.fasta.gz"),
        r = os.path.join(UNIDIR, "tmp", "uniref50.entries")
    output:
        # fn = os.path.join(UNIDIR, "tmp", "uniref50.split.{filenumber}.fa")
        fn = temp(os.path.join(UNIDIR, "tmp", "uniref50.split.{filenumber}.fa")) if REMOVETMP else os.path.join(UNIDIR, "tmp", "uniref50.split.{filenumber}.fa")
    params:
        n = NUMFILES-1,
        d = os.path.join(UNIDIR, "tmp"),
    shell:
        """
        RECORDS=$(cat {input.r});
        RPF=$((RECORDS/{params.n}));
        RPF=$((RPF+10));
        gunzip -c {input.u} | \
        awk -v rpf=$RPF -v fname={output.fn} 'BEGIN {{n=0;c=0}} \
                /^>/ {{ \
                    if(n%rpf==0){{ \
                        c++; \
                        file=sprintf("{params.d}/uniref50.split.%d.fa",c); \
                    }}; \
                    n++; \
                }} \
                {{ if (file==fname) {{ print >> file; }}; fflush(); }};'
        """

rule assign_to_taxonomy:
    input:
        u = os.path.join(UNIDIR, "tmp", "uniref50.split.{filenumber}.fa"),
        r = os.path.join(OUTPUTDIR, "ranks.tsv")
    output:
        # outdir = directory(os.path.join(OUTPUTDIR, RANK, "{filenumber}"))
        outdir = temp(directory(os.path.join(OUTPUTDIR, RANK, "{filenumber}"))) if REMOVETMP else directory(os.path.join(OUTPUTDIR, RANK, "{filenumber}"))
    log:
        l = os.path.join(LOGDIR, "rewriting_uniref.{filenumber}.log")
    run:
        rewrite_output(input.u, input.r, RANK, output.outdir, log.l)

rule concat_taxonomy:
    input:
        expand(os.path.join(OUTPUTDIR, RANK, "{count}"), count=range(1, NUMFILES))
    output:
        #o = touch("concat_taxonomy_complete")
        o = temp(touch("concat_taxonomy_complete")) if REMOVETMP else touch("concat_taxonomy_complete")
    params:
        odir = os.path.join(OUTPUTDIR, RANK)
    log:
        l = os.path.join(LOGDIR, "concatenating.log")
    shell:
        """
        for F in $(find {params.odir} -type f -printf "%f\n" | sort -u); do 
            echo `date`: {params.odir}/*/$F >> {log.l};
            cat {params.odir}/*/$F > {params.odir}/$F;
        done
        """

rule dont_remove_tempfiles_too_soon:
    """
    This is somewhat of a pseudo-rule to prevent removing the 
    temporary results files too soon. 
    We need to wait until concat_taxonomy is complete
    for _all_ data before we remove _any_ of the directories
    """
    input:
        "concat_taxonomy_complete",
        expand(os.path.join(OUTPUTDIR, RANK, "{count}"), count=range(1, NUMFILES))
    output:
        touch("Complete")

