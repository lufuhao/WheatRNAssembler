SYNOPSIS:

perl $0 [Options]
Version: LUFUHAO20160204

Requirements:
	Linux: java, gzip, gunzip, cat, zcat, parallel, 
	Programs: samtools, bgzip, tabix, CDBtools, bam2fastq
		Trinity, mira, cap3, cdhit
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin, Statistics::Basic
	         File::Copy, Bio::DB::Sam, File::Path, File::Copy, Storable, Data::Dumper
	         FuhaoPerl5Lib
	Script: bam_multiextract.sh

Descriptions:
	Assembly homoeolog-aware ESTs for polyploid transcriptome. Before running this script, 
	You need to do following steps:
	1. 1st homoeolog-UNaware assembly
	2. Mapping reads back
	3. Call variation
	4. Determine anceltral alleles for each ploidy, and/or subgenome info for each read
	5. Calculate EST abundance [ using express, RSEM, etc ] for polyplody and the ancestrals
	6. Cluster ESTs [ wcd, peace, etc ]
	
	This script will do:
	1. Extract BAM and 1st assembly for each EST cluster
	2. Determine ploidy based on the read subgenome info/ancestrol EST abundance
	3. Find which allele belong to which subgenome by read subgenome info/ancestral alleles
	4. Classify all the reads into each ploidy and shared group
	5. Assemble each ploidy with shared reads
	6. Final clean

Options:
	--help|-h
	    Print this help/usage;
	--reference|-r  <Fasta>
	    [Msg] Sequence file in fasta;
	--cluster|-c    <File>
	    [Msg] Cluster file, 1 line 1 cluster, first column is cluster number
	    1 cluster may have >=1 sequence ID in reference fasta;
	--list          <ClusterLineNo1-ClusterLineNo2>
	    [Opt] Only use cluters from ClusterNo1 to ClusterNo2;
	    OR provide a file containing all cluster numbers, 1 number perl line
	--bam           <bam/bamlist file>
	    [Msg] BAM of reads to reference.fa
	--vcf           <variation.vcf.gz file>
	    [Msg] VCF file called from freebayes, INDEL '-' is not supported
	--fpkm          <FPKM configure file>
	    [Msg] EST abundance for each EST cluster
	--allele        <Allele configure file>
	    [Msg] Anceltral alleles for each ploidy
	--numreads      <file_totalreads>
	    [Msg] Total reads number for each ReadGroup referred in BAM file
	--taggenome     <Genome tag>
	    [Msg] Genome tag in annotated BAM, default: zg
	--tagchrom      <Chromosome tag>
	    [Msg] Chromosome tag in BAM, default: zc

SAMtools
	--tabixpath     <[Opt] /path/to/tabix if not in PATH, default: tabix>
	--bgzippath     <[Opt] /path/to/bgzip if not in PATH, default: bgzip>
	--samtoolspath  <[Opt] /path/to/samtools if not in PATH, default: samtools>

VCF
	--vcfqual       <INT>
	    [Opt] minimum VCF quality, default: 20
	--minmapq       <INT>
	    [Opt] minimum mapping quality, default: 0

HapCompass
	--javapath      <[Opt] /path/to/java if not in PATH, default: java>
	--javamem       <[Opt] memory for hapcompass, default: 2G>
	--hapcompasspath  <[Msg] /path/to/hapcompass.jar>

bam2fastq
	--bam2fastqpath <[Opt] /path/to/bam2fastq if not in PATH, default: bam2fastq>

Trinity
	--trinitypath   <[Opt] /path/to/Trinity if not in PATH, default: Trinity>
	--maxinsert     <[Opt] INT, default 800>
	--cpus          <[Opt] INT, default:1>

MIRA4
	--mira4path     <[Opt] /path/to/mira if not in PATH, default: mira>
	--threads       <[Opt] INT, default:1>

CAP3
	--cap3path      <[Opt] /path/to/cap3 if not in PATH, default: cap3>

LOG_Output
	--logcluster    <FILE_OUT>
	    [Opt] Cluster running LOG output, default: 0.cluster.log
	--logfpkm       <FILE_OUT>
	    [Opt] Reference FPKM LOG output, default: 0.reference.fpkm.log
	--logallele     <FILE_OUT>
	    [Opt] Allele assignment LOG output, default: 0.allele.log
	--loggeno       <FILE_OUT>
	    [Opt] Geno assignment LOG output, default: 0.geno.log
	--logcfpkm      <FILE_OUT>
	    [Opt] Cluster FPKM LOG output, default: 0.cluster.FPKM.log

MISC
	--minexp        <INT>
	    <[Opt] Minimum EST abundance value
	--clean
	    [Opt] Clean temporary folder and files. 
	    But maybe problematic when running on Linux Cluster or mounted NFS
	--debug
	--verbose
	    [Opt] Detailed output for trouble-shooting;

Example:
	perl $0 \\
	    -r ref.fa -c est.cluster --bam reads.1stasm.bam --vcf bam.vcf \\
	    --fpkm bam.express --allele ancestral.alleles --list 1-100 \\
	    --hapcompasspath "/path/to/hapcompass.jar" \\
	    --numreads "/path/to/bam.totalreads" \\
	    --trinitypath trinity --threads 5 --cpus 5

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

Copyright
	Copyright (c) 2015-2018 Fu-Hao Lu
########################################################

File format

--fpkm
#reference\tABD\tAB\A\D
contig1	10	6	5	4
contig2	15	0	0	6
...



--allele
#reference\tposition\tAallele\tBallele\tDallele
#Allele: only number accespted referring to VCF gile
contig1	98	0	1	0
contig1	206	1	1	0
......



--totalreads
Total	10000000
KroDevGrainPool1	500
