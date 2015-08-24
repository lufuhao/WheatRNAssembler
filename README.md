SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150603

Requirements:
	Programs: hapcompass.jar hc2vcf.jar, gzip, gunzip, cat, zcat, 
			samtools, freebayes, vcf-tools, express, parallel
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin, Statistics::Basic
	         File::Copy

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--reference|-r	<Fasta>
		[Msg] Sequence file in fasta;
	--cluster|-c	<File>
		[Msg] Cluster file, 1 line 1 cluster, 
		1 cluster may have >=1 sequence ID in reference fasta;
	--list	<ClusterLineNo1-ClusterLineNo2>
		[Opt] Only use cluters from line1 to line2;
	--bam1	<AABBDD.bam.list file>
	--vcf1	<AABBDD.vcf.gz file>
	--fpkm	<FPKM configure file>
	--allele	<Allele configure file>
	--numreads	<file_totalreads>
	--taggenome	<[Opt] Genome tag in AABBDD BAM, default: zg>
	--tagchrom	<[Opt] Chromosome tag in AABBDD BAM, default: ac>
	--vardelimiter	<[Opt] VCF genotype delimiter, default: '/'>

Freebayes
	--freebayespath	<[Opt] /path/to/freebayes if not in PATH>
	--freebayesparallelpath <[Opt] /path/to/freebayes-parallel if not in PATH>
	--freebayes_mincov	<[Opt] minimum coverage, default: 3>
	--freebayes_minalt	<[Opt] minimum alternative count, default: 3>
	--freebayes_parallel	###run freebayes in parallel
	--fastabin	<[Opt] bin size for freebayes-parallel fragment, default: 200>
	
SAMtools
	--tabixpath	<[Opt] /path/to/tabix if not in PATH>
	--bgzippath	<[Opt] /path/to/bgzip if not in PATH>
	--samtoolspath	<[Opt] /path/to/samtools if not in PATH>

VCFtools
	--vcfmergepath	<[Opt] /path/to/vcf-merge if not in PATH>
	--vcfconcatpath	<[Opt] /path/to/vcf-concat if not in PATH>
	--vcfsortpath	<[Opt] /path/to/vcf-sort if not in PATH>
	--vcfqual	<[Opt] minimum VCF quality, default: 20>
	--minmapq	<[Opt] minimum mapping quality, default: 0>

HapCompass
	--javapath	<[Opt] /path/to/java if not in PATH>
	--hapcompasspath	<[Opt] /path/to/hapcompass.jar>
	--hc2vcfpath	<[Opt] /path/to/hc2vcf.jar>

Trinity	
	--trinitypath <[Opt] /path/to/Trinity if not in PATH>
	--maxinsert <[Opt] Int: default 800>

Running LOG
	--logcluster	<[Opt] Cluster running LOG>
	--logfpkm	<[Opt] Reference FPKM LOG>
	--logallele <[Opt] Allele assignment LOG>
	--loggeno	<[Opt] Geno assignment LOG>
	--logcfpkm	<[Opt] Cluster FPKM LOG>

MISC
	--clean	<clean cluster log>
	--phred	<phred64|phred33>
	--threads	<[Opt] Int, default:1>
	--cpus	<[Opt] Int, default:1>
	--debug
	--verbose
		Detailed output for trouble-shooting;
	--version|v!
		Print current SCRIPT version;

Example:
	perl $0 

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
