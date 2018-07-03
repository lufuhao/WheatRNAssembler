#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use Bio::DB::Sam;
use Scalar::Util;
use File::Path qw/remove_tree/;
use File::Copy qw/move/;
use FuhaoPerl5Lib::BamKit;
use FuhaoPerl5Lib::VcfKit;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::FileKit;
use FuhaoPerl5Lib::FastaKit;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150603

Requirements:
	Programs: hapcompass.jar hc2vcf.jar, gzip, gunzip, cat, zcat, 
			samtools, freebayes, vcf-tools, express, parallel
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin, Statistics::Basic

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
	--bam2	<AABB.bam.list file>
	--bam3	<AA.bam.list file>
	--bam4	<DD.bam.list file>
	--vcf1	<AABBDD.vcf.gz file>
	--vcf2	<AABB.vcf.gz file>
	--vcf3	<AA.vcf.gz file>
	--vcf4	<DD.vcf.gz file>
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

HapTree
	--haptreepath	<[Opt] /path/to/haptree if not in PATH>
	
MIRA4
	--mirapath	<[Opt] /path/to/mira if not in PATH>
	
CD-HIT-EST
	--cdhitestpath	<[Opt] /path/to/cd-hit-est if not in PATH>

Bowtie2
	--bowtie2buildpath	<[Opt] /path/to/bowtie2-build if not in PATH>
	--bowtie2path	<[Opt] /path/to/bowtie2 if not in PATH>

Running LOG
	--logcluster	<[Opt] Cluster running LOG>
	--logfpkm	<[Opt] express LOG>
	--logallele <[Opt] Allele assignment LOG>

MISC
	--phred	<phred64|phred33>
	--threads	<[Opt] Int, default:1>
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

EOH
#Express
#	--express_path	<[Opt] /PATH/to/express>
#	--express_insert	<[Opt] mean insert size, default: 250>
#	--express_stdev	<[Opt] insert size stdev, default: 300>
#	--express_readlen	<[Opt] read_length, default: 150>
#HapCompass
#	--javapath	<[Opt] /path/to/java if not in PATH>
#	--hapcompasspath	<[Opt] /path/to/hapcompass.jar>
#	--hc2vcfpath	<[Opt] /path/to/hc2vcf.jar>

###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $numthreads, $debug, $ver);
#Global
my ($reference, $file_cluster, $list);
my ($file_bam_aabbdd, $file_bam_aabb, $file_bam_aa, $file_bam_dd);###filelist
my ($file_vcf_aabbdd, $file_vcf_aabb, $file_vcf_aa, $file_vcf_dd);
my ($geno_delimiter);###vcf format
my ($bam_genome_tag, $bam_chromo_tag);###Bam
#CDBtools
my ($path_cdbfasta, $path_cdbyank);
#express
my ($path_express, $express_frag_len_mean, $express_frag_len_stddev, $express_max_read_len, $minimum_fpkm);###express
#freebayes
my ($path_freebayes, $freebayes_min_coverage, $freebayes_min_alternative_count, $vcfqual, $min_mapq,$freebayes_parallel, $path_freebayesparallel, $freebayes_fastabin);###freebayes
#samtools
my ($path_samtools, $path_tabix, $path_bgzip);###samtools
#HapCompass
my ($path_hapcompassjar, $path_hc2vcf, $path_java);###HapCompass
#vcf-tools
my ($path_vcfmerge, $path_vcfconcat, $path_vcfsort);#VCFtools
#haptree
my ($path_haptree);
#MIRA4
my ($path_mira4);###MIRA4
#CDHIT
my ($path_cdhitest);
#Bowtie2
my ($pathbowtie2build, $pathbowtie2p, $maxinsert, $fq_qual_score);
#LOG
my ($cluster_log, $fpkm_log, $allele_log);###LOG
my ($totalreads_AABBDD, $totalreads_AABB, $totalreads_AA, $totalreads_DD);

GetOptions(
	"help|h!" => \$help,
	"reference|r:s" => \$reference,
	"cluster|c:s" => \$file_cluster,
	"list:s" => \$list,
	"bam1:s" => \$file_bam_aabbdd,
	"bam2:s" => \$file_bam_aabb, 
	"bam3:s" => \$file_bam_aa, 
	"bam4:s" => \$file_bam_dd,
	"vcf1:s" => \$file_vcf_aabbdd, 
	"vcf2:s" => \$file_vcf_aabb, 
	"vcf3:s" => \$file_vcf_aa, 
	"vcf4:s" => \$file_vcf_dd,
	"taggenome:s" => \$bam_genome_tag,
	"tagchrom:s" => \$bam_chromo_tag,
	"vardelimiter:s" => \$geno_delimiter,
	"cdbfastapath:s" => \$path_cdbfasta,
	"cdbyankpath:s" => \$path_cdbyank,
#	"express_path:s" => \$path_express,
#	"express_insert:i" => \$express_frag_len_mean,
#	"express_stdev:i" => \$express_frag_len_stddev,
#	"express_readlen:i" => \$express_max_read_len,
	"freebayespath:s" => \$path_freebayes,
	"freebayesparallelpath:s" => \$path_freebayesparallel,
	"tabixpath:s" => \$path_tabix,
	"bgzippath:s" => \$path_bgzip,
	"freebayes_mincov:i" => \$freebayes_min_coverage,
	"freebayes_minalt:i" => \$freebayes_min_alternative_count,
	"freebayes_parallel!" => \$freebayes_parallel,
	"fastabin:i" =>	\$freebayes_fastabin,
	"vcfqual:i" => \$vcfqual,
	"minmapq:i" => \$min_mapq,
	"samtoolspath:s" => \$path_samtools,
#	"javapath:s" => \$path_java,
#	"hapcompasspath:s" => \$path_hapcompassjar,
#	"hc2vcfpath:s" => \$path_hc2vcf,
	"haptreepath:s" => \$path_haptree,
	"vcfmergepath:s" => \$path_vcfmerge,
	"vcfconcatpath:s" => \$path_vcfconcat,
	"vcfsortpath:s" => \$path_vcfsort,
	"mirapath:s" => \$path_mira4,
	"cdhitestpath:s" => \$path_cdhitest,
	"bowtie2buildpath:s" => \$pathbowtie2build,
	"bowtie2path:s" => \$pathbowtie2p,
	"logcluster:s" => \$cluster_log,
	"logfpkm:s" => \$fpkm_log,
	"logallele:s" => \$allele_log,
	"phred:s" => \$fq_qual_score,
	"maxinsert:i" => \$maxinsert,
	"threads|t:i" => \$numthreads,
	"debug|" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
my $RootDir=$Bin;
my $cluster_start=0;
my $cluster_end;
($cluster_start, $cluster_end)=split(/-/, $list);
$debug=0 unless (defined $debug);
$numthreads=1 unless (defined $numthreads);
$geno_delimiter=',' unless (defined $geno_delimiter);
$bam_genome_tag='zg' unless (defined $bam_genome_tag);
$bam_chromo_tag='zc'unless (defined $bam_chromo_tag);
$fq_qual_score='phred33' unless (defined $fq_qual_score);
my $cmd_success_code=1;
my $cmd_failure_code=0;
### CDBtools
$path_cdbfasta='cdbfasta' unless (defined $path_cdbfasta);
$path_cdbyank='cdbyank' unless (defined $path_cdbyank);
###Express
$path_express='express' unless (defined $path_express);
$express_frag_len_mean=250 unless (defined $express_frag_len_mean);
$express_frag_len_stddev=300 unless (defined $express_frag_len_stddev);
$express_max_read_len=150 unless (defined $express_max_read_len);
$minimum_fpkm=0 unless (defined $minimum_fpkm);
###Freebayes
$path_freebayes='freebayes' unless (defined $path_freebayes);
$freebayes_min_coverage=3 unless (defined $freebayes_min_coverage);
$freebayes_min_alternative_count=3 unless (defined $freebayes_min_alternative_count);
$vcfqual=20 unless (defined $vcfqual);
$min_mapq=0 unless (defined $min_mapq);
$freebayes_parallel=0 unless (defined $freebayes_parallel);
$path_freebayesparallel='freebayes-parallel' unless (defined $path_freebayesparallel);
$freebayes_fastabin=200 unless (defined $freebayes_fastabin);
#samtools
$path_samtools='samtools' unless (defined $path_samtools);
$path_tabix='tabix' unless (defined $path_tabix);
$path_bgzip='bgzip' unless (defined $path_bgzip);
#HapCompass
$path_hapcompassjar="$RootDir/utils/hapcompass/hapcompass.jar" unless (defined $path_hapcompassjar);
$path_hc2vcf="$RootDir/utils/hapcompass/hc2vcf.jar" unless (defined $path_hc2vcf);
$path_java='java' unless (defined $path_java);
#VCFtools
$path_vcfmerge='vcf-merge' unless (defined $path_vcfmerge);
$path_vcfconcat='vcf-concat' unless (defined $path_vcfconcat);
$path_vcfsort='vcf-sort'  unless (defined $path_vcfsort);
#HapTree
$path_haptree='haptree' unless (defined $path_haptree);
#MIRA4
$path_mira4='mira' unless (defined $path_mira4);
#CDHITEST
$path_cdhitest='cd-hit-est' unless (defined $path_cdhitest);
#Bowtie2
$pathbowtie2build='bowtie2-build' unless (defined $pathbowtie2build);
$pathbowtie2p='bowtie2' unless (defined $pathbowtie2p);
$maxinsert=800 unless (defined $maxinsert);
#LOG
$cluster_log='0.cluster.log' unless (defined $cluster_log);
#print "### ClusterLog: $cluster_log\n";### For test ###
$fpkm_log='0.fpkm.log' unless (defined $fpkm_log);
$allele_log='0.allele.log' unless (defined $allele_log);
### Number of total reads
$totalreads_AABBDD=939286414 unless (defined $totalreads_AABBDD);
$totalreads_AABB=618630534 unless (defined $totalreads_AABB);
$totalreads_AA=622332862 unless (defined $totalreads_AA);
$totalreads_DD=1059289085 unless (defined $totalreads_DD);
my @totalreads=($totalreads_AABBDD, $totalreads_AABB, $totalreads_AA, $totalreads_DD);



### input and output ################################################
#Check important file inpput

die "Error: undefined reference file\n" unless (defined $reference);
($reference)=AddFilePath($reference);
die "Error: not-found reference file\n" unless (-s $reference);
die "Please specify cluster file\n" unless (defined $file_cluster and -s $file_cluster);
die "Please specify the file of bam list for AABBDD\n" unless (defined $file_bam_aabbdd and -s $file_bam_aabbdd);
die "Please specify the file of bam list for AABB\n" unless (defined $file_bam_aabb and -s $file_bam_aabb);
die "Please specify the file of bam list for AA\n" unless (defined $file_bam_aa and -s $file_bam_aa);
die "Please specify the file of bam list for DD\n" unless (defined $file_bam_dd and -s $file_bam_dd);

#Remove last-run cluster log 
unlink $cluster_log if (-e $cluster_log);
unlink $fpkm_log if (-e $fpkm_log);

#read AABBDD bam files
my @bam_AABBDD=();
if ($file_bam_aabbdd=~/\.bam$/i) {
	push (@bam_AABBDD, $file_bam_aabbdd);
}
else {
	close BAMAABBDD if (defined fileno(BAMAABBDD));
	open (BAMAABBDD, "<$file_bam_aabbdd") || die "InputOutputError: can not open AABBDD $file_bam_aabbdd\n";
	while (my $line1=<BAMAABBDD>) {
		chomp $line1;
		if (defined $line1 and $line1 ne '' and -s $line1) {
			push (@bam_AABBDD, $line1);
		}
		else {
			die "InputOutputError: bam $line1 in AABBDD $file_bam_aabbdd not exist or empty\n";
		}
	}
	close BAMAABBDD;
}
die "InputOutputError: empty AABBDD bam files" if (scalar(@bam_AABBDD)<1);
print "InputOutputReport: Reading AABBDD bams: ".scalar(@bam_AABBDD)."\n" if ($debug or $verbose);
foreach (@bam_AABBDD) {
	print "---> $_\n" if ($debug or $verbose);
	unless (-s "$_.bai") {
		if (! &IndexBam($_, $path_samtools)) {
			die "InputOutputError: AABBDD BAM index error\n";
		}
	}
}



# read AABB bamfiles
my @bam_AABB=();
if ($file_bam_aabb=~/\.bam$/i) {
	push (@bam_AABB, $file_bam_aabb);
}
else {
	close BAMAABB if (defined fileno(BAMAABB));
	open (BAMAABB, "<$file_bam_aabb") || die "InputOutputError: can not open AABB $file_bam_aabb\n";
	while (my $line2=<BAMAABB>) {
		chomp $line2;
		if (defined $line2 and $line2 ne '' and -s $line2) {
			push (@bam_AABB, $line2);
		}
		else {
			die "InputOutputError: bam $line2 in AABB $file_bam_aabb not exist or empty\n";
		}
	}
	close BAMAABB;
}
die "InputOutputError: empty AABB bam files" if (scalar(@bam_AABB)<1);
print "InputOutputReport: Reading AABB bams: ".scalar(@bam_AABB)."\n" if ($debug or $verbose);
foreach (@bam_AABB) {
	print "---> $_\n" if ($debug or $verbose);
	unless (-s "$_.bai") {
		if (! &IndexBam($_, $path_samtools)) {
			die "InputOutputError: AABB BAM index error\n";
		}
	}
}

# read AA bam files
my @bam_AA=();
if ($file_bam_aa=~/\.bam$/i) {
	push (@bam_AA, $file_bam_aa);
}
else {
	close BAMAA if (defined fileno(BAMAA));
	open (BAMAA, "<$file_bam_aa") || die "InputOutputError: can not open AA $file_bam_aa\n";
	while (my $line3=<BAMAA>) {
		chomp $line3;
		if (defined $line3 and $line3 ne '' and -s $line3) {
			push (@bam_AA, $line3);
		}
		else {
			die "InputOutputError: bam $line3 in AA $file_bam_aa not exist or empty\n";
		}
	}
	close BAMAA;
}
die "InputOutputError: empty AA bam files" if (scalar(@bam_AA)<1);
print "InputOutputReport: Reading AA bams: ".scalar(@bam_AA)."\n" if ($debug or $verbose);
foreach (@bam_AA) {
	print "---> $_\n" if ($debug or $verbose);
	unless (-s "$_.bai") {
		if (! &IndexBam($_, $path_samtools)) {
			die "InputOutputError: AA BAM index error\n";
		}
	}
}

# read DD BAM files
my @bam_DD=();
if ($file_bam_dd=~/\.bam$/i) {
	push (@bam_DD, $file_bam_dd);
}
else {
	close BAMDD if (defined fileno(BAMDD));
	open (BAMDD, "<$file_bam_dd") || die "InputOutputError: can not open DD $file_bam_dd\n";
	while (my $line4=<BAMDD>) {
		chomp $line4;
		if (defined $line4 and $line4 ne '' and -s $line4) {
			push (@bam_DD, $line4);
		}
		else {
			die "InputOutputError: bam $line4 in DD $file_bam_dd not exist or empty\n";
		}
	}
	close BAMDD;
}
die "InputOutputError: empty DD bam files" if (scalar(@bam_DD)<1);
print "InputOutputInfo: Reading DD bams: ".scalar(@bam_DD)."\n" if ($debug or $verbose);
foreach (@bam_DD) {
	print "---> $_\n" if ($debug or $verbose);
	unless (-s "$_.bai") {
		if (! &IndexBam($_, $path_samtools)) {
			die "InputOutputError: DD BAM index error\n";
		}
	}
}

###input vcf
if (defined $file_vcf_aabbdd) {
	if (! -s and $file_vcf_aabbdd) {
		die "MainError: can not find AABBDD file: $file_vcf_aabbdd\n";
	}
	else {
		if ($file_vcf_aabbdd=~/\.vcf$/i) {
			if (! &BgzipVcf($file_vcf_aabbdd, "$file_vcf_aabbdd.gz", $path_bgzip)) {
				die "InputOutputError: AABBDD VCF bgzip\n";
			}
			elsif (! &IndexVcf("$file_vcf_aabbdd.gz", $path_tabix)) {
				die "InputOutputError: AABBDD VCF index\n";
			}
			else {
				$file_vcf_aabbdd.='.gz';
			}
		}
		elsif ($file_vcf_aabbdd=~/\.vcf\.gz$/i) {
			unless (-s "$file_vcf_aabbdd.tbi") {
				if (! &IndexVcf($file_vcf_aabbdd, $path_tabix)) {
					die "InputOutputError: AABBDD VCF index2\n";
				}
			}
		}
		else {
			die "InputOutputError: unknown AABBDD VCF format\n";
		}
	}
}

if (defined $file_vcf_aabb) {
	if (! -s $file_vcf_aabb) {
		die "MainError: can not find AABB file: $file_vcf_aabb\n";
	}
	else {
		if ($file_vcf_aabb=~/\.vcf$/i) {
			if (! &BgzipVcf($file_vcf_aabb, "$file_vcf_aabb.gz", $path_bgzip)) {
				die "InputOutputError: AABB VCF bgzip\n";
			}
			elsif (! &IndexVcf("$file_vcf_aabb.gz", $path_tabix)) {
				die "InputOutputError: AABB VCF index\n";
			}
			else {
				$file_vcf_aabb.='.gz';
			}
		}
		elsif ($file_vcf_aabb=~/\.vcf\.gz$/i) {
			unless (-s "$file_vcf_aabb.tbi") {
				if (! &IndexVcf($file_vcf_aabb, $path_tabix)) {
					die "InputOutputError: AABB VCF index2\n";
				}
			}
		}
		else {
			die "InputOutputError: unknown AABB VCF format\n";
		}
	}
}
if (defined $file_vcf_aa) {
	if (! -s $file_vcf_aa) {
		die "MainError: can not find AA file: $file_vcf_aa\n";
	}
	else {
		if ($file_vcf_aa=~/\.vcf$/i) {
			if (! &BgzipVcf($file_vcf_aa, "$file_vcf_aa.gz", $path_bgzip)) {
				die "InputOutputError: AA VCF bgzip\n";
			}
			elsif (! &IndexVcf("$file_vcf_aa.gz", $path_tabix)) {
				die "InputOutputError: AA VCF index\n";
			}
			else {
				$file_vcf_aa.='.gz';
			}
		}
		elsif ($file_vcf_aa=~/\.vcf\.gz$/i) {
			unless (-s "$file_vcf_aa.tbi") {
				if (! &IndexVcf($file_vcf_aa, $path_tabix)) {
					die "InputOutputError: AA VCF index2\n";
				}
			}
		}
		else {
			die "InputOutputError: unknown AA VCF format\n";
		}
	}
}
if (defined $file_vcf_dd) {
	if (! -s $file_vcf_dd) {
		die "MainError: can not find DD file: $file_vcf_dd\n";
	}
	else {
		if ($file_vcf_dd=~/\.vcf$/i) {
			if (! &BgzipVcf($file_vcf_dd, "$file_vcf_dd.gz", $path_bgzip)) {
				die "InputOutputError: DD VCF bgzip\n";
			}
			elsif (! &IndexVcf("$file_vcf_dd.gz", $path_tabix)) {
				die "InputOutputError: DD VCF index\n";
			}
			else {
				$file_vcf_dd.='.gz';
			}
		}
		elsif ($file_vcf_dd=~/\.vcf\.gz$/i) {
			unless (-s "$file_vcf_dd.tbi") {
				if (! &IndexVcf($file_vcf_dd, $path_tabix)) {
					die "InputOutputError: DD VCF index2\n";
				}
			}
		}
		else {
			die "InputOutputError: unknown DD VCF format\n";
		}
	}
}

if ($verbose) {
	print "InputOutputReport: VCF inputs: \n";
	if (defined $file_vcf_aabbdd) {print "AABBDD VCF1: $file_vcf_aabbdd\n";} else {print "AABBDD VCF1: null\n";}
	if (defined $file_vcf_aabb) {print "AABB VCF2: $file_vcf_aabb\n";} else {print "AABB VCF2: null\n";}
	if (defined $file_vcf_aa) {print "AA VCF3: $file_vcf_aa\n";} else {print "AA VCF3: null\n";}
	if (defined $file_vcf_dd) {print "DD VCF4: $file_vcf_dd\n";} else {print "DD VCF4: null\n";}
}



### Main ############################################################
my $RunDir=getcwd;
#Prepare temporary folder
if (! -d "$RunDir/AABBDD") {
	mkdir ("$RunDir/AABBDD", 0766) || die "(Main))Error: can not create folder $RunDir/AABBDD\n";
}
unlink glob "$RunDir/AABBDD/*";
if (! -d "$RunDir/AABB") {
	mkdir ("$RunDir/AABB", 0766) || die "(Main)Error: can not create folder AABB\n";
}
unlink glob "$RunDir/AABB/*";
if (! -d "$RunDir/AA") {
	mkdir ("$RunDir/AA", 0766) || die "(Main)Error: can not create folder AA\n";
}
unlink glob "$RunDir/AA/*";
if (! -d "$RunDir/DD") {
	mkdir ("$RunDir/DD", 0766) || die "(Main)Error: can not create folder DD\n";
}
unlink glob "$RunDir/DD/*";

# index fasta
my ($test_cdbfasta, $fasta_index)=&CdbFasta($reference, $path_cdbfasta);
if ($test_cdbfasta == $cmd_failure_code ) {
	die "(Main)Error: can not index fasta: $reference\n";
}
print "Fasta Index: $fasta_index\n";
# Load clustered fasta id list, 1 cluster perl line, saparated by space
close CLUSTER if (defined fileno(CLUSTER));
open (CLUSTER, "<$file_cluster") || die "(Main)Error: can not open file: $file_cluster\n";
# Output cluster succeeds or not. Format: 
close CLUSTERLOG if (defined fileno(CLUSTERLOG));
open (CLUSTERLOG, ">$cluster_log") || die "(Main)Error: can not write file: $cluster_log\n";
# output fpkm calculated
close FPKMLOG if (defined fileno(FPKMLOG));
open (FPKMLOG, ">$fpkm_log") || die "(Main)Error: can not write file: $fpkm_log\n";
print FPKMLOG "###ClusterID\n###UnigeneID\tAABBDD\tAABB\tAA\tDD\n";
# output allele log
close ALLELELOG if (defined fileno(ALLELELOG));
open (ALLELELOG, ">$allele_log") || die "(Main)Error: can not write file: $allele_log\n";
my $cluster_linenum=0;###for quick find which line/cluster has problem
my @cluster_seqids=();###
while (my $cluster_line=<CLUSTER>) {
	chomp $cluster_line;
	$cluster_linenum++;
	if ($cluster_linenum<$cluster_start) {
		next;
	}
	elsif ($cluster_linenum>=$cluster_start) {
		if (defined $cluster_end and $cluster_linenum>$cluster_end) {
			next;
		}
	}
	print STDOUT "\n\n\n\n\n##### Prcessing Cluster $cluster_linenum ###\n";
	print STDERR "\n\n\n\n\n##### Prcessing Cluster $cluster_linenum ###\n";
	@cluster_seqids=(); ###Empty this in case of abnormal duplicates
	@cluster_seqids=split(/\s+/, $cluster_line);
##COMMENT: Check if empty line
	if (scalar(@cluster_seqids)<1) {
		print STDERR "Main1Warnings: line $cluster_linenum in $file_cluster ignored as empty\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tNumIDinCluster<1\n";
		next;
	}
##COMMENT: extract fasta reference
	if ( ! -d "$RunDir/Clust$cluster_linenum") {
		unless (mkdir ("$RunDir/Clust$cluster_linenum", 0766)) {
			print STDERR "(Main2)Error: create folder $RunDir/Clust$cluster_linenum\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tMkdirCluster\n";
			next;
		}
	}
	unless ( -s "$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa" and "$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa.fai") {
		if (&CdbYank($fasta_index, "$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, $path_cdbyank) == $cmd_failure_code) {
			print STDERR "(Main2)Error: failed to extract ClusterID $cluster_linenum: @cluster_seqids\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tCdbYank\n";
			next;
		}
		if (&IndexFasta("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", $path_samtools) == $cmd_failure_code) {
			print STDERR "(Main2)Error: failed to index ClusterID $cluster_linenum: @cluster_seqids\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tIndexFasta\n";
			next;
		}
	}
##COMMENT: Empty folder
	unlink glob "$RunDir/AABBDD/*";
	unlink glob "$RunDir/AABB/*";
	unlink glob "$RunDir/AA/*";
	unlink glob "$RunDir/DD/*";
##COMMENT: Extract bam from @bam_AABBDD, @bam_AABB, @bam_AA, @bam_DD
	my $test_aabbdd_exporessed=0;
	my @cluster_bam_files=();
	my $aabbdd_vcf_obj=[];
	###AABBDD BAM
	if (&GetBam(\@bam_AABBDD, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.nogroup.bam")) {
		print STDERR "(Main3)Error: AABBDD BAM extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGetBamAABBDD\n";
		next;
	}
	### AABBDD FPKM
	 (my $TEtestcode_aabbdd, $test_aabbdd_exporessed)=&TestExpressed("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.nogroup.bam", $totalreads[0]);
	if ($TEtestcode_aabbdd) {
		print STDERR "(Main3)Error: Calculate AABBDD FPKM error\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tTestExpressAABBDD\n";
		next;
	}
	if ($test_aabbdd_exporessed ==0) {
		print STDERR "(Main3)Error: AABBDD FPKM 0\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tAABBDD_FPKM_0\n";
		next;
	}
	### AABBDD VCF
	if (defined $file_vcf_aabbdd and -s $file_vcf_aabbdd) {
		if (! &ExtractVcf($file_vcf_aabbdd, $cluster_line, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz")) {
			print STDERR "(Main3)Error: AABBDD VCF extraction failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tVcfExtractAABBDD\n";
			next;
		}
	}
	else {
		unless (-s "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz" and -s "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz.tbi") {
			if (&RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.nogroup.bam", 3, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", 0)) {
				print STDERR "(Main3)Error: freebayes running AABBDD failed\n";
				print CLUSTERLOG $cluster_line."\tFail\t1\tFreeBayesAABBDD\n";
				next;
			}
		}
	}
	unless (-s "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf") {
		if (! &exec_cmd_return("gunzip -c $RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz > $RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf")) {
			print STDERR "(Main3)Error: gunzip uncompress $RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFreeBayesAABBDDgunzip\n";
			next;
		}
		elsif (! -s "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf") {
			print STDERR "(Main3)Error: gunzip output $RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFreeBayesAABBDDgunzip_output\n";
			next;
		}
	}
	(my $test_readvcf_aabbdd, $aabbdd_vcf_obj)=&ReadVcf("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", $vcfqual);
	if (! $test_readvcf_aabbdd) {
		print STDERR "(Main3)Error: Reading AABBDD VCF failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\treadvcfAABBDD\n";
		next;
	}
	my $test_AABBDD_have_variation=0;
	foreach (keys %{$aabbdd_vcf_obj}) {
		if (scalar(keys %{${$aabbdd_vcf_obj}{$_}} >0)) {
			$test_AABBDD_have_variation=1;
			last;
		}
	}
	if ($debug) {
		print "$cluster_linenum: no vaiations detected, direct assembly\n";
	}




	### AABB BAM
	my $test_aabb_expressed=0;
	my $aabb_vcf_obj=[];
###BLOCKAABB
	BLOCKAABB: {
	if (&GetBam(\@bam_AABB, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.nogroup.bam")) {
		print STDERR "(Main4)Error: AABB BAM extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGetBamAABB\n";
		next;
	}
	### AABB FPKM
	(my $TEtestcode_aabb, $test_aabb_expressed)=&TestExpressed("$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.nogroup.bam", $totalreads[1]);
	if ($TEtestcode_aabb) {
		print STDERR "(Main4)Error: Calculate AABB FPKM error\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tTestExpressAABB\n";
		next;
	}
	if ($test_aabb_expressed ==1) {
		#run freebayes for AABB
		if (defined $file_vcf_aabb and -s $file_vcf_aabb) {
			if (! &ExtractVcf($file_vcf_aabb, $cluster_line, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz")) {
				print STDERR "(Main4)Error: AABB VCF extraction failed\n";
				print CLUSTERLOG $cluster_line."\tFail\t1\tVcfExtractAABB\n";
				next;
			}
		}
		else {
			unless (-s "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz" and -s "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz.tbi") {
				if (&RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.nogroup.bam", 2, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz")) {
					print STDERR "(Main4)Error: freebayes running AABB failed\n";
					print CLUSTERLOG $cluster_line."\tFail\t1\tFreeBayesAABB\n";
					next;
				}
			}
		}
		(my $test_readvcf_aabb, $aabb_vcf_obj)=&ReadVcf("$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz", $vcfqual);
		if (! $test_readvcf_aabb) {
			print STDERR "(Main4)Error: Reading AABB VCF failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\treadvcfAABB\n";
			next;
		}
	}}###BLOCKAABB
	
	
	
	
	### AA bam
	my $test_aa_expressed=0;
	my $aa_vcf_obj=[];
###BLOCKAA
	BLOCKAA: {
	if (&GetBam(\@bam_AA, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.nogroup.bam")) {
		print STDERR "(Main5)Error: AA BAM extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGetBamAA\n";
		next;
	}
	### AA FPKM
	(my $TEtestcode_aa, $test_aa_expressed)=&TestExpressed("$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.nogroup.bam", $totalreads[2]);
	if ($TEtestcode_aa) {
		print STDERR "(Main5)Error: Calculate AA FPKM error\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tTestExpressAA\n";
		next;
	}
	if ($test_aa_expressed ==1) {
		#run_freebayes for AA
		if (defined $file_vcf_aa and -s $file_vcf_aa) {
			if (! &ExtractVcf($file_vcf_aa, $cluster_line, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz")) {
				print STDERR "(Main5)Error: AA VCF extraction failed\n";
				print CLUSTERLOG $cluster_line."\tFail\t1\tVcfExtractAA\n";
				next;
			}
		}
		else {
			unless (-s "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz" and -s "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz.tbi") {
				if (&RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.nogroup.bam", 2, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz")) {
					print STDERR "(Main5)Error: freebayes running AA failed\n";
					print CLUSTERLOG $cluster_line."\tFail\t1\tFreeBayesAA\n";
					next;
				}
			}
		}
		(my $test_readvcf_aa, $aa_vcf_obj)=&ReadVcf("$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz", $vcfqual);
		if (! $test_readvcf_aa) {
			print STDERR "(Main5)Error: Reading AA VCF failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\treadvcfAA\n";
			next;
		}
	}}###BLOCKAA
	
	
	
	### DD bam
	my $test_dd_expressed=0;
	my $dd_vcf_obj=[];
###BLOCKDD
	BLOCKDD: {
	if (&GetBam(\@bam_DD, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.nogroup.bam")) {
		print STDERR "(Main6)Error: DD BAM extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGetBamDD\n";
		next;
	}
	(my $TEtestcode_dd, $test_dd_expressed)=&TestExpressed("$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.nogroup.bam", $totalreads[3]);
	if ($TEtestcode_dd) {
		print STDERR "(Main6)Error: Calculate DD FPKM error\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tTestExpressDD\n";
		next;
	}
	if ($test_dd_expressed ==1) {
	#run_freebayes for DD
		if (defined $file_vcf_dd and -s $file_vcf_dd) {
			if (! &ExtractVcf($file_vcf_dd, $cluster_line, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz")) {
				print STDERR "(Main6)Error: DD VCF extraction failed\n";
				print CLUSTERLOG $cluster_line."\tFail\t1\tVcfExtractDD\n";
				next;
			}
		}
		else {
			unless (-s "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz" and -s "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz.tbi") {
				if (&RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.nogroup.bam", 2, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz")) {
					print STDERR "(Main6)Error: freebayes running DD failed\n";
					print CLUSTERLOG $cluster_line."\tFail\t1\tFreeBayesDD\n";
					next;
				}
			}
		}
		(my $test_readvcf_dd, $dd_vcf_obj)=&ReadVcf("$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz", $vcfqual);
		if (! $test_readvcf_dd) {
			print STDERR "(Main6)Error: Reading DD VCF failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\treadvcfDD\n";
			next;
		}
	}}###BLOCKDD
	


##COMMENT: Group alleles into subgenome
	my $test_bb_expressed=0;
	my ($test_readbam, $aabbdd_bam_obj)=&ReadSam("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.bam", "Clust$cluster_linenum/ref.$cluster_linenum.fa", 1);
	if (! $test_readbam) {
		print STDERR "(Main7)Error: Reading AABBDD BAM failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\treadsamAABBDD\n";
		next;
	}
	(my $test_groupvcf, my $rnaseq_allele_hashindex, my $fixed_allele_hashindex, my $allele2readids_hashindex, my $readid2genome_hashindex, my $sharedReadid_hashindex, $test_bb_expressed)=&GroupVcf($aabbdd_bam_obj, $aabbdd_vcf_obj, $aabb_vcf_obj, $aa_vcf_obj, $dd_vcf_obj);
	if ($test_groupvcf) {
		print STDERR "(Main7)Error: GroupVcf SUB failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGroupVcf\n";
		next;
	}
	my %FixedAllele_1rnaseq=%{$rnaseq_allele_hashindex}; undef $rnaseq_allele_hashindex;
	my %FixedAllele_2genome=%{$fixed_allele_hashindex}; undef $fixed_allele_hashindex;



##COMMENT: decide ploidy here
	my $ploidy=$test_aa_expressed+$test_dd_expressed+$test_bb_expressed;
	my $geno_ploidy='';
	$geno_ploidy.='A' if ($test_aa_expressed>0);
	$geno_ploidy.='B' if ($test_bb_expressed>0);
	$geno_ploidy.='D' if ($test_dd_expressed>0);
	$geno_ploidy='ABD' if ($ploidy==0 or $geno_ploidy eq '');



##COMMENT: Fill in blanks and prepare vcf for
	my $test_run_HapCompass=0;
	my %FixedAllele_3fillin=();
	my %fixedGeno_4fillin=();
	if ($ploidy==1) {
		$test_run_HapCompass=0;
		%FixedAllele_3fillin=%FixedAllele_2genome;
	}
	elsif ($ploidy==3 or $ploidy==2 or $ploidy==0) {
		(my $test_fiv, $test_run_HapCompass, my $fixedgeno_hashindex, my $fillinAlleles_hashindex)=&FillInVariations(\%FixedAllele_2genome,$aabbdd_vcf_obj, $geno_ploidy);
		if ($test_fiv==0) {###Fillin succeeds
			%FixedAllele_3fillin=%{$fillinAlleles_hashindex}; undef $fillinAlleles_hashindex;
			%fixedGeno_4fillin=%{$fixedgeno_hashindex}; undef $fixedgeno_hashindex;
		}
		elsif ($test_fiv==1) {
			print STDERR "(Main8)Error: FillIn1 SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin1\n";
			next;
		}
		elsif ($test_fiv==2) {
			print STDERR "(Main8)Error: FillIn2 SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin2\n";
			next;
		}
		elsif ($test_fiv==1) {
			print STDERR "(Main8)Error: FillIn3 SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin3\n";
			next;
		}
		else {
			print STDERR "(Main8)Error: FillIn SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin4\n";
			next;
		}
	}
	else {
		print STDERR "(Main8)Error: unknown ploidy\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tPloidy\n";
		next;
	}
##COMMENT:
	if (! &RunHaptree("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.haptree.vcf", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.reads", \%fixedGeno_4fillin, $allele2readids_hashindex, $path_haptree)) {
		print STDERR "(Main9)Error: HapTree running\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tHaptree\n";
	}
	foreach my $chrom (keys %{$aabbdd_vcf_obj}) {
		print "(Main9)Test: AlleleTesting: chromosome $chrom\n" if ($debug);
		foreach my $posit (keys %{${$aabbdd_vcf_obj}{$chrom}}) {
			print "(Main9)Test: AlleleTesting: positions $posit\n" if ($debug);
			my $printline=$cluster_linenum."\t".$chrom."\t".$posit."\t".${${${$aabbdd_vcf_obj}{$chrom}}{$posit}}[2]."\t";
			#%FixedAllele_1rnaseq
			$printline.=(exists ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'A'}) ? ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'A'} : 'xxxxx';
			$printline.="\t";
			$printline.=(exists ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'B'}) ? ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'B'} : 'xxxxx';
			$printline.="\t";
			$printline.=(exists ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'D'}) ? ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'D'} : 'xxxxx';
			$printline.="\t";
			#%FixedAllele_2genome
			$printline.=(exists ${${$FixedAllele_2genome{$chrom}}{$posit}}{'A'}) ? ${${$FixedAllele_2genome{$chrom}}{$posit}}{'A'} : 'xxxxx';
			$printline.="\t";
			$printline.=(exists ${${$FixedAllele_2genome{$chrom}}{$posit}}{'B'}) ? ${${$FixedAllele_2genome{$chrom}}{$posit}}{'B'} : 'xxxxx';
			$printline.="\t";
			$printline.=(exists ${${$FixedAllele_2genome{$chrom}}{$posit}}{'D'}) ? ${${$FixedAllele_2genome{$chrom}}{$posit}}{'D'} : 'xxxxx';
			$printline.="\t";
			#%FixedAllele_3fillin
			$printline.=(exists ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'A'}) ? ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'A'} : 'xxxxx';
			$printline.="\t";
			$printline.=(exists ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'B'}) ? ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'B'} : 'xxxxx';
			$printline.="\t";
			$printline.=(exists ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'D'}) ? ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'D'} : 'xxxxx';
			$printline.="\t";
			#%fixedGeno_4fillin
			my $fillin_geno=(exists ${$fixedGeno_4fillin{$chrom}}{$posit}) ? ${$fixedGeno_4fillin{$chrom}}{$posit} : 'xxxxx';
			$printline.=$fillin_geno."\t";
			#$test_run_HapCompass
			$printline.=$test_run_HapCompass."\n";
			print ALLELELOG $printline;
		}
	}

	
	
	
	
	
	
	print CLUSTERLOG $cluster_line."\tSucceed\t0\t0\n";
}
close CLUSTER;
close CLUSTERLOG;
close FPKMLOG;
close ALLELELOG;



#####################################################################
###                         sub functions                         ###
#####################################################################



### extract batch of sequence alignment from Bams, and index
### &GetBam(bam_arr, sequence_arr, output)
### Global:$cluster_linenum, $path_samtools
### Dependency:
### Note: 
sub GetBam {
	my ($EBfiles_bam_index, $EBseq_obj, $EBoutput, $EBoutnogroup)=@_;
	
	my $EBsubinfo='SUB(GetBam)';
	
	return 1 if (scalar(@{$EBseq_obj})<1);###return error if no sequences
	return 0 if (-s $EBoutput and -s "$EBoutput.bai" and -s $EBoutnogroup and -s "$EBoutnogroup.bai");
	
##COMMENT: extract bam from each
	my $EBsam_seq=join (' ', @{$EBseq_obj});
	my $EBi=0;
	my @EBtemp_bams=();
#	print STDERR "${EBsubinfo}Info: sequence to be extracted: $EBsam_seq\n";###for test###
	foreach my $EBind_bam (@{$EBfiles_bam_index}) {
		$EBi++;
		my $EBbamout='';
		if (scalar(@{$EBfiles_bam_index}) ==1) {
			$EBbamout="$EBoutput.multihead.bam";
		}
		else {
			$EBbamout="$EBoutput.$EBi.bam";
		}
		my $EBcmd="$path_samtools view -bh $EBind_bam $EBsam_seq > $EBbamout";
		if (! &exec_cmd_return($EBcmd)) {
			print STDERR "${EBsubinfo}Error: bam extract $EBsam_seq (Cluster: $cluster_linenum) from $EBind_bam failed\n";
			return 1;
		}
		else {
			push (@EBtemp_bams, "$EBoutput.$EBi.bam");
		}
	}
##COMMENT: merge bams 
	unless (scalar(@{$EBfiles_bam_index}) ==1) {
		my $EBmerge_bams=join(' ', @EBtemp_bams);
		my $EBcmd2="$path_samtools merge $EBoutput.multihead.bam $EBmerge_bams";
		if (! &exec_cmd_return($EBcmd2)) {
			print STDERR "${EBsubinfo}Error: bam merge (Cluster: $cluster_linenum) from $EBmerge_bams failed\n";
			return 1;
		}
		else {
			if (-s "$EBoutput.multihead.bam") {
				map {unlink($_)} @EBtemp_bams;###delete temporary files ### For test ###
			}
			else{
				print STDERR "${EBsubinfo}Error: bam merge (Cluster: $cluster_linenum) from $EBmerge_bams not found\n";
				return 1;
			}
		}
	}
	if (! &SamCleanHeader("$EBoutput.multihead.bam", "$EBoutput.withRG.bam", "$EBoutput.noRG.bam", $path_samtools)) {
		print STDERR "${EBsubinfo}Error: bam header clean (Cluster: $cluster_linenum) failed\n";
		return 1;
	}
	else {
		if ( -s "$EBoutput.multihead.bam" and -s "$EBoutput.withRG.bam") {
			unlink ("$EBoutput.multihead.bam") if (-s "$EBoutput.multihead.bam");
			if (! &SortBam("$EBoutput.withRG.bam", "$EBoutput.withRG.st.bam", $path_samtools)) {
				print STDERR "${EBsubinfo}Error: BAM sort running +RG error: $EBoutput.withRG.bam\n";
				return 1;
			}
			elsif (! -s "$EBoutput.withRG.st.bam") {
				print STDERR "${EBsubinfo}Error: BAM sort output +RG error: $EBoutput.withRG.bam\n";
				return 1;
			}
			else {
				rename ("$EBoutput.withRG.st.bam", $EBoutput);
				if (! &IndexBam("$EBoutput", $path_samtools)) {
					print STDERR "${EBsubinfo}Error: BAM index running +RG error: $EBoutput\n";
					return 1;
				}
			}
			if (! &SortBam("$EBoutput.noRG.bam", "$EBoutput.noRG.st.bam", $path_samtools)) {
				print STDERR "${EBsubinfo}Error: BAM sort running +RG error: $EBoutput.noRG.bam\n";
				return 1;
			}
			elsif (! -s "$EBoutput.noRG.st.bam") {
				print STDERR "${EBsubinfo}Error: BAM sort output +RG error: $EBoutput.noRG.bam\n";
				return 1;
			}
			else {
				rename ("$EBoutput.noRG.st.bam", $EBoutnogroup);
				if (! &IndexBam($EBoutnogroup, $path_samtools)) {
					print STDERR "${EBsubinfo}Error: BAM index running +RG error: $EBoutnogroup\n";
					return 1;
				}
			}
		}
		else{
			print STDERR "${EBsubinfo}Error: bam headerclean output (Cluster: $cluster_linenum) not found\n";
			return 1;
		}
	}
	return 0;
}



### Evaluate Ploidy by diploid
### &TestExpressed(\@BAMs, \@total_reads)
### Global: @totalreads, $path_samtools
### Dependency:
### Note:
sub TestExpressed_old {
	my ($TEbam_arrindex, $TEseqids_arrindex)=@_;
	
	my $TEsubinfo='SUB(EvaluatePloidy)';
	my @TEtest_expressed=();
	
	if (scalar(@{$TEbam_arrindex}) <1 or scalar(@{$TEbam_arrindex}) != scalar(@totalreads)) {
		print STDERR "${TEsubinfo}Error: empty or invalid bam files\n";
		return 1;
	}
	
	for (my $TEnum=0; $TEnum<scalar(@{$TEbam_arrindex}); $TEnum++) {
		my ($TEtestcode, $TEfpkm_hashindex)=&CalcFPKM(${$TEbam_arrindex}[$TEnum], $totalreads[$TEnum], 0, $path_samtools);
		if (! $TEtestcode) {
			print STDERR "${TEsubinfo}Error: Calculate FPKM error\n";
			return 1;
		}
		my $EPexpressed=0;
		foreach (@{$TEseqids_arrindex}) {
			if (exists ${$TEfpkm_hashindex}{$_} and ${$TEfpkm_hashindex}{$_}>=$minimum_fpkm) {
				$EPexpressed=1;
			}
		}
		push (@TEtest_expressed, $EPexpressed);
	}
	
	return (0, \@TEtest_expressed);
}



### Test Expressed or not
### &TestExpressed(BAMin, $total_reads)
### Global: $path_samtools, $minimum_fpkm
### Dependency:
### Note:
sub TestExpressed {
	my ($TEbamfile, $TEtotalreads)=@_;
	
	my $TEsubinfo='SUB(TestExpressed)';
	
	my ($TEtestcode, $TEfpkm_hashindex)=&CalcFPKM($TEbamfile, $TEtotalreads, 0, $path_samtools);
	if (! $TEtestcode) {
		print STDERR "${TEsubinfo}Error: Calculate FPKM error\n";
		return 1;
	}
	my $EPexpressed=0;
	foreach (keys %{$TEfpkm_hashindex}) {
		if (exists ${$TEfpkm_hashindex}{$_} and exists ${${$TEfpkm_hashindex}{$_}}{'unknown'} and ${${$TEfpkm_hashindex}{$_}}{'unknown'}>=$minimum_fpkm) {
			$EPexpressed=1;
			last;
		}
	}
	
	return (0, $EPexpressed);
}

###RunFreebayes, return merged vcf and raw AABBDD vcf
###&RunFreebayes(ref_fasta, seqids_arrindex, bam, ploidy, output.vcf.gz, guided_vcf)
###Global: $freebayes_min_coverage, $freebayes_min_alternative_count, $min_mapq,$express_frag_len_mean, $express_frag_len_stddev, $path_freebayesparallel, $numthreads, $freebayes_fastabin, $path_vcfsort
###Dependancy: 
sub RunFreebayes {
	my ($RFfile_reference, $RFseqids_arrindex, $RFfile_bam, $RFploidy, $RFoutput, $RFguide_vcf)=@_;
	unlink $RFoutput if (-e $RFoutput);
	unlink "$RFoutput.tbi" if (-e "$RFoutput.tbi");
	my $RFcmd='';
	my @RFseqids=@{$RFseqids_arrindex};
	if (scalar(@RFseqids)<1) {
		print "SUB(RunFreebayes)Error: no seqs for freebayes";
		return 1;
	}
	if (! defined $RFguide_vcf) {
		print "SUB(RunFreebayes)Error: undefined guide VCF for freebayes";
		return 1;
	}
	unless (-d "$RunDir/Freebayes") {
		unless (mkdir ("$RunDir/Freebayes", 0766)) {
			print STDERR "SUB(RunFreebayes)Error: can not create folder Freebayes\n";
			return 1;
		}
	}
	unlink glob "$RunDir/Freebayes/*";
	my $RFcmd_freebayes=$path_freebayes;
	if ($freebayes_parallel) {
		my $RFtest_regioncmd=&CreateFastaRegion("$RFfile_reference.fai", $freebayes_fastabin, "$RFfile_reference.$freebayes_fastabin.region");
		if (! $RFtest_regioncmd) {
			print STDERR "SUB(RunFreebayes)Error: CreateFastaRegion for $RFfile_reference\n";
			return 1;
		}
		else {
			$RFcmd_freebayes=$path_freebayesparallel." $RFfile_reference.$freebayes_fastabin.region $numthreads ";
		}
	}
	if ($RFguide_vcf eq '0') {
		$RFcmd="$RFcmd_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy --pooled-discrete --min-mapping-quality $min_mapq --genotype-qualities $RFfile_bam | $path_vcfsort | $path_bgzip > $RFoutput";
		if (! &exec_cmd_return($RFcmd)) {
			print STDERR "SUB(RunFreebayes)Error1:  freebayes running error\n";
			return 1;
		}
		if (! -s $RFoutput) {
			print STDERR "SUB(RunFreebayes)Error1: freebayes output error\n";
			return 1;
		}
	}
	elsif (-s $RFguide_vcf) {
		if (scalar(@RFseqids) ==1) {
			$RFcmd="$RFcmd_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy --pooled-discrete --min-mapping-quality $min_mapq --variant-input $RFguide_vcf --only-use-input-alleles --region $RFseqids[0] --genotype-qualities  $RFfile_bam | $path_vcfsort | $path_bgzip > $RFoutput";
			if (! &exec_cmd_return($RFcmd)) {
				print STDERR "SUB(RunFreebayes)Error2:  freebayes running error\n";
				return 1;
			}
			if (! -s $RFoutput) {
				print STDERR "SUB(RunFreebayes)Error2: freebayes output error\n";
				return 1;
			}
		}
		else {#freebayes can not tolerate more than 1 seqs using --variant-input, so need to separate the seqids
			my @RFtemp_vcf=();
			foreach my $RFind_seq (@RFseqids) {
				$RFcmd="$RFcmd_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy --pooled-discrete --min-mapping-quality $min_mapq --variant-input $RFguide_vcf --only-use-input-alleles --region $RFind_seq --genotype-qualities  $RFfile_bam | $path_vcfsort | $path_bgzip > $RunDir/Freebayes/$RFind_seq.vcf.gz";
				if (! &exec_cmd_return($RFcmd)) {
					print STDERR "SUB(RunFreebayes)Error3:  freebayes running error\n";
					return 1;
				}
				if (! -s "$RunDir/Freebayes/$RFind_seq.vcf.gz") {
					print STDERR "SUB(RunFreebayes)Error3: freebayes output error\n";
					return 1;
				}
				push (@RFtemp_vcf, "$RunDir/Freebayes/$RFind_seq.vcf.gz");
			}
			$RFcmd="$path_vcfconcat ".join(' ', @RFtemp_vcf)." | $path_bgzip > $RFoutput";
			if (! &exec_cmd_return($RFcmd)) {
				print STDERR "SUB(RunFreebayes)Error4:  freebayes running error\n";
				return 1;
			}
			if (! -s $RFoutput) {
				print STDERR "SUB(RunFreebayes)Error4: freebayes output error\n";
				return 1;
			}
			unlink (@RFtemp_vcf);###delete temporary files
		}
	}
	else {
		print STDERR "SUB(RunFreebayes)Error2:  unknown Freebayes guide VCF error\n";
		return 1;
	}
	my $RFcmd2="$path_tabix -p vcf $RFoutput";
	if (! &exec_cmd_return($RFcmd2)) {
		print STDERR "SUB(RunFreebayes)Error: freebayes index error for $RFoutput\n";
		return 1;
	}
	elsif (! -e "$RFoutput.tbi") {
		print STDERR "SUB(RunFreebayes)Error: freebayes index non-exist for $RFoutput\n";
		return 1;
	}
	else {
		return 0;
	}
}



###Asigning SNP allele 
###&AssignVariationAllele()
###Global: $debug, $verbose
###Dependancy: &ExtractAllele
sub AssignVariationAllele {
	my ($AVAaabbdd_vcf_index, $AVAaabb_vcf_index, $AVAaa_vcf_index, $AVAdd_vcf_index)=@_;
	my %AVAaabbdd_vcf_obj=%{$AVAaabbdd_vcf_index}; undef $AVAaabbdd_vcf_index;
	my %AVAaabb_vcf_obj=%{$AVAaabb_vcf_index}; undef $AVAaabbdd_vcf_index;
	my %AVAaa_vcf_obj=%{$AVAaa_vcf_index};undef $AVAaa_vcf_index;
	my %AVAdd_vcf_obj=%{$AVAdd_vcf_index}; undef $AVAdd_vcf_index;
	my $AVAtest_AA_expressed=1;
	$AVAtest_AA_expressed=0 if (scalar(keys %AVAaa_vcf_obj)==0);
	my $AVAtest_DD_expressed=1;
	$AVAtest_DD_expressed=0 if (scalar(keys %AVAdd_vcf_obj)==0);
	my $AVAtest_bb_expressed=0;
#Format: %AVAgenome2allele=(chr => (pos => ('A' => allele, 'B' => Allele, 'D' => allele)))
	my %AVAgenome2allele=();
#Format: %AVAallele2genome=(chr => (pos => (allele => 'A', allele => 'B', allele => 'D')))
	my %AVAallele2genome=();
	foreach my $AVAchrom (keys %AVAaabbdd_vcf_obj) {
		foreach my $AVApos (sort {$a<=>$b} keys %{$AVAaabbdd_vcf_obj{$AVAchrom}}) {
			print "SUB(AssignVariationAllele)Test: Chr: $AVAchrom\tPos: $AVApos\t Arr: @{${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}\n" if ($debug or $verbose);
			my ($AVAaabbdd_genohash_index, $AVAaabbdd_genoarr_index)=&ExtractAllele(${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[2], 1);
			if (scalar(@$AVAaabbdd_genoarr_index) ==1) {
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}=${$AVAaabbdd_genoarr_index}[0];
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}=${$AVAaabbdd_genoarr_index}[0];
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}=${$AVAaabbdd_genoarr_index}[0];
				${${$AVAallele2genome{$AVAchrom}}{$AVApos}}{${$AVAaabbdd_genoarr_index}[0]}='ABD';
			}
			else {
				my $AVAaa_allele='';
				my $AVAbb_allele='';
				my $AVAdd_allele='';
###COMMENT: defined dd allele
				if (exists ${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}) {
					my ($AVAdd_genohash_index, $AVAdd_genoarr_index)=&ExtractAllele(${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[2], 1);
					if (scalar(@{$AVAdd_genoarr_index})==1 and ${$AVAdd_genoarr_index}[0]=~/^\d+$/) {
						$AVAdd_allele=${$AVAdd_genoarr_index}[0];
						if (exists ${$AVAaabbdd_genohash_index}{$AVAdd_allele} and ${$AVAdd_genohash_index}{$AVAdd_allele} eq ${$AVAaabbdd_genohash_index}{$AVAdd_allele}) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}=$AVAdd_allele;
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{$AVAdd_allele}.='D';
						}
						else {
							print STDERR "SUB(AssignVariationAllele)Warning: non AABBDD allele $AVAdd_allele in DD at $AVAchrom:$AVApos\n";
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'});
##COMMENT: defined AA allele
				if (exists ${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}) {
					my ($AVAaa_genohash_index, $AVAaa_genoarr_index)=&ExtractAllele(${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[2], 1);
					if (scalar(@$AVAaa_genoarr_index)==1 and $AVAaa_allele=~/^\d+$/) {
						$AVAaa_allele=${$AVAaa_genoarr_index}[0];
						if (exists ${$AVAaabbdd_genohash_index}{$AVAaa_allele} and ${$AVAaa_genohash_index}{$AVAaa_allele} eq ${$AVAaabbdd_genohash_index}{$AVAaa_allele}) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}=$AVAaa_allele;
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{$AVAaa_allele}.='A';
						}
						else {
							print STDERR "SUB(AssignVariationAllele)Warning: non AABBDD allele $AVAaa_allele in AA at $AVAchrom:$AVApos\n";
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'});
##COMMENT: defined B allele---Difficult
				if ($AVAtest_AA_expressed==1 and $AVAaa_allele=~/^\d$/ and exists ${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}) {
					my ($AVAaabb_genohash_index, $AVAaabb_genoarr_index)=&ExtractAllele(${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[2], 1);
					if (scalar(@{$AVAaabb_genoarr_index})==2 and exists ${$AVAaabb_genohash_index}{$AVAaa_allele}) {
						foreach (@{$AVAaabb_genoarr_index}) {
							$AVAbb_allele=$_ unless ($AVAaa_allele eq $_);
						}
						if (exists ${$AVAaabbdd_genohash_index}{$AVAbb_allele} and ${$AVAaabb_genohash_index}{$AVAbb_allele} eq ${$AVAaabbdd_genohash_index}{$AVAbb_allele} and $AVAdd_allele ne $AVAbb_allele) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}=$AVAbb_allele;
							$AVAtest_bb_expressed=1;
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'});
##COMMENT: output debug info
				if ($debug) {
					print "chr: $AVAchrom\tPos: $AVApos\tRef: ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[0]\tVar: ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[1]\tAABBDDgeno ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[2]\tAABBgeno: ";
					if (exists ${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}) {
						print ${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tAA: ";
					if (exists ${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}) {
						print ${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tDD: ";
					if (exists ${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}) {
						print ${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tRNAphased: ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}/${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}/${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}\n";
				}
			}
		}
	}
	return (\%AVAgenome2allele, $AVAtest_bb_expressed);
}



###Extract allele based VCF ref, var, geno
###&ExtractAllele(ref, var, geno)
###Example: &ExtractAllele('AA', 'TT,GG', 0/1/1, 1/2)
###Global: 
###Dependency:
###Note: 
sub ExtractAllele {
	my ($EAref, $EAvar, $EAgeno, $EAreturn_code)=@_;
#Format: %EAgeno_hash=(0 => 'AA', 1 => 'TT', 2 => 'GG')
	my %EAgeno_hash=();
	$EAgeno_hash{0}=$EAref;
	my @EAvar_arr=split(/$geno_delimiter/,$EAvar);
	for (my $EAi=0; $EAi<scalar(@EAvar_arr);$EAi++) {
		$EAgeno_hash{$EAi+1}=$EAvar_arr[$EAi];
	}
#Format: %EAgeno2_hash=(0 => count, 1 => count, ...)
	my %EAgeno2_hash=();
#Format: @EAgeno_arr=(0, 1, 2)
	my @EAgeno_arr=();
	if ($EAgeno eq '.') {
		push (@EAgeno_arr, '?');
	}
	elsif ($EAgeno =~m/\//) {
		@EAgeno_arr=split(/\//, $EAgeno);
		foreach (@EAgeno_arr) {
			if ($_=~m/^\d$/) {
				$EAgeno2_hash{$_}++;
			}
			else {
				print STDERR "SUB(ExtractAllele)Error: non number genotype\n";
			}
		}
	}
	else {
		print STDERR "SUB(ExtractAllele)Error: unknown genotype\n";
	}
	my @EAgeno2_arr=sort {$a<=>$b} keys %EAgeno2_hash;
	if ($EAreturn_code ==1 ) {
		return (\%EAgeno_hash, \@EAgeno2_arr);
	}
	elsif ($EAreturn_code ==2) {
		return \@EAgeno_arr;
	}
	elsif ($EAreturn_code ==3) {
		return (\%EAgeno2_hash, \@EAgeno_arr);
	}
}



### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependancy:
###Note: 1=Bio::DB::Sam objective; 2=reference name array, 3=to be defined
sub ReadSam {
	my ($RSbam_file, $RSref_file, $RSret_code)=@_;
	
	my $RSsubinfo='SUB(ReadSam)';
	
	unless ( defined $RSbam_file and -s $RSbam_file) {
		print STDERR "${RSsubinfo}Error: invalid BAM file\n";
		return 0;
	}
	unless (defined $RSref_file and -s $RSref_file) {
		print STDERR "${RSsubinfo}Error: invalid Fasta file\n";
		return 0;
	}
	unless (defined $RSret_code and $RSret_code=~/^[12]{1}$/) {
		print STDERR "${RSsubinfo}Error: invalid return code\n";
		return 0;
	}
	
	my $RSsam_obj=Bio::DB::Sam->new(  -bam => "$RSbam_file", 
							-fasta => "$RSref_file",
							-expand_flags => 1,
							-split_splices => 1,
							-autoindex => 1,
						);
	if ($RSret_code==1) {
		return (1, $RSsam_obj);
	}
	elsif ($RSret_code==2) {
		my @arr=$RSsam_obj->seq_ids;
		return (1, @arr);
	}
	elsif ($RSret_code==3) {
		##do sth
	}
}



###group SNPs based on read length
###&GroupVcf($ReadSam_obj, AABBDD_Vcf_obj, AABB_vcf_obj, AA_vcf_obj, DD_vcf_obj)
###Global: $debug, $verbose
###Dependancy:&ReadVcf, &ReadVariantType, $min_mapq, $bam_genome_tag
sub GroupVcf {
	my ($GVaabbdd_samobj, $GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj)=@_;
	my $GVbb_express=0;
	(my $GVfixAlleles_index, $GVbb_express)=&AssignVariationAllele($GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj);
	my $GVtest_failed=0;###Running Control ###PAUSE
	my %GVvcf=%{$GVaabbdd_vcfobj}; undef $GVaabbdd_vcfobj; #ReadVcf object
#Format: %GVfixAlleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %GVfixAlleles_rnaseq=%{$GVfixAlleles_index};
#Format: %GVshare_readIDs = (chr => (readid1 => 1, readid2 => 1)), and then exclude those IDs assigned to alleles
	my %GVshare_readIDs=();
#Format: %GVreadid_at_pos=(chr => (pos => @readids)); All reads at each position
	my %GVreadid_at_pos=();
#Format: %GVreadid_by_allele=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
	my %GVreadid_by_allele=();
#Formar: %GVreadid_by_allele=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
#Format: %GVreadid2genome=(readID=> ('A' =>1, 'B' =>1, 'C'=1))
	my %GVreadid2genome=();
#Format: %GVreadid2chrmosome=(readID => ('1AL' => 1; 3B => 1, ...)
	my %GVreadid2chrmosome=();

	foreach my $GVchrom (@cluster_seqids) {

###COMMENT: retrieve all alignments mapped to this reference sequence
		my @GVchr_alignments=$GVaabbdd_samobj->get_features_by_location(-seq_id => "$GVchrom");
		foreach (@GVchr_alignments) {
			if ($_->qual >= $min_mapq) {#filter mapping quality
				${$GVshare_readIDs{$GVchrom}}{$_->name}++;
			}
		}
		my @GVpositions=sort {$a<=>$b} (keys %{$GVvcf{$GVchrom}});
		print "SUB(GroupVcf)Test: Reference2: $GVchrom\n" if ($debug or $verbose);
		print "SUB(GroupVcf)Test: Number of variations on $GVchrom: ", scalar(@GVpositions), "\n" if ($debug or $verbose);
###COMMENT: retrieve alignment at each position, uniform read ID assigned to each subgenome
		my %GVpos=();
		foreach my $GVind_pos (@GVpositions) {
			print "SUB(GroupVcf)Test:\t".$GVchrom, "\t", "Pos: $GVind_pos\tRef:${${$GVvcf{$GVchrom}}{$GVind_pos}}[0]\tVar:${${$GVvcf{$GVchrom}}{$GVind_pos}}[1]\tGen:${${$GVvcf{$GVchrom}}{$GVind_pos}}[2]", "\n" if ($debug or $verbose);
			@{$GVpos{$GVind_pos}}= $GVaabbdd_samobj->get_features_by_location(-seq_id => "$GVchrom", -start => "$GVind_pos", -end => "$GVind_pos");
			my $GVaa_allele='?'; my $GVbb_allele='?';my $GVdd_allele='?';
			$GVaa_allele=${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'A'} if (exists ${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'A'});
			$GVbb_allele=${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'B'} if (exists ${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'B'});
			$GVdd_allele=${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'C'} if (exists ${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'C'});
			my %GVallele_assign=();
#Format: %GVallele_assign=('A' => (allele1 =>1, allele2 =>1), 'B' => (allele1 =>1, allele2 =>1), 'D' => (allele1 =>1, allele2 =>1))
			foreach my $GVind_align (@{$GVpos{$GVind_pos}}) {
				if ($GVind_align->qual >= $min_mapq) {
					#$GVread_group: this read contains which allele genotype (0(ref)/1/2/3/.)
					my $GVread_group=&ReadVariantType($GVind_pos, ${${$GVvcf{$GVchrom}}{$GVind_pos}}[0], ${${$GVvcf{$GVchrom}}{$GVind_pos}}[1], ${${$GVvcf{$GVchrom}}{$GVind_pos}}[2], $GVind_align);
					my $GVgenome_chr=$GVind_align->get_tag_values("$bam_chromo_tag");
					if ($GVgenome_chr=~/\d\w{1,2}/) {
						while ($GVgenome_chr=~/(\d[ABDLS]{1,2})/g) {
							${${$GVreadid2chrmosome{$GVind_align->name}}}{$1}++;
						}
					}
					if ($GVread_group=~/^\d{1,2}$/) {
						push (@{${$GVreadid_at_pos{$GVchrom}{$GVind_pos}}}, $GVind_align->name);
						unless (exists $GVreadid2genome{$GVind_align->name}) {
							my $GVgenome_tagvalue=$GVind_align->get_tag_values("$bam_genome_tag");
							if ($GVgenome_tagvalue=~m/A/i) {
								${$GVreadid2genome{$GVind_align->name}}{'A'}=1;
							}
							if ($GVgenome_tagvalue=~m/B/i) {
								${$GVreadid2genome{$GVind_align->name}}{'B'}=1;
							}
							if ($GVgenome_tagvalue=~m/D/i) {
								${$GVreadid2genome{$GVind_align->name}}{'D'}=1;
							}
						}
						push (@{${${$GVreadid_by_allele{$GVchrom}}{$GVind_pos}}{$GVread_group}}, $GVind_align->name);
						if (exists ${$GVreadid2genome{$GVind_align->name}}{'A'}) {
							${$GVallele_assign{'A'}}{$GVread_group}++;
						}
						if (exists ${$GVreadid2genome{$GVind_align->name}}{'B'}) {
							${$GVallele_assign{'B'}}{$GVread_group}++;
						}
						if (exists ${$GVreadid2genome{$GVind_align->name}}{'D'}) {
							${$GVallele_assign{'D'}}{$GVread_group}++;
						}
					}
				}
				delete ${$GVshare_readIDs{$GVchrom}}{$GVind_align->name} if (exists ${$GVshare_readIDs{$GVchrom}}{$GVind_align->name});
			}
##COMMENT: retrieve best allele for each sungenome based on genome-info
			foreach my $GVsubgenome (keys %GVallele_assign) {
				my $GVmax=0;
				my $GVbest_allele='';
				my $GVbest_count=1;
				foreach my $GValleleat_pos (keys %{$GVallele_assign{$GVsubgenome}}) {
					if (${$GVallele_assign{$GVsubgenome}}{$GValleleat_pos} > $GVmax) {
						$GVbest_count=1;
						$GVbest_allele=$GValleleat_pos;
						$GVmax=${$GVallele_assign{$GVsubgenome}}{$GValleleat_pos};
						print "SUB(GroupVcf)Test: BestCount: $GVbest_count\n" if ($debug or $verbose);###for test###
					}
					elsif (${$GVallele_assign{$GVsubgenome}}{$GValleleat_pos} == $GVmax){
						$GVbest_count++;
					}
				}
				$GVbest_allele='?' if ($GVbest_count != 1);
				print "SUB(GroupVcf)Test: ".$GVsubgenome."\t".$GVbest_allele."\n" if ($debug or $verbose);
				if ($GVsubgenome eq 'A' and $GVaa_allele ne '?') {
					${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'A'}=$GVbest_allele;
				}
				elsif ($GVsubgenome eq 'B' and $GVbb_allele ne '?') {
					${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'B'}=$GVbest_allele;
					$GVbb_express=1;
				}
				elsif ($GVsubgenome eq 'D' and $GVdd_allele ne '?') {
					${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'D'}=$GVbest_allele;
				}
			}
		}
		print "SUB(GroupVcf)Test: number of share reads ids on $GVchrom: ".scalar(keys %{$GVshare_readIDs{$GVchrom}})."\n" if ($debug or $verbose);
	}
	return ($GVtest_failed, \%GVfixAlleles_rnaseq, $GVfixAlleles_index, \%GVreadid_by_allele, \%GVreadid2chrmosome, \%GVshare_readIDs, $GVbb_express);
}



### Read variant type based on the vcf genotypes
### &ReadVariantType(Pos, RefAllele, VarAllele, GenoType, Bio::DB::Sam alignment)
### Global: $debug
### Dependancy: Bio::DB::Sam, $geno_delimiter
### Note:
sub ReadVariantType {
	my ($RVTpos, $RVTref, $RVTvar, $RVTgeno, $RVTsam_align, )=@_;
	my $RVTref_seqid=$RVTsam_align->seq_id;
	my $RVTreadid=$RVTsam_align->name;
	my $RVTquery_seq=$RVTsam_align->query->dna;
	my $RVTreadstart=$RVTsam_align->start;
	my $RVTcigar=$RVTsam_align->cigar_str;
	my $RVTmdstring=$RVTsam_align->get_tag_values('MD');
	my $returnString='';
	print "\n\n\nSUB(ReadVariantType)Input: \n\tReference ID: ".$RVTref_seqid."\n\tVCF_pos: ".$RVTpos."\n\tGenotype: ".$RVTgeno."\n\tReadID: ".$RVTreadid."\n\tRead Seq: ".$RVTquery_seq."\n\tReadStart: ".$RVTreadstart."\n\tRefAllele: ".$RVTref."\n\tVarAllele: ".$RVTvar."\n\tCIGAR: ".$RVTcigar."\n\tMDstring: ".$RVTmdstring."\n" if ($debug);###test###
##COMMENT: Separate genotype, Assign number to each genotype: ref=0, var1=1, var2=2 ...
	my %RVTallele2geno=();
	$RVTallele2geno{$RVTref}=0;
#	print "SUB(ReadVariantType)Test: ref: $RVTref, $RVTallele2geno{$RVTref}\n";###test###
	my $RVTi=1;
	my @RVTvars=split(/$geno_delimiter/, $RVTvar);
#	print "VariantAllele: $RVTvar, SplitNumber: ".scalar(@RVTvars)."\n" if ($debug); ### For test ###
	my $RVTcat_is_deletion=0;
	foreach my $RVind_var (@RVTvars) {
		if ($RVind_var=~ /^[a-zA-Z]+$/) {
			$RVTallele2geno{$RVind_var}=$RVTi;
		}
		elsif ($RVind_var eq '.') {
			$RVTcat_is_deletion=1;
		}
		$RVTi++;
	}
#	map {print "Allele: $_, Geno: $RVTallele2geno{$_}\n" } keys %RVTallele2geno; ### For Test ###
###COMMENT: verify CIGAR length = query length
	my $RVTcigarOperations = &SplitCigar($RVTcigar);
	my $EVTcigar_cal_length=0;
	foreach (@$RVTcigarOperations){#calcular cigar length
		return 'L' unless (defined $_->[0] and defined $_->[1]);
		if ($_->[1] =~/^[MZIS=X]{1}$/) {
			$EVTcigar_cal_length+=$_->[0];
		}
	}
	unless ($EVTcigar_cal_length == length($RVTquery_seq)) {#if legnth(cigar) != length(read_seq), exit
		return 'L';
	}

#	foreach (keys %RVTallele2geno) {print "SUB(ReadVariantType)Test: ref2: $_, $RVTallele2geno{$_}\n";} ###test###
##COMMENT: Determine SNP/InDel vcf position for a read
	my $RVTrefPos = $RVTreadstart;
	my $RVTseqPos = 0;
	my $RVTthis_read_var='';
	my $EVTcapture_start_pos=0;
	my $RVTcapture_end_pos=0;
	my $RVTtest_var_start=0;
	my $RVTtest_var_end=0;
	my $RVTref_var_end=$RVTpos+length($RVTref);
	my $RVT5endloss=0;
	if ($RVTpos<$RVTreadstart) {
		$RVTpos=$RVTreadstart;
		$RVT5endloss=1;
	}
	foreach my $operation (@$RVTcigarOperations) {
		my $cig_length = $operation->[0];
		my $cig_op = $operation->[1];
		print "SUB(ReadVariantType)Test: Cigar: ($cig_length, $cig_op)\n" if ($debug);
		my ($RVTadd_ref, $RVTadd_seq)=(0, 0);
		my $RVTcigar_type=0;
		if ($cig_op =~ /^D$/) {
			$RVTadd_ref=$cig_length;
			$RVTcigar_type=3;
#			print "Deletion\n";###for test###
		}
		elsif ($cig_op =~ /^I$/) {
			$RVTadd_seq=$cig_length;
			$RVTcigar_type=2;
#			print "Insertion\n";###for test###
		}
		elsif ($cig_op =~ /^[M=]{1}$/) {
			$RVTadd_ref= $cig_length;
			$RVTadd_seq= $cig_length;
			$RVTcigar_type=1;
#			print "match\n";###for test###
		}
		elsif ($cig_op =~ /^S$/) {
			$RVTadd_seq=$cig_length;
			$RVTcigar_type=4;
		}
		elsif ($cig_op =~ /^N$/) {
			$RVTadd_ref= $cig_length;
			$RVTcigar_type=5;
		}
		elsif ($cig_op =~ /^P$/) {
			$RVTcigar_type=6;
		}
		elsif ($cig_op =~ /^H$/) {
			$RVTcigar_type=7;
		}
		else {
			print STDERR "SUB(ReadVariantType)Error: unrecognized Cigar State: $cig_op (Cigar: ". $RVTsam_align->cigar_str. " at position $RVTpos of Reference sequence ".$RVTsam_align->seq_id."\n";
			return '?';
		}
		$EVTcapture_start_pos=$RVTseqPos if ($RVTtest_var_start==0);
		$RVTcapture_end_pos=$RVTseqPos;
		if (($RVTrefPos <= $RVTpos) and ($RVTpos <($RVTrefPos+$RVTadd_ref))) {
			if (($RVTrefPos+$RVTadd_ref)<$RVTref_var_end) {
				if ($RVTtest_var_start==0) {
					if ($RVTtest_var_start==0 and $RVTcigar_type==1) {
						$EVTcapture_start_pos+=($RVTpos-$RVTrefPos);
						$RVTcapture_end_pos+=$RVTadd_seq;
						$RVTtest_var_start=1;
#						print "Starts1: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
					}
					else {
						return 'unknown2';
					}
				}
				else {
					if ($RVTcigar_type=~/^[23]{1}$/) {
						$RVTcapture_end_pos+=$RVTadd_seq;
#						print "Starts2: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
					}
					else {
						return 'unknown3';
					}
				}
			}
			else {
				if ($RVTtest_var_start==0 and $RVTcigar_type==1) {
					$EVTcapture_start_pos+=$RVTpos-$RVTrefPos;
					$RVTcapture_end_pos=$EVTcapture_start_pos+$RVTref_var_end-$RVTpos;
					$RVTtest_var_end=1;
#					print "Starts3: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
				}
				else {
					return 'Unknown1';
				}
			}
		}
		elsif ($RVTrefPos > $RVTpos) {
			if ($RVTcigar_type=~/^[4567]{1}$/) {
				return 'Nreads';
			}
			elsif (($RVTrefPos+$RVTadd_ref)<=$RVTref_var_end) {#insertion&deletion
				 $RVTcapture_end_pos+=$RVTadd_seq;
#				 print "Starts4: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
			}
			elsif ($RVTref_var_end>$RVTrefPos and $RVTref_var_end<($RVTrefPos+$RVTadd_ref)) {
					$RVTcapture_end_pos+=$RVTref_var_end-$RVTrefPos;
					$RVTtest_var_end=1;
#					print "Starts5: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
			}
		}
		last if ($RVTtest_var_end==1);
		$RVTrefPos+=$RVTadd_ref;
		$RVTseqPos+=$RVTadd_seq;
	}
	print "SUB(ReadVariantType)Test: Read pos: $RVTseqPos\n" if ($debug);
	print "SUB(ReadVariantType)Test: ReadLength: ".length($RVTquery_seq)."\tSubstr(Start - End): ".$EVTcapture_start_pos." - ".$RVTcapture_end_pos."\n" if ($debug);###For test###

	if ($RVTcapture_end_pos>length($RVTquery_seq)) {
		print STDERR "SUB(ReadVariantType)Test: Readlength: ".length($RVTquery_seq)."\tSubstrEnd: $RVTcapture_end_pos\n";
		print STDERR "SUB(ReadVariantType)Warnings: use read ($RVTreadid) length at $RVTpos (Ref:Var:Gonotype=$RVTref : $RVTvar : $RVTgeno) of $RVTref_seqid instead\n";###For test###
		$RVTcapture_end_pos=length($RVTquery_seq);
	}
	unless ($EVTcapture_start_pos < length($RVTquery_seq) and $EVTcapture_start_pos < $RVTcapture_end_pos) {
		print STDERR "SUB(ReadVariantType)Error: Unknown2\n";###For test###
		return "unknown";
	}
	my $RVTcapture_string='';
	if (($RVTcapture_end_pos-$EVTcapture_start_pos) >0) {
		$RVTcapture_string=substr($RVTquery_seq, $EVTcapture_start_pos, $RVTcapture_end_pos-$EVTcapture_start_pos);
	}
	elsif (($RVTcapture_end_pos-$EVTcapture_start_pos) ==0 and $RVT5endloss==1) {
		$RVTcapture_string='.';
	}
	if ($RVTcapture_string eq '') {
		print STDERR "SUB(ReadVariantType)Error: substr failed\n";###For test###
		return "unknown2";
	}
	print "SUB(ReadVariantType)Test: Captured String: $RVTcapture_string\n" if ($debug); ### For test ###
	my @returnarr=();
	if ($RVTtest_var_end) {
		if ($RVT5endloss==0) {
			if (exists $RVTallele2geno{$RVTcapture_string}) {
				return $RVTallele2geno{$RVTcapture_string};
			}
			else {
				return 'NotExist';
			}
		}
		else {
			foreach my $RVTind_allele (keys %RVTallele2geno) {
				if ($RVTind_allele=~/$RVTcapture_string$/) {
					push (@returnarr, $RVTallele2geno{$RVTind_allele});
				}
			}
		}
	}
	else {
		foreach my $RVTind_allele (keys %RVTallele2geno) {
			if ($RVT5endloss==0) {
				if ($RVTind_allele=~/^$RVTcapture_string/) {
					push (@returnarr, $RVTallele2geno{$RVTind_allele});
				}
			}
			else {
				if ($RVTind_allele=~/$RVTcapture_string/) {
					push (@returnarr, $RVTallele2geno{$RVTind_allele});
				}
			}
		}
	}
#	print "@returnarr\n";###test
	return (scalar(@returnarr)==1) ? $returnarr[0]:'M';	###Modify  
}



### FillInVariations
###&FillInVariations()
###Global: $debug
###Dependency:&ExtractAllele
###Note
sub FillInVariations {
	my ($FIVfixallele_hashindex, $FIVvcf_obj, $FIVgenoploidy)=@_;
##COMMENT: Test SUB(FillInVariations) successful or not 
##Format: 0=successful, 1=ExtraAlleles, 2=InsufficientAlleles, 3=Unknown
	my $FIVtest_cmd=0;
	my $FIVneed_runcompass=0;###Test need to run next step compass (if there is any unfixed allele) or not (if all alleles are fixed);
##Format: %GVfixAlleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %FIVfillin_alleles=%{$FIVfixallele_hashindex}; undef $FIVfixallele_hashindex;
##Format: $FIVreadvcf{chr}->{pos}=(ref, var, gen/gen2, ...);
	my %FIVreadvcf=%{$FIVvcf_obj}; undef $FIVfixallele_hashindex;
##Format: $FIVfixed_geno{chr}->{pos}='0/1/0'; ('A/B/D')
	my %FIVfixed_geno=();
##COMMENT: decide which subgenome to output
	my ($FIVout_aa, $FIVout_bb, $FIVout_dd)=(0, 0, 0);
	if ($FIVgenoploidy=~/A/) {
		$FIVout_aa=1;
#		print "AA\n"; ### For test ###
	}
	if ($FIVgenoploidy=~/B/) {
		$FIVout_bb=1;
#		print "BB\n"; ### For test ###
	}
	if ($FIVgenoploidy=~/D/) {
		$FIVout_dd=1;
#		print "DD\n"; ### For test ###
	}
	my $FIVploidy=$FIVout_aa+$FIVout_bb+$FIVout_dd;
	print "SUB(FillInVariations)Test: Total ploidy: $FIVploidy\n" if ($debug or $verbose);
##COMMENT: Fill in unknown alleles		
	FIVBLOCK1: {foreach my $FIVchrom (keys %FIVreadvcf) {
		my @FIVpositions=sort {$a<=>$b} (keys %{$FIVreadvcf{$FIVchrom}});
		print "SUB(FillInVariations)Test: Reference: $FIVchrom\n" if ($debug or $verbose);
		print "SUB(FillInVariations)Test: Number of variations on $FIVchrom: ", scalar(@FIVpositions), "\n"  if ($debug or $verbose);
		foreach my $FIVpos (@FIVpositions) {
			print "SUB(FillInVariations)Test: pos: $FIVpos\n" if ($debug or $verbose);
			print "SUB(FillInVariations)Test: @{${$FIVreadvcf{$FIVchrom}}{$FIVpos}}\n" if ($debug or $verbose);
			my ($FIVref, $FIVvar, $FIVgenos)=@{${$FIVreadvcf{$FIVchrom}}{$FIVpos}};
			my ($FIVallele_hashindex, $FIVgeno_arrayindex)=&ExtractAllele($FIVref, $FIVvar, $FIVgenos, 3);
			my %FIVallele_count=%{$FIVallele_hashindex}; undef $FIVfixallele_hashindex;
			my @FIVgeno=@{$FIVgeno_arrayindex};
			my $FIVnum_allele_notsure=0;
			my ($FIVaa_alleles,$FIVbb_alleles,$FIVdd_alleles)=('?', '?', '?');
			$FIVaa_alleles=${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'} if (exists ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'});
			$FIVbb_alleles=${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'} if (exists ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'});
			$FIVdd_alleles=${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'} if (exists ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'});
			my ($FIVmissing_aa, $FIVmissing_bb, $FIVmissing_dd)=(0, 0, 0);
			my %FIVallele_taken=();
			if ($FIVout_aa ==1) {
				if ($FIVaa_alleles ne '?' and exists $FIVallele_count{$FIVaa_alleles}) {
					$FIVallele_count{$FIVaa_alleles}--;
					$FIVallele_taken{$FIVaa_alleles}++;
				}
				else {
					$FIVnum_allele_notsure++;
					$FIVmissing_aa=1;
				}
			}
			if ($FIVout_bb ==1) {
				if ($FIVbb_alleles ne '?' and exists $FIVallele_count{$FIVbb_alleles}) {
					$FIVallele_count{$FIVbb_alleles}--;
					$FIVallele_taken{$FIVbb_alleles}++;
				}
				else {
					$FIVnum_allele_notsure++;
					$FIVmissing_bb=1;
				}
			}
			if ($FIVout_dd ==1) {
				if ($FIVdd_alleles ne '?' and exists $FIVallele_count{$FIVdd_alleles}) {
					$FIVallele_count{$FIVdd_alleles}--;
					$FIVallele_taken{$FIVdd_alleles}++;
				}
				else {
					$FIVnum_allele_notsure++;
					$FIVmissing_dd=1;
				}
			}
			my @FIVallele_left=();
			foreach (keys %FIVallele_count) {
				push (@FIVallele_left, $_) if ($FIVallele_count{$_}>0);
			}
			if (scalar(@FIVallele_left)<1 and $FIVnum_allele_notsure>0) {
				$FIVtest_cmd=2;
				print STDERR "SUB(FillInVariations)Error: no left allele at $FIVchrom:$FIVpos\n";###For test###
				last FIVBLOCK1;
			}
			my @FIVarr_geno_fix=();
			if ($FIVnum_allele_notsure==0) {
				push (@FIVarr_geno_fix, ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'}) if ($FIVout_aa==1);
				push (@FIVarr_geno_fix, ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'}) if ($FIVout_bb==1);
				push (@FIVarr_geno_fix, ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'}) if ($FIVout_dd==1);
			}
			elsif ($FIVnum_allele_notsure==1) {
				if (scalar(@FIVallele_left)!=1) {
					my @FIVallele_left_new=();
					foreach (@FIVallele_left) {
						push (@FIVallele_left_new, $_) unless (exists $FIVallele_taken{$_});
					}
					if (scalar(@FIVallele_left_new)!=1) {
						$FIVtest_cmd=1;
						print STDERR "SUB(FillInVariations)Error: extra at $FIVchrom: $FIVpos\n";###For test###
						print "AllelesLeft: @FIVallele_left_new\nAlleleNum: $FIVnum_allele_notsure\n";###For test###
						last FIVBLOCK1;
					}
					else {
						@FIVallele_left=@FIVallele_left_new;
					}
				}
				${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'}=shift @FIVallele_left if ($FIVmissing_aa ==1);
				${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'}=shift @FIVallele_left if ($FIVmissing_bb ==1);
				${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'}=shift @FIVallele_left if ($FIVmissing_dd ==1);
				push (@FIVarr_geno_fix, ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'}) if ($FIVout_aa==1);
				push (@FIVarr_geno_fix, ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'}) if ($FIVout_bb==1);
				push (@FIVarr_geno_fix, ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'}) if ($FIVout_dd==1);
			}
			elsif ($FIVnum_allele_notsure==2) {
				if (scalar(@FIVallele_left)==1){###Check if only one allele left
					$FIVallele_left[1]=$FIVallele_left[0];
				}
				if (scalar(@FIVallele_left)!=2) {###check if two allele left, return error if not
					$FIVtest_cmd=1;
					print STDERR "SUB(FillInVariations)Error: extra at $FIVchrom: $FIVpos\n";###For test###
					print "AllelesLeft: @FIVallele_left\nAlleleNum: $FIVnum_allele_notsure\n";###For test###
					last FIVBLOCK1;
				}
				my ($FIVallele_aa,$FIVallele_bb, $FIVallele_dd)=('?', '?', '?');
				if ($FIVout_aa==1) {
					if ($FIVmissing_aa ==1 and scalar(@FIVallele_left)>0) {
						$FIVallele_aa=shift @FIVallele_left;
						$FIVneed_runcompass=1;
					}
					elsif (exists ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'}) {
						$FIVallele_aa=${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'A'};
					}
					else {
						$FIVtest_cmd=3;
						print STDERR "SUB(FillInVariations)Error: AA allele at $FIVchrom: $FIVpos\n";###For test###
						last FIVBLOCK1;
					}
					push (@FIVarr_geno_fix, $FIVallele_aa);
				}
				if ($FIVout_bb==1) {
					if ($FIVmissing_bb ==1 and scalar(@FIVallele_left)>0) {
						$FIVallele_bb=shift @FIVallele_left;
						$FIVneed_runcompass=1;
					}
					elsif (exists ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'}) {
						$FIVallele_bb=${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'B'};
					}
					else {
						$FIVtest_cmd=3;
						print STDERR "SUB(FillInVariations)Error: BB allele at $FIVchrom: $FIVpos\n";###For test###
						last FIVBLOCK1;
					}
					push (@FIVarr_geno_fix, $FIVallele_bb);
				}
				if ($FIVout_dd==1) {
					if ($FIVmissing_dd ==1 and scalar(@FIVallele_left)>0) {
						$FIVallele_dd=shift @FIVallele_left;
						$FIVneed_runcompass=1;
					}
					elsif (exists ${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'}) {
						$FIVallele_dd=${${$FIVfillin_alleles{$FIVchrom}}{$FIVpos}}{'D'};
					}
					else {
						$FIVtest_cmd=3;
						print STDERR "SUB(FillInVariations)Error: DD allele at $FIVchrom: $FIVpos\n";###For test###
						last FIVBLOCK1;
					}
					push (@FIVarr_geno_fix, $FIVallele_dd);
				}
			}
			elsif ($FIVnum_allele_notsure==3) {
				@FIVarr_geno_fix=@FIVgeno;
				$FIVneed_runcompass=1;
			}
			${$FIVfixed_geno{$FIVchrom}}{$FIVpos}=join('/', @FIVarr_geno_fix);
			print "Fixed Alleles: ".${$FIVfixed_geno{$FIVchrom}}{$FIVpos}."\t$FIVchrom: $FIVpos\n" if ($debug);
		}
	}}###FIVBLOCK1;
	return ($FIVtest_cmd, $FIVneed_runcompass, \%FIVfixed_geno, \%FIVfillin_alleles);
}



### Read haptree phased alleles
### ReadHaptreeOut()
### Global: 
### Dependency: 
### Note: 
sub ReadHaptreeOut {
	my ($RHOphasefile, $RHOgeno_ploidy, $RHOfixalleles_hashindex)=@_;
	
	my $RHOsubinfo='SUB(ReadHaptreeOut)';
	my $RHOtest_cmd=0;
	unless (-s $RHOphasefile) {
		print STDERR "${RHOsubinfo}Error: invalid HapTree output: $RHOphasefile\n";
		return 1;
	}
### Initialising ploidy
	my $RHOnum_ploidy=length($RHOgeno_ploidy);
	my ($RHOaa_expressed, $RHObb_expressed, $RHOdd_expressed)=(0, 0, 0);
	my @RHOgenomes=();
	if ($RHOgeno_ploidy=~/A/) {
		$RHOaa_expressed=1;
		push (@RHOgenomes, 'A');
	}
	if ($RHOgeno_ploidy=~/B/) {
		$RHObb_expressed=1;
		push (@RHOgenomes, 'B');
	}
	if ($RHOgeno_ploidy=~/D/) {
		$RHOdd_expressed=1;
		push (@RHOgenomes, 'D');
	}
	if ($RHOnum_ploidy != $RHOaa_expressed + $RHObb_expressed + $RHOdd_expressed or ($RHOnum_ploidy <2 or $RHOnum_ploidy > 3)) {
		print STDERR "${RHOsubinfo}Error: ploidy (=$RHOnum_ploidy) error\n";
		return 1;
	}
##FORMAT: %RHOfixalleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %RHOfixalleles=%{$RHOfixalleles_hashindex}; undef $RHOfixalleles_hashindex;
	close PHASED if (defined fileno(PHASED));
	unless (open PHASED, "<", $RHOphasefile) {
		print STDERR "${RHOsubinfo}Error: can not $RHOphasefile\n";
		return 1;
	}
##FORMAT: file format
#BLOCK 	1000	1	1	1
#snp1	1000	0	0	1	
#snp2	1028	0	0	1	
#snp3	1118	1	1	0	
#snp4	1143	0	1	1	
#snp5	1320	0	0	1	
#snp6	1497	1	1	0	
#snp7	1518	1	0	0	
#snp8	1545	0	1	1	
#snp9	1677	1	1	0	
#snp10	1769	1	0	1	
#snp11	1795	0	1	1
	my $RHOblocknum=0;
	my %RHOtemp_block=();
	my %RHOfinal_fixallele=();
	while (my $RHOline=<PHASED>) {
		if ($RHOline=~/^BLOCK/) {
			$RHOblocknum++;
			next;
		}
		chomp $RHOline;
		my @RHOarr1=split(/\t/, $RHOline);
		next if (scalar(@RHOarr1) != ($RHOnum_ploidy+2));
		my ($RHOchrom, $RHOposit)=('', '');
		if ($RHOarr1[0]=~/^([\w+])_(\d+)$/) {
			$RHOchrom=$1;
			$RHOposit=$2;
		}
		else {
			print STDERR "${RHOsubinfo}Error: can not $RHOphasefile\n";
			return 1;
		}
		for (my $RHOi=2; $RHOi<scalar(@RHOarr1); $RHOi++) {
			${${${$RHOtemp_block{$RHOblocknum}}{$RHOchrom}}{$RHOposit}}{$RHOi-1}=$RHOarr1[$RHOi];
		}
	}
	close PHASED;
### count each group to each sungenome in each block
##FORMAT: %RHOtemp_assign=( block => (colnum => (A/B/D => total+count)));
	my %RHOtemp_assign=();
	my %RHOcorrelation=();
	foreach my $RHOblock (sort {$a <=> $b} (keys %RHOtemp_block)) {
		my %RHOfixarr=();
		my %RHOphase=();
		foreach my $RHOchrom2 (keys %{$RHOtemp_block{$RHOblock}}) {
			foreach my $RHOposit2 (keys %{${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}) {
				for (my $RHOi2=1; $RHOi2<=$RHOnum_ploidy; $RHOi2++) {
					push (@{$RHOphase{$RHOi2}}, ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
					if ($RHOaa_expressed) {
						if (exists ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'A'}) {
							${${$RHOtemp_assign{$RHOblock}}{$RHOi2}}{'A'}++ if (${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'A'} eq ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
							push (@{$RHOfixarr{'A'}}, ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'A'});
						}
						else {
							print STDERR "${RHOsubinfo}Error: Final AA assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
							return 1;
						}
					}
					if ($RHObb_expressed) {
						if (exists ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'B'}) {
							${${$RHOtemp_assign{$RHOblock}}{$RHOi2}}{'B'}++ if (${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'B'} eq ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
							push (@{$RHOfixarr{'B'}}, ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'B'});
						}
						else {
							print STDERR "${RHOsubinfo}Error: Final BB assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
							return 1;
						}
					}
					if ($RHOdd_expressed) {
						if (exists ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'D'}) {
							${${$RHOtemp_assign{$RHOblock}}{$RHOi2}}{'D'}++ if (${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'D'} eq ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
							push (@{$RHOfixarr{'D'}}, ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'D'});
						}
						else {
							print STDERR "${RHOsubinfo}Error: Final DD assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
							return 1;
						}
					}
				}
			}
		}
### Correlation
		foreach my $RHOind_fixallele (keys %RHOfixarr) {
			foreach my $RHOind_phasing (keys %RHOphase) {
				my ($RHOcorr_retcode, $RHOcov, $RHOcor, $RHOnum)=&CallCorrelation(\@{$RHOfixarr{$RHOind_fixallele}}, \@{$RHOphase{$RHOind_phasing}});
				${${$RHOcorrelation{$RHOblock}}{$RHOind_phasing}}{$RHOind_fixallele}=$RHOcor;
			}
		}
	}
### assign each group to each sungenome in each block
	my $RHOtestcmd_count=0;
	RHOBLOCK1: {foreach my $RHOind_block2 (sort {$a <=> $b} (keys %RHOtemp_assign)) {
		my $test_corelation=0;
		my %RHOgroup2genome=();
		my %RHOgenome2group=();
		RHOBLOCK2: {foreach my $RHOind_group (sort {$a <=> $b} (keys %{$RHOtemp_assign{$RHOind_block2}})) {
			my $RHObest_group='unknown';
			my $RHOmax_count=0;
			my $RHOnum_max=0;
			foreach my $RHOgenotype (@RHOgenomes) {
				if (exists ${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype}){
					if (${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype} > $RHOmax_count) {
						$RHOnum_max=1;
						$RHOmax_count=${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype};
						$RHObest_group=$RHOgenotype;
					}
					elsif (${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype} == $RHOmax_count) {
						$RHOnum_max++;
						$RHObest_group.="$RHOgenotype";
					}
				}
			}
			if ($RHOnum_max !=1 or $RHObest_group !~ m/^[ABD]{1}$/) {
				$test_corelation=1;
			}
			$RHOgroup2genome{$RHOind_group}=$RHObest_group;
			if (exists $RHOgenome2group{$RHObest_group}) {
				
				$test_corelation=1;
			}
			else {
				$RHOgenome2group{$RHObest_group}=$RHOind_group;
			}
		}}###RHOBLOCK2
		
	}}###RHOBLOCK1
}



### running mira
### &RunMira(fastq, manifest_name, output_fasta, seq_prefix, seq_desc)
### Global: $numthreads $freebayes_min_alternative_count
### Dependency: 
### Note: Be sure to chdir back in case of any error
sub RunMira {
	my ($RMfastq, $RMmira_manifest, $RMassembly_fasta, $RMseq_prefix, $RMseq_desc)=shift;
	my $RMtestcmd=0;
	my $RMreturn_fasta;
	if (! -s $RMfastq) {
		print STDERR "SUB(RunMira)Error: can not find fastq input\n";
		return 1;
	}
	my $RMcurpath=getcwd;
	my $RMrundir=&RetrieveDir($RMmira_manifest);
	unless (-d $RMrundir) {
		print STDERR "SUB(RunMira)Error: can not retrieve manifest path for $RMmira_manifest\n";
		return 1;
	}
	close MANIFEST if (defined fileno(MANIFEST));
	unless (open (MANIFEST, ">$RMmira_manifest")) {
		print STDERR "SUB(RunMira)Error: can not open manifest file: $RMmira_manifest\n";
		return 1;
	}
	my $RMproject='Ta';
	my $RMjob='denovo,est,accurate';
	my $RMparameters="COMMON_SETTINGS -GENERAL:number_of_threads=$numthreads -NAG_AND_WARN:check_template_problems=no:check_maxreadnamelength=no -CO:mark_repeats=yes:assume_snp_instead_repeat=yes:name_prefix=$RMseq_prefix -OUT:output_result_caf=no:output_result_tcs=no:output_result_maf=no SOLEXA_SETTINGS -CO:min_reads_per_group=$freebayes_min_alternative_count -AS:minimum_reads_per_contig=10";
	print MANIFEST "project = $RMproject\njob = $RMjob\nparameters = $RMparameters\n###Readgroup\n";
	print MANIFEST "readgroup = wheat\nautopairing\ndata = $RMfastq\ntechnology = solexa\ntemplate_size = 50 1000 autorefine\nsegment_placement = ---> <---\n";
	close MANIFEST;
	unless (chdir $RMrundir) {
		print STDERR "SUB(RunMira)Error: can not chdir to manifest path: $RMrundir\n";
		return 1;
	}
	if (! &exec_cmd_return("$path_mira4 $RMmira_manifest")) {
		print STDERR "SUB(RunMira)Error: MIRA return non-zero code\n";
		chdir $RMcurpath;
		return 1;
	}
	my @RMfasta_files=glob ("$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}*out*unpadded*fasta");#, "$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}_LargeContigs_out.unpadded.fasta");
	if (scalar(@RMfasta_files)<1) {
		print STDERR "SUB(RunMira)Error: MIRA output empty\n";
		return 1;
	}
	elsif (scalar(@RMfasta_files)==1) {
		$RMreturn_fasta=shift @RMfasta_files;
	}
	else {
		$RMreturn_fasta=$RMrundir."/merge.fa";
		close RMOUT unless (defined fileno(RMOUT));
		unless (open (RMOUT, ">$RMreturn_fasta")) {
			print STDERR "SUB(RunMira)Error: merge fasta1: open $RMreturn_fasta\n";
			chdir $RMcurpath;
			return 1;
		}
		foreach (@RMfasta_files) {
			unless (open (RMFA, "<$_")) {
				print STDERR "SUB(RunMira)Error: merge fasta2: open $_\n";
				chdir $RMcurpath;
				return 1;
			}
			while (<RMFA>) {
				print RMOUT $_;
			}
			close RMFA;
		}
		close RMOUT;
	}
	unless (rename ($RMreturn_fasta, $RMassembly_fasta) and -s $RMassembly_fasta) {
		print STDERR "SUB(RunMira)Error: move fasta error: $RMreturn_fasta to $RMassembly_fasta\n";
		return 1;
	}
	return 0;	
}



### run cd-hit-est to derep the fasta file MIRA4 assembled
### &CHE(input.fasta, output.fasta)
### Global: $path_cdhitest
### Dependency:
### Note:
sub CdHitEst {
	my ($CHEfastain, $CHEfastaout)=@_;
	unless (defined $CHEfastain and -s $CHEfastain) {
		print STDERR "SUB(CdHitEst)Error: fasta input for CDHIT not defined or exists\n";
		return 1;
	}
	unless (defined $CHEfastaout) {
		print STDERR "SUB(CdHitEst)Error: fasta output for CDHIT not defined\n";
		return 1;
	}
	if (-e $CHEfastaout) {
		print STDERR "SUB(CdHitEst)Warnings: fasta output for CDHIT existed but deleted\n";
		unlink $CHEfastaout;
	}
	if (! &exec_cmd_return("$path_cdhitest -i $CHEfastain -o $CHEfastaout -c 1.00 -n 10 -T 0 -r 1 -d 0 -M 30000")) {
		print STDERR "SUB(CdHitEst)Error: CDHIT running error: $CHEfastain\n";
		return 1;
	}
	elsif ( -s $CHEfastaout) {
		return 0;
	}
	else {
		print STDERR "SUB(CdHitEst)Error: CDHIT output error: $CHEfastain\n";
		return 1;
	}
}



###Detect a reference type
###&IsReference($var)
###Global: None
###Dependancy: Scalar::Util
sub IsReference {
	my $x=shift;
	use Scalar::Util 'reftype';
	my $reftype = reftype $x;
	if ( !defined $reftype ) {
    	print "\$x is not a reference.\n";
	}
	elsif ( $reftype eq 'HASH' ) {
    	print "hash\n";# do something with %$x
    	my @arr1=keys %$x;
    	print scalar(@arr1)."\n";
	}
	elsif ( $reftype eq 'ARRAY' ) {
		print "arr\n";# do something with @$x
	}
	elsif ( $reftype eq 'SCALAR' ) {
		print "scalar\n";# do something with $$x
	}
	else {
		print "unknown\n";# do something else
	}
}



###Test if value or values(delimited by /) contain 0
###&IsZeroIn($str)
###Global:none
###Dependancy: $geno_delimiter
sub IsZeroIn {
	my $IZIstr=shift;
	my @IZIarr=split(/\//, $IZIstr);
	my $ISIzero=1;
	foreach (@IZIarr) {
		$ISIzero=0 if ($_ ==0);
	}
	return $ISIzero;
}



### run bowtie2-build to index reference sequences
### RunBowtie2Index(reference.fasta, index_name)
### Global: $pathbowtie2build
### Dependency: &exec_cmd_return
### Note:
sub RunBowtie2Index {
	my ($RBIreference, $RBIindex)=@_;
	unless (defined $RBIreference and $RBIreference ne '' and -s $RBIreference) {
		print STDERR "SUB(RunBowtie2Index)Error: invalid fasta file for bowtie2-build\n";
		return 1;
	}
	unless (defined $RBIindex and $RBIindex ne '') {
		print STDERR "SUB(RunBowtie2Index)Error: invalid index name for bowtie2-build\n";
		return 1;
	}
##clean existing index
	unlink "$RBIindex.1.bt2" if (-e "$RBIindex.1.bt2");
	unlink "$RBIindex.2.bt2" if (-e "$RBIindex.2.bt2");
	unlink "$RBIindex.3.bt2" if (-e "$RBIindex.3.bt2");
	unlink "$RBIindex.4.bt2" if (-e "$RBIindex.4.bt2");
	unlink "$RBIindex.rev.1.bt2" if (-e "$RBIindex.rev.1.bt2");
	unlink "$RBIindex.rev.2.bt2" if (-e "$RBIindex.rev.2.bt2");
## run bowtie2-build
	if (! &exec_cmd_return("$pathbowtie2build -f $RBIreference $RBIindex")) {
		print STDERR "SUB(RunBowtie2Index)Error: bowtie2-build running error\n";
		return 1;
	}
	elsif (! -s "$RBIindex.1.bt2" or ! -s "$RBIindex.2.bt2" or ! -s "$RBIindex.3.bt2" or ! -s "$RBIindex.4.bt2" or ! -s "$RBIindex.rev.1.bt2" or ! -s "$RBIindex.rev.2.bt2") {
		print STDERR "SUB(RunBowtie2Index)Error: bowtie2-build output error\n";
		return 1;
	}
	return 0;
}


### Run bowtie2 and map read to reference, and index
### RunBowtie2p(reference, fq)
### Global: $pathbowtie2p, $fq_qual_score, $numthreads
### Dependency: 
### Note: $fq_qual_score='phred33'; $maxinsert=800
sub RunBowtie2p {
	my ($RBPindex, $RBPreadgroup, $RBPrep, $RBPfq_pair1, $RBPfq_pair2, $RBPbamout)=@_;
	unless (defined $RBPfq_pair1 and -s $RBPfq_pair1 and defined $RBPfq_pair2 and -s $RBPfq_pair2) {
		print STDERR "SUB(RunBowtie2p)Error: invalid bowtie2 R1 and R2\n";
		return 1;
	}
	(my $RBPbamout_base=$RBPbamout)=~s/\.\w+$//;
#Par${tissue}${RBPrep}_rep${RBPrep}.vs.$bt2index.sort
	unlink $RBPbamout if (-e $RBPbamout_base);
	unlink "$RBPbamout.bai" if (-e "$RBPbamout_base.bai");
	if (! &exec_cmd_return("$pathbowtie2p -q --$fq_qual_score --threads $numthreads --maxins $maxinsert --rg-id Par${RBPreadgroup}${RBPrep} --rg \"SM:Par${RBPreadgroup}\" --rg 'PL:ILLUMINA' --rg \"LB:Par${RBPreadgroup}\" -x $RunDir/4.mapping/$RBPindex -1 $RBPfq_pair1 -2 $RBPfq_pair2 | $path_samtools view -b -h -S - | $path_samtools sort - $RBPbamout_base" )) {
		print STDERR "SUB(RunBowtie2p)Error: bowtie2 running error\n";
		return 1;
	}
	elsif (! -s "Par${RBPreadgroup}${RBPrep}_rep${RBPrep}.vs.$RBPindex.sort.bam") {
		print STDERR "SUB(RunBowtie2p)Error: bowtie2 output error\n";
		return 1;
	}
	if (! &IndexBam("Par${RBPreadgroup}${RBPrep}_rep${RBPrep}.vs.$RBPindex.sort.bam")) {
		print STDERR "SUB(RunBowtie2p)Error: BAM index error\n";
		return 1;
	}
	return 0;
}
