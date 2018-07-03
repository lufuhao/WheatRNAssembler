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
use FuhaoPerl5Lib::BamKit qw/ReduceReadNameLength SamCleanHeader SortBam IndexBam CalcFPKM ReadSam/;
use FuhaoPerl5Lib::VcfKit qw/ReadVcf RunFreebayes ReadVariantType/;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::FileKit;
use FuhaoPerl5Lib::FastaKit qw/CdbFasta CdbYank IndexFasta RenameFasta RunFqTrinity/;
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

HapCompass
	--javapath	<[Opt] /path/to/java if not in PATH>
	--hapcompasspath	<[Opt] /path/to/hapcompass.jar>
	--hc2vcfpath	<[Opt] /path/to/hc2vcf.jar>

Trinity	
	--trinitypath <[Opt] /path/to/Trinity if not in PATH>
	--maxinsert <[Opt] Int: default 800>

Running LOG
	--logcluster	<[Opt] Cluster running LOG>
	--logfpkm	<[Opt] express LOG>
	--logallele <[Opt] Allele assignment LOG>

MISC
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

EOH



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
my ($file_totalreads);
#CDBtools
my ($path_cdbfasta, $path_cdbyank);
#express
#my ($path_express, $express_frag_len_mean, $express_frag_len_stddev, $express_max_read_len);###express
my ($minimum_fpkm, $maxinsert);
#freebayes
my ($path_freebayes, $freebayes_min_coverage, $freebayes_min_alternative_count, $vcfqual, $min_mapq, $freebayes_parallel, $path_freebayesparallel, $freebayes_fastabin);###freebayes
#samtools
my ($path_samtools, $path_tabix, $path_bgzip);###samtools
#HapCompass
my ($path_hapcompassjar, $path_hc2vcf, $path_java, $javamem);###HapCompass
#vcf-tools
my ($path_vcfmerge, $path_vcfconcat, $path_vcfsort);#VCFtools
#haptree
#my ($path_haptree);
#MIRA4
#my ($path_mira4);###MIRA4
my ($path_trinity, $num_cpus);
#CDHIT
#my ($path_cdhitest);
#Bowtie2
#my ($pathbowtie2build, $pathbowtie2p);
my $fq_qual_score;
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
	"numreads:s" => \$file_totalreads,
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
	"javapath:s" => \$path_java,
	"javamem:s" => \$javamem,
	"hapcompasspath:s" => \$path_hapcompassjar,
	"hc2vcfpath:s" => \$path_hc2vcf,
#	"haptreepath:s" => \$path_haptree,
	"vcfmergepath:s" => \$path_vcfmerge,
	"vcfconcatpath:s" => \$path_vcfconcat,
	"vcfsortpath:s" => \$path_vcfsort,
#	"mirapath:s" => \$path_mira4,
	"trinitypath:s" => \$path_trinity,
#	"cdhitestpath:s" => \$path_cdhitest,
#	"bowtie2buildpath:s" => \$pathbowtie2build,
#	"bowtie2path:s" => \$pathbowtie2p,
	"logcluster:s" => \$cluster_log,
	"logfpkm:s" => \$fpkm_log,
	"logallele:s" => \$allele_log,
	"phred:s" => \$fq_qual_score,
	"maxinsert:i" => \$maxinsert,
	"threads|t:i" => \$numthreads,
	"cpus:i" => \$num_cpus,
	"debug|" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
print "Setting up defaults\n"; ### For test ###
my $RootDir=$Bin;
my $cluster_start=0;
my $cluster_end;
($cluster_start, $cluster_end)=split(/-/, $list);
$debug=0 unless (defined $debug);
$numthreads=1 unless (defined $numthreads);
$num_cpus=1 unless (defined $num_cpus);
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
#$path_express='express' unless (defined $path_express);
#$express_frag_len_mean=250 unless (defined $express_frag_len_mean);
#$express_frag_len_stddev=300 unless (defined $express_frag_len_stddev);
#$express_max_read_len=150 unless (defined $express_max_read_len);
$minimum_fpkm=0 unless (defined $minimum_fpkm);
$maxinsert=800 unless (defined $maxinsert);
###Freebayes
$path_freebayes='freebayes' unless (defined $path_freebayes);
$freebayes_min_coverage=3 unless (defined $freebayes_min_coverage);
$freebayes_min_alternative_count=3 unless (defined $freebayes_min_alternative_count);
$vcfqual=20 unless (defined $vcfqual);
$min_mapq=0 unless (defined $min_mapq);
$freebayes_parallel=0 unless (defined $freebayes_parallel);
$path_freebayesparallel='freebayes-parallel' unless (defined $path_freebayesparallel);
$freebayes_fastabin=200 unless (defined $freebayes_fastabin);
#my $freebayes_additional_cmd=" --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --min-mapping-quality $min_mapq --genotype-qualities --pooled-discrete ";
my $freebayes_additional_cmd=" --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --min-mapping-quality $min_mapq --pooled-discrete ";
#samtools
$path_samtools='samtools' unless (defined $path_samtools);
$path_tabix='tabix' unless (defined $path_tabix);
$path_bgzip='bgzip' unless (defined $path_bgzip);
#HapCompass
$path_hapcompassjar="$RootDir/utils/hapcompass/hapcompass.jar" unless (defined $path_hapcompassjar);
$path_hc2vcf="$RootDir/utils/hapcompass/hc2vcf.jar" unless (defined $path_hc2vcf);
$path_java='java' unless (defined $path_java);
$javamem='2G' unless (defined $javamem and $javamem=~/^\d+G$/);
#VCFtools
$path_vcfmerge='vcf-merge' unless (defined $path_vcfmerge);
$path_vcfconcat='vcf-concat' unless (defined $path_vcfconcat);
$path_vcfsort='vcf-sort'  unless (defined $path_vcfsort);
#HapTree
#$path_haptree='haptree' unless (defined $path_haptree);
#MIRA4
#$path_mira4='mira' unless (defined $path_mira4);
$path_trinity='Trinity' unless (defined $path_trinity);
my $trinity_addcmd='--max_memory '.$javamem.' --run_as_paired --CPU '.$num_cpus.' --group_pairs_distance '.$maxinsert.' --full_cleanup --min_kmer_cov 3 --min_glue 3';
#CDHITEST
#$path_cdhitest='cd-hit-est' unless (defined $path_cdhitest);
#my $cdhit_addcmd=' -c 1.00 -n 10 -T 0 -r 1 -d 0 -M 30000 ';
#Bowtie2
#$pathbowtie2build='bowtie2-build' unless (defined $pathbowtie2build);
#$pathbowtie2p='bowtie2' unless (defined $pathbowtie2p);
#LOG
$cluster_log='0.cluster.log' unless (defined $cluster_log);
#print "### ClusterLog: $cluster_log\n";### For test ###
$fpkm_log='0.fpkm.log' unless (defined $fpkm_log);
$allele_log='0.allele.log' unless (defined $allele_log);
### Number of total reads
print "Reading Total number of reads\n";
my %reads_per_rg=();
close NUMREADS if (defined fileno(NUMREADS));
unless (open (NUMREADS, "<$file_totalreads")) {
	die "Error: can not open --numreads file\n";
}
while (my $read_num_line=<NUMREADS>) {
#	print "REadline: $read_num_line";### For test ###
	chomp $read_num_line;
	next if ($read_num_line=~/^#/);
	my @read_num_arr=split(/\s+/, $read_num_line);
	next unless (scalar(@read_num_arr)==2);
	next unless (defined $read_num_arr[1] and $read_num_arr[1] =~/^\d+$/);
#	print "Term: $read_num_arr[0]\t TotalReads: $read_num_arr[1]\n"; ### For test ###
	if ($read_num_arr[0] eq 'TotalAABBDD') {
		$totalreads_AABBDD=$read_num_arr[1];
		die "Error: wrong AABBDD read sum number\n" unless (defined $totalreads_AABBDD and $totalreads_AABBDD);
	}
	elsif ($read_num_arr[0] eq 'TotalAABB') {
		$totalreads_AABB=$read_num_arr[1];
		die "Error: wrong AABB read sum number\n" unless (defined $totalreads_AABB and $totalreads_AABB);
	}
	elsif ($read_num_arr[0] eq 'TotalAA') {
		$totalreads_AA=$read_num_arr[1];
		die "Error: wrong AA read sum number\n" unless (defined $totalreads_AA and $totalreads_AA);
	}
	elsif ($read_num_arr[0] eq 'TotalDD') {
		$totalreads_DD=$read_num_arr[1];
		die "Error: wrong DD read sum number\n" unless (defined $totalreads_DD and $totalreads_DD);
	}
	else {
		$reads_per_rg{$read_num_arr[0]}=$read_num_arr[1];
	}
}
close NUMREADS;

my @totalreads=($totalreads_AABBDD, $totalreads_AABB, $totalreads_AA, $totalreads_DD);

if (1) {### For test ###
	print "### SUM of reads number ###\n";
	print "AABBDD: $totalreads_AABBDD\n";
	print "AABB:   $totalreads_AABB\n";
	print "AA:     $totalreads_AA\n";
	print "DD:     $totalreads_DD\n";
	print "\n### SUM of reads group ###\n";
	map {print $_."\t".$reads_per_rg{$_}."\n"} keys %reads_per_rg;
	print "### SUM of reads number END ###\n";
}



### input and output ################################################
#Check important file inpput
print "Checking input and output\n";### For test ###
die "Error: undefined reference file\n" unless (defined $reference);
($reference)=AddFilePath($reference);
die "Error: not-found reference file\n" unless (-s $reference);
die "Please specify cluster file\n" unless (defined $file_cluster and -s $file_cluster);
die "Please specify the file of bam list for AABBDD\n" unless (defined $file_bam_aabbdd and -s $file_bam_aabbdd);
die "Please specify the file of bam list for AABB\n" unless (defined $file_bam_aabb and -s $file_bam_aabb);
die "Please specify the file of bam list for AA\n" unless (defined $file_bam_aa and -s $file_bam_aa);
die "Please specify the file of bam list for DD\n" unless (defined $file_bam_dd and -s $file_bam_dd);
die "Please specify the file for number of reads\n" unless (defined $file_totalreads and -s $file_totalreads);

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
	mkdir ("$RunDir/AABBDD", 0766) || die "(Main0))Error: can not create folder $RunDir/AABBDD\n";
}
unlink glob "$RunDir/AABBDD/*";
if (! -d "$RunDir/AABB") {
	mkdir ("$RunDir/AABB", 0766) || die "(Main0)Error: can not create folder AABB\n";
}
unlink glob "$RunDir/AABB/*";
if (! -d "$RunDir/AA") {
	mkdir ("$RunDir/AA", 0766) || die "(Main0)Error: can not create folder AA\n";
}
unlink glob "$RunDir/AA/*";
if (! -d "$RunDir/DD") {
	mkdir ("$RunDir/DD", 0766) || die "(Main0)Error: can not create folder DD\n";
}
unlink glob "$RunDir/DD/*";

#(Main1) index fasta
my ($test_cdbfasta, $fasta_index)=&CdbFasta($reference, $path_cdbfasta);
if ($test_cdbfasta == $cmd_failure_code ) {
	die "(Main1)Error: can not index fasta: $reference\n";
}
print "Fasta Index: $fasta_index\n";
# Load clustered fasta id list, 1 cluster perl line, saparated by space
close CLUSTER if (defined fileno(CLUSTER));
open (CLUSTER, "<$file_cluster") || die "(Main1)Error: can not open file: $file_cluster\n";
# Output cluster succeeds or not. Format: 
close CLUSTERLOG if (defined fileno(CLUSTERLOG));
open (CLUSTERLOG, ">$cluster_log") || die "(Main1)Error: can not write file: $cluster_log\n";
# output fpkm calculated
close FPKMLOG if (defined fileno(FPKMLOG));
open (FPKMLOG, ">$fpkm_log") || die "(Main1)Error: can not write file: $fpkm_log\n";
print FPKMLOG "###ClusterID\n###UnigeneID\tAABBDD\tAABB\tAA\tDD\n";
# output allele log
close ALLELELOG if (defined fileno(ALLELELOG));
open (ALLELELOG, ">$allele_log") || die "(Main1)Error: can not write file: $allele_log\n";
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
	next if ($cluster_line=~/^#/);
	print STDOUT "\n\n\n\n\n##### Prcessing Cluster $cluster_linenum ###\n$cluster_line\n";
	print STDERR "\n\n\n\n\n##### Prcessing Cluster $cluster_linenum ###\n";
	@cluster_seqids=(); ###Empty this in case of abnormal duplicates
	@cluster_seqids=split(/\s+/, $cluster_line);
##COMMENT: Check if empty line
	if (scalar(@cluster_seqids)<1) {
		print STDERR "(Main2)Warnings: line $cluster_linenum in $file_cluster ignored as empty\n";
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
	print "(Main3)Step: AABBDD BAMs\n"; ### For test ###
	if (! &GetBam(\@bam_AABBDD, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.nogroup.bam")) {
		print STDERR "(Main3)Error: AABBDD BAM extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGetBamAABBDD\n";
		next;
	}
	if (! &ReduceReadNameLength("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.RRNL.bam", $path_samtools)) {
		print STDERR "(Main12)Error: RRNL\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tRRNL\n";
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
			if (! &RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.nogroup.bam", 3, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", 0, $freebayes_additional_cmd, 0, 0, 1, 200, $path_freebayes, $path_freebayesparallel, $path_vcfsort, $path_vcfconcat, $path_tabix, $path_bgzip, $path_samtools)) {
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
	print "(Main4)Step: AABB BAMs\n"; ### For test ###
	if (! &GetBam(\@bam_AABB, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.nogroup.bam")) {
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
				if (! &RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.nogroup.bam", 2, "$RunDir/Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", $freebayes_additional_cmd, 0, 0, 1, 200, $path_freebayes, $path_freebayesparallel, $path_vcfsort, $path_vcfconcat, $path_tabix, $path_bgzip, $path_samtools)) {
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
	print "\n(Main5)Step: AA BAMs\n"; ### For test ###
	if (! &GetBam(\@bam_AA, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.nogroup.bam")) {
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
				if (! &RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.nogroup.bam", 2, "$RunDir/Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", $freebayes_additional_cmd, 0, 0, 1, 200, $path_freebayes, $path_freebayesparallel, $path_vcfsort, $path_vcfconcat, $path_tabix, $path_bgzip, $path_samtools)) {
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
	print "\n(Main6)Step: DD BAMs\n"; ### For test ###
	if (! &GetBam(\@bam_DD, \@cluster_seqids, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.group.bam", "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.nogroup.bam")) {
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
### &RunFreebayes($ref.fasta, \@seqids, $bam, $ploidy, $output.vcf.gz, $guided.vcf/0, $additions_cmd, [$run_parallel(1/0)], [Parallel_region file], [$threads], [$fastabin], [path_freebayes], [path_freebayesparallel], [path_vcfsort], [path_concat], [path_tabix], [path_bgzip], [path_samtools])
				if (! &RunFreebayes("$RunDir/Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.nogroup.bam", 2, "$RunDir/Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", $freebayes_additional_cmd, 1, 0, 1, 0, $path_freebayes, $path_freebayesparallel, $path_vcfsort, $path_vcfconcat, $path_tabix, $path_bgzip, $path_samtools)) {
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
	print "\n(Main7)Step: GroupVcf\n"; ### For test ###
	my %polymorphic_sites=();
	my $test_bb_expressed=0;
	my ($test_readbam, $aabbdd_bam_obj)=&ReadSam("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.RRNL.bam", "Clust$cluster_linenum/ref.$cluster_linenum.fa", 1);
	unless ($test_readbam) {
		print STDERR "(Main7)Error: Reading AABBDD BAM failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\treadsamAABBDD\n";
		next;
	}
	(my $test_groupvcf, my $rnaseq_allele_hashindex, my $fixed_allele_hashindex, my $fragments_hashindex, my $readidsum_hashindex, my $readid_alleles_hashindex, $test_bb_expressed)=&GroupVcf($aabbdd_bam_obj, $aabbdd_vcf_obj, $aabb_vcf_obj, $aa_vcf_obj, $dd_vcf_obj);
	if ($test_groupvcf) {
		print STDERR "(Main7)Error: GroupVcf SUB failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tGroupVcf\n";
		next;
	}
	my %FixedAllele_1rnaseq=%{$fixed_allele_hashindex}; undef $fixed_allele_hashindex;
#Format: %FixedAllele_2genome=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %FixedAllele_2genome=%{$rnaseq_allele_hashindex}; undef $rnaseq_allele_hashindex;
#Formar: %readsum=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
#								'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#								'shared' => 1/0,
#								$chrom => ($pos => (allele => quality))
#								RG => (1 => (seq1=>sequence, qual1=> quality), 2=> (seq1=>sequence, qual1=> quality), U => ...)
#							)
#					);
	my %readsum=%{$readidsum_hashindex}; undef $readidsum_hashindex;
#Format: %fragments=(readid => (chrom => (pos => (geno1 => quality, geno2 => quality)));
	my %fragments=%{$fragments_hashindex}; undef $fragments_hashindex;
#Format: %readids_by_alleles=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
	my %readids_by_alleles=%{$readid_alleles_hashindex}; undef $readid_alleles_hashindex;
	my $num_sites_to_phase=0;
	foreach my $chrom (keys %FixedAllele_2genome) {
		my @posits=keys %{$FixedAllele_2genome{$chrom}};
		foreach my $pos (@posits) {
			my %morphisms=();
			if (exists ${${$FixedAllele_2genome{$chrom}}{$pos}}{'A'}) {
				$morphisms{${${$FixedAllele_2genome{$chrom}}{$pos}}{'A'}}++;
#				print "KeyA: ".${${$FixedAllele_2genome{$chrom}}{$pos}}{'A'}."\n";### For test ###
			}
			if (exists ${${$FixedAllele_2genome{$chrom}}{$pos}}{'B'}) {
				$morphisms{${${$FixedAllele_2genome{$chrom}}{$pos}}{'B'}}++;
#				print "KeyB: ".${${$FixedAllele_2genome{$chrom}}{$pos}}{'B'}."\n";
			}
			if (exists ${${$FixedAllele_2genome{$chrom}}{$pos}}{'D'}) {
				$morphisms{${${$FixedAllele_2genome{$chrom}}{$pos}}{'D'}}++;
#				print "KeyD: ".${${$FixedAllele_2genome{$chrom}}{$pos}}{'D'}."\n";
			}
			my @morphism=keys %morphisms;
			if (scalar(@morphism)==1 and $morphism[0]=~/^\d+$/) {
				if (exists ${${$readids_by_alleles{$chrom}}{$pos}}{$morphism[0]}) {
					print "(Main7)Info: Chrom:Pos $chrom:$pos no phlymorphism: @morphism\n";### For test ###
					foreach (@{${${$readids_by_alleles{$chrom}}{$pos}}{$morphism[0]}}) {
						${$readsum{$_}}{'shared'}=1;
					}
#					print "Main(7)Info: Deleting FixedAllele_2genome keys: Chrom:Pos: $chrom:$pos\n";### For test ###
					delete ${$FixedAllele_2genome{$chrom}}{$pos};
#					if (exists ${$FixedAllele_2genome{$chrom}}{$pos}) {### For test ###
#						print "Main(7)Info: hash not deleted $chrom:$pos\n";
#					}
				}
			}
			else {
				print "(Main7)Info: Chrom:Pos $chrom:$pos has phlymorphism: @morphism\n";
				$polymorphic_sites{$chrom}{$pos}++;
				$num_sites_to_phase++;
			}
		}
	}
	if (1) {
		print "(Main7)Info: number of sites to be phased: $num_sites_to_phase\n";
#		foreach my $chrom (keys %FixedAllele_2genome) {
#			print "(main7)Test: Chrom: $chrom\n";
#			foreach my $pos (keys %{$FixedAllele_2genome{$chrom}}) {
#				print "(main7)Test: Chrom: $chrom Pos: $pos\n";
#			}
#		}
	}


##COMMENT: (Main8) decide ploidy here
	print "\n(Main8)Step: Ploidy\n"; ### For test ###
	my $test_run_phase=0;
	my $assembly_desc='';
	my $ploidy=1;
	my $geno_ploidy='';
	my %genomecount=();
	my %chromcount=();
	my $totalmappedreads=0;
	my $total_notag=0;
	my $total_tagged=0;
	my @final_ploidy=();
	foreach my $ind_read (keys %readsum) {
		if (exists ${$readsum{$ind_read}}{'gen'}) {
			foreach (keys %{${$readsum{$ind_read}}{'gen'}}) {
				$genomecount{$_}++;
			}
			$total_tagged++;
		}
		else {
			$total_notag++;
		}
		
		if (exists ${$readsum{$ind_read}}{'chr'}) {
			foreach (keys %{${$readsum{$ind_read}}{'chr'}}) {
				$chromcount{$_}++;
			}
		}
		if (exists ${$readsum{$ind_read}}{'mapped'}) {
			$totalmappedreads++;
		}
	}
	if (1) {### For test ###
		print "##### Genome summary #####\nTotal mapped reads: $totalmappedreads\nTagged: $total_tagged\nNotag: $total_notag\n\n";
		map {print "Genome: $_\tCount: $genomecount{$_}\n"} (sort (keys %genomecount)); 
		print "\n";
		map {print "Chrom: $_\tCount: $chromcount{$_}\n"} (sort (keys %chromcount));
		print "##### Genome summary #####\n";
		print "(Main8)Info: Ploidy: $ploidy\nAA: $test_aa_expressed\nBB: $test_bb_expressed\nDD: $test_dd_expressed\nPloid string: $geno_ploidy\n";
	}
	my $inferred_genoploidy='';
	my $inferred_chrom='';
	if (($total_tagged/$totalmappedreads) >=0.4) {
		foreach (keys %genomecount) {
			$inferred_genoploidy.=$_ if (($genomecount{$_}/$totalmappedreads) >=0.1);
			push (@final_ploidy, $_);
		}
		print "(Main8)Info: Inferred geno: $inferred_genoploidy\n";
		unless ($inferred_genoploidy =~/^[ABD]+$/) {
			print STDERR "(Main8)Error: unwanted \$inferred_genoploidy: $inferred_genoploidy\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tinferred_genoploidy\n";
			next;
		}
		my $num_chr_assigned=0;
		map {my $gen=substr($_, 1, 1); if ($inferred_genoploidy=~/$gen/) {$num_chr_assigned++;}} (keys %chromcount);
		print "(Main8)Info: Num_InferChrom: $num_chr_assigned\n";
		if ($num_chr_assigned>=1) {
			foreach (keys %chromcount) {
				my $gen=substr($_, 1, 1); 
				next unless ($inferred_genoploidy=~/$gen/);
				if ($chromcount{$_} >= ($totalmappedreads/$num_chr_assigned)) {
					$inferred_chrom.=$_;
					print "(Main8)Info: ChromInfer: $_\t$chromcount{$_}\n";
				}
			}
			print "(Main8)Info: Inferred geno: $inferred_chrom\n";
		}

		
	}
	else {
		if ($num_sites_to_phase==0) {
			$ploidy=1;
			$test_run_phase=0;
		}
		elsif ($num_sites_to_phase>0) {
			if ($num_sites_to_phase==1) {
				$test_run_phase=0;
			}
			else {
				$test_run_phase=1;
			}
			$ploidy=$test_aa_expressed+$test_dd_expressed+$test_bb_expressed;
			$geno_ploidy.='A' if ($test_aa_expressed>0);
			$geno_ploidy.='B' if ($test_bb_expressed>0);
			$geno_ploidy.='D' if ($test_dd_expressed>0);
			$geno_ploidy='ABD' if ($ploidy==0 or $geno_ploidy eq '');
		}
	}



##COMMENT: (Main9) Fill in genotypes accroding to the corrected alleles
#Formar: %polymorphismreadids=(readid => ('gen' => ('A' => ++, 'B' => ++, 'D' => ++), 
#										'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#										'shared' => 1/0,
#										$chrom => ($pos => (allele => quality))
#										)
#							);
	print "\n(Main9)Step: Fillin\n"; ### For test ###
	my %polymorphismreadids=();
	my %FixedAllele_3fillin=();
##Format: $fixedGeno_4fillin{chr}->{pos}='0/1/0'; ('A/B/D')
	my %fixedGeno_4fillin=();
	my %final_geno=();
	if ($ploidy==1) {
		$test_run_phase=0;
		if (0) {### For test ###
			print "(Main9)Info: ploidy==1; Skip SNP phasing\n";
		}
		%final_geno=%FixedAllele_2genome;
	}
	elsif ($ploidy==3 or $ploidy==2 or $ploidy==0) {
		(my $test_fiv, $test_run_phase, my $fixedgeno_hashindex, my $fillinAlleles_hashindex)=&FillInVariations(\%FixedAllele_2genome,$aabbdd_vcf_obj, $geno_ploidy, %readsum);
		if ($test_fiv==0) {###Fillin succeeds
			%FixedAllele_3fillin=%{$fillinAlleles_hashindex}; undef $fillinAlleles_hashindex;
			%fixedGeno_4fillin=%{$fixedgeno_hashindex}; undef $fixedgeno_hashindex;
		}
		elsif ($test_fiv==1) {
			print STDERR "(Main9)Error: FillIn1 SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin1\n";
			next;
		}
		elsif ($test_fiv==2) {
			print STDERR "(Main9)Error: FillIn2 SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin2\n";
			next;
		}
		elsif ($test_fiv==1) {
			print STDERR "(Main9)Error: FillIn3 SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin3\n";
			next;
		}
		else {
			print STDERR "(Main9)Error: FillIn SUB failed\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFillin4\n";
			next;
		}
	}
	else {
		print STDERR "(Main9)Error: unknown ploidy\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tPloidy\n";
		next;
	}
	
	if (0) {### Test %final_geno ###
		foreach my $chrom (keys %final_geno) {
			print "(Main9)Test: Chrom: $chrom\n";
			foreach my $pos (keys %{$FixedAllele_2genome{$chrom}}) {
				print "(Main9)Test: Chrom: $chrom Pos: $pos\n";
			}
		}
	}
##COMMENT: Main(10) SNP phasing
=oldhaptree
	if (! &RunHaptree("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.haptree.vcf", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.reads", \%fixedGeno_4fillin, $allele2readids_hashindex, $path_haptree)) {
		print STDERR "(Main10)Error: HapTree running\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tHaptree\n";
	}
=cut
	print "\n(Main10)Step: phasing\n"; ### For test ###
	if ($test_run_phase) {
		unless (&HapcompassVcf("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.hapcompass.vcf", \%fixedGeno_4fillin, \%polymorphic_sites)) {
			print STDERR "(Main10)Error: prepare HapCompass VCF error\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tHapcompassVcf\n";
			next;
		}
		unless (&HapcompassFragments(\%fragments, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.fragments")) {
			print STDERR "(Main10)Error: prepare HapCompass Fragments error\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFragments\n";
			next;
		}
		unless (&RunHapCompass(0, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.fragments", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.hapcompass.vcf", $ploidy, "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.fragments.hapcompass", "$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.phased.vcf", ' ', $path_java, $path_hapcompassjar, $path_hc2vcf)) {
			print STDERR "(Main10)Error: HapCompass running error\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tHapCompass\n";
			next;
		}
#		%final_geno=(); ## tobe defined
	}
	else {
		print "(Main10)Info: no need to run HapCompass\n";
	}
	foreach my $chrom (keys %{$aabbdd_vcf_obj}) {
		print "(Main10)Test: AlleleTesting: chromosome $chrom\n" if ($debug);
		foreach my $posit (keys %{${$aabbdd_vcf_obj}{$chrom}}) {
			print "(Main10)Test: AlleleTesting: positions $posit\n" if ($debug);
			my $printline=$cluster_linenum."\t".$chrom."\t".$posit."\t".${${${$aabbdd_vcf_obj}{$chrom}}{$posit}}[2]."\t";
			#%FixedAllele_1rnaseq
			if (exists $FixedAllele_1rnaseq{$chrom} and exists ${$FixedAllele_1rnaseq{$chrom}}{$posit}) {
				$printline.=(exists ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'A'}) ? ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'A'} : 'undef';
				$printline.="\t";
				$printline.=(exists ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'B'}) ? ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'B'} : 'undef';
				$printline.="\t";
				$printline.=(exists ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'D'}) ? ${${$FixedAllele_1rnaseq{$chrom}}{$posit}}{'D'} : 'undef';
				$printline.="\t";
			}
			else {
				$printline.='undef'."\t".'undef'."\t".'undef'."\t";
			}
			#%FixedAllele_2genome
			if (exists $FixedAllele_2genome{$chrom} and exists ${$FixedAllele_2genome{$chrom}}{$posit}) {
				$printline.=(exists ${${$FixedAllele_2genome{$chrom}}{$posit}}{'A'}) ? ${${$FixedAllele_2genome{$chrom}}{$posit}}{'A'} : 'xxxxx';
				$printline.="\t";
				$printline.=(exists ${${$FixedAllele_2genome{$chrom}}{$posit}}{'B'}) ? ${${$FixedAllele_2genome{$chrom}}{$posit}}{'B'} : 'xxxxx';
				$printline.="\t";
				$printline.=(exists ${${$FixedAllele_2genome{$chrom}}{$posit}}{'D'}) ? ${${$FixedAllele_2genome{$chrom}}{$posit}}{'D'} : 'xxxxx';
				$printline.="\t";
			}
			else {
				$printline.='undef'."\t".'undef'."\t".'undef'."\t";
			}
			#%FixedAllele_3fillin
			if (exists $FixedAllele_3fillin{$chrom} and exists ${$FixedAllele_3fillin{$chrom}}{$posit}) {
				$printline.=(exists ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'A'}) ? ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'A'} : 'undef';
				$printline.="\t";
				$printline.=(exists ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'B'}) ? ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'B'} : 'undef';
				$printline.="\t";
				$printline.=(exists ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'D'}) ? ${${$FixedAllele_3fillin{$chrom}}{$posit}}{'D'} : 'undef';
				$printline.="\t";
			}
			else {
				$printline.='undef'."\t".'undef'."\t".'undef'."\t";
			}
			#%fixedGeno_4fillin
			if (exists $fixedGeno_4fillin{$chrom}) {
				my $fillin_geno=(exists ${$fixedGeno_4fillin{$chrom}}{$posit}) ? ${$fixedGeno_4fillin{$chrom}}{$posit} : 'xxxxx';
				$printline.=$fillin_geno."\t";
			}
			else {
				$printline.='undef'."\t";
			}
			#$test_run_phase
			$printline.=$test_run_phase."\n";
			print ALLELELOG $printline;
		}
	}

	if (0) {### Test %final_geno ###
		foreach my $chrom (keys %final_geno) {
			print "(Main10)Test: Chrom: $chrom\n";
			foreach my $pos (keys %{$FixedAllele_2genome{$chrom}}) {
				print "(Main10)Test: Chrom: $chrom Pos: $pos\n";
				if (defined ${$FixedAllele_2genome{$chrom}}{$pos}){
					print "NOT defined\n";
				}
				else {
					print "Defined\n";
				}
			}
		}
	}



### (Main11) Get phlymorphic readids
#	%polymorphismreadids=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
	print "\n(Main11)Step: Get allelic reads\n"; ### For test ###
	if ($num_sites_to_phase>0 ) {
		foreach my $chrom (keys %final_geno) {
#Format: %FixedAllele_2genome=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))		
#		my @arr=keys %{$final_geno{$chrom}}; print "(Main11)Test: Chrom Pos: @arr\n";### For test ###
			foreach my $pos (keys %{$final_geno{$chrom}}) {
				unless (defined ${$final_geno{$chrom}}{$pos}) {### For test ###
					print "(Main11)Error: undefined final_geno at Chrom:Pos $chrom:$pos\n";
					next;
				}
				my @keys=keys %{${$final_geno{$chrom}}{$pos}};
				print "(Main11)Test: @keys\n";### For test ###
				if (scalar(@keys)<2 or scalar(@keys) >3) {
					print STDERR "(Main11)Error: InvalidGenoCall\n";
					print CLUSTERLOG $cluster_line."\tFail\t1\tInvalidGenoCall\n";
					next;
				}
				foreach my $ind_geno (@keys) {
					if ($ind_geno=~/^[ABD]{1}$/) {
						my $ind_allele=${${$final_geno{$chrom}}{$pos}}{$ind_geno};
						if (exists ${${$readids_by_alleles{$chrom}}{$pos}}{$ind_allele}) {
							foreach my $read_id (@{${${$readids_by_alleles{$chrom}}{$pos}}{$ind_allele}}) {
								${${$polymorphismreadids{$read_id}}{'gen'}}{$ind_geno}++;
							}
						}
						else {
							print STDERR "(Main11)Error: AllelesNonExist: Chrom:Pos:Alleles $chrom:$pos:$ind_allele\n";
							print CLUSTERLOG $cluster_line."\tFail\t1\tAllelesNonExist\n";
							next;
						}
					}
					else {
						print STDERR "(Main11)Error: NoABDGenoCall\n";
						print CLUSTERLOG $cluster_line."\tFail\t1\tNoABDGenoCall\n";
						next;
					}
				}
			}
		}
	}
	if (1) {
		print "(Main11)Test: Total mapped reads: ".$totalmappedreads."\n";
		print "(Main11)Test: Total allelic reads: ".scalar(keys %polymorphismreadids)."\n";
	}



### (Main12) BaM2fastq2
	print "\n(Main12)Step: bam2fastq\n"; ### For test ###
	my ($test_bam2fastq, $fastqfiles, $rgfpkms)=&Bam2FastQ2("$RunDir/Clust$cluster_linenum/AABBDD.$cluster_linenum.group.RRNL.bam", $inferred_genoploidy, \%readsum, \%polymorphismreadids, "$RunDir/Clust$cluster_linenum/reads", 1, 0, 1, $path_samtools);
	if (! $test_bam2fastq) {
		print STDERR "(Main12)Error: Bam2FastQ2\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tBam2FastQ2\n";
		next;
	}

### (Main13) Trinity
#use trinity instead of MIRA4 because MIRA4 does not support high coverage reads and poor performance during testing
	print "\n(Main13)Step: Trinity\n"; ### For test ###
	my %fasta_assembled=();
	foreach (keys %{$fastqfiles}) {
		unless (defined ${$fastqfiles}{$_} and -s ${$fastqfiles}{$_}) {
			print STDERR "(Main13)Error: Geno: $_ fastq file empty or not exists: ${$fastqfiles}{$_}\n";
			next;
		}
#		if (! &RunMira(${$fastqfiles}{$_}, "$RunDir/Clust$cluster_linenum/mira4.$_.manifest", "$RunDir/Clust$cluster_linenum/mira4.$_.fasta", "Clust${cluster_linenum}_", "Geno: $inferred_genoploidy Chrom: $inferred_chrom", 3, $path_mira4, 1)) {
#			print STDERR "(Main13)Error: MIRA-$_\n";
#			print CLUSTERLOG $cluster_line."\tFail\t1\tMIRA-$_\n";
#			next;
#		}
#		$fasta_assembled{$_}="$RunDir/Clust$cluster_linenum/mira4.$_.fasta");
		unless (RunFqTrinity(${$fastqfiles}{$_}, "$RunDir/Clust$cluster_linenum/Clust$cluster_linenum.trinity.$_.fasta", $trinity_addcmd, $path_trinity)) {
			print STDERR "(Main13)Error: Trinity-$_\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tTrinity-$_\n";
			next;
		}
		$fasta_assembled{$_}="$RunDir/Clust$cluster_linenum/Clust$cluster_linenum.trinity.$_.fasta";
	}
	unless (chdir $RunDir) {
		print STDERR "(Main13)Error: can not chdir to : $RunDir\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\tChDirAfterTrinity\n";
		next;
	}

=un-necessary
### (Main14) CD-HIT
	my %fasta_cdhit=();
	foreach (keys %fasta_assembled) {
		unless (&CdHitEst($fasta_assembled{$_}, "$fasta_assembled{$_}.cdhit.fasta", $cdhit_addcmd, $path_cdhitest)) {
			print STDERR "(Main14)Error: cdhit: $_\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tCdHitEst\n";
			next;
		}
		else {
			push (@fasta_cdhit, "$_.cdhit.fasta");
			unlink $_ if (-e $_);
		}
	}
=cut
### (Main14) rename fastq
	print "\n(Main14)Step: FastaRenamer\n"; ### For test ###
	foreach (keys %fasta_assembled) {
		unless (RenameFasta($fasta_assembled{$_}, "$fasta_assembled{$_}.rename.fa", "Cluster${cluster_linenum}_${_}_", 9, "Geno=$inferred_genoploidy Chrom=$inferred_chrom")) {
			print STDERR "(Main14)Error: FastaRenamer: $_\n";
			print CLUSTERLOG $cluster_line."\tFail\t1\tFastaRenamer\n";
			next;
		}
	}
### merge fasta
	
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
### Return: 1=success, 0=failed, 2=BAM no alignments
sub GetBam {
	my ($EBfiles_bam_index, $EBseq_obj, $EBoutput, $EBoutnogroup)=@_;
	
	my $EBsubinfo='SUB(GetBam)';
	my $EBsam_seq=join (' ', @{$EBseq_obj});
	my $EBi=0;
	my @EBtemp_bams=();
	
	return 0 if (scalar(@{$EBseq_obj})<1);###return error if no sequences
	return 1 if (-s $EBoutput and -s "$EBoutput.bai" and -s $EBoutnogroup and -s "$EBoutnogroup.bai");
	
##COMMENT: extract bam from each
#	print STDERR "${EBsubinfo}Info: sequence to be extracted: $EBsam_seq\n";###for test###
	foreach my $EBind_bam (@{$EBfiles_bam_index}) {
		$EBi++;
		my $EBbamout;
		if (scalar(@{$EBfiles_bam_index}) ==1) {
			$EBbamout="$EBoutput.multihead.bam";
		}
		else {
			$EBbamout="$EBoutput.$EBi.bam";
		}
		my $EBcmd="$path_samtools view -b -h -F 4 $EBind_bam $EBsam_seq > $EBbamout";
		if (! &exec_cmd_return($EBcmd)) {
			print STDERR "${EBsubinfo}Error: bam extract $EBsam_seq (Cluster: $cluster_linenum) from $EBind_bam failed\n";
			return 0;
		}
		else {
			push (@EBtemp_bams, "$EBoutput.$EBi.bam");
		}
	}
##COMMENT: merge bams 
	if (scalar(@{$EBfiles_bam_index}) > 1) {
		my $EBmerge_bams=join(' ', @EBtemp_bams);
		my $EBcmd2="$path_samtools merge $EBoutput.multihead.bam $EBmerge_bams";
		if (! &exec_cmd_return($EBcmd2)) {
			print STDERR "${EBsubinfo}Error: bam merge (Cluster: $cluster_linenum) from $EBmerge_bams failed\n";
			return 0;
		}
		else {
			if (-s "$EBoutput.multihead.bam") {
				map {unlink($_)} @EBtemp_bams;###delete temporary files ### For test ###
				
			}
			else{
				print STDERR "${EBsubinfo}Error: bam merge (Cluster: $cluster_linenum) from $EBmerge_bams not found\n";
				return 0;
			}
		}
	}
##COMMENT: cleanheader
	my $EBtest_cleanheader=&SamCleanHeader("$EBoutput.multihead.bam", "$EBoutput.withRG.bam", "$EBoutput.noRG.bam", $EBseq_obj,  $path_samtools);
	if ($EBtest_cleanheader==0) {
		print STDERR "${EBsubinfo}Error: bam header clean (Cluster: $cluster_linenum) failed\n";
		return 0;
	}
	elsif ($EBtest_cleanheader==2) {
		print STDERR "${EBsubinfo}Warnings: no alignments in BAM\n";
		return 2;### Note this ###
	}
	else {
		if ( -s "$EBoutput.noRG.bam" and -s "$EBoutput.withRG.bam") {
			unlink ("$EBoutput.multihead.bam") if (-s "$EBoutput.multihead.bam");
			if (! &SortBam("$EBoutput.withRG.bam", "$EBoutput.withRG.st.bam", $path_samtools, ' ')) {
				print STDERR "${EBsubinfo}Error: BAM sort running +RG error: $EBoutput.withRG.bam\n";
				return 0;
			}
			unlink "$EBoutput.withRG.bam" if (-s "$EBoutput.withRG.bam");
			if (! &MoveFile("$EBoutput.withRG.st.bam", $EBoutput)) {
				print STDERR "${EBsubinfo}Error: BAM moving: $EBoutput.withRG.bam\n";
				return 0;
			}
			if (! &IndexBam("$EBoutput", $path_samtools)) {
				print STDERR "${EBsubinfo}Error: BAM index running +RG error: $EBoutput\n";
				return 0;
			}
			if (! &SortBam("$EBoutput.noRG.bam", "$EBoutput.noRG.st.bam", $path_samtools)) {
				print STDERR "${EBsubinfo}Error: BAM sort running +RG error: $EBoutput.noRG.bam\n";
				return 0;
			}
			unlink "$EBoutput.noRG.bam" if (-s "$EBoutput.noRG.bam");
			if (! &MoveFile("$EBoutput.noRG.st.bam", $EBoutnogroup)) {
				print STDERR "${EBsubinfo}Error: BAM moving: $EBoutput.withRG.bam\n";
				return 0;
			}
			if (! &IndexBam($EBoutnogroup, $path_samtools)) {
				print STDERR "${EBsubinfo}Error: BAM index running +RG error: $EBoutnogroup\n";
				return 0;
			}
		}
		else{
			print STDERR "${EBsubinfo}Error: bam headerclean output (Cluster: $cluster_linenum) not found\n";
			return 0;
		}
	}
	return 1;
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


###Asigning SNP allele 
###&AssignVariationAllele()
###Global: $debug, $verbose
###Dependancy: &ExtractAllele
sub AssignVariationAllele {
	my ($AVAaabbdd_vcfindex, $AVAaabb_vcf_index, $AVAaa_vcf_index, $AVAdd_vcf_index)=@_;

	my $AVAtest_AA_expressed=1;
	$AVAtest_AA_expressed=0 if (scalar(keys %{$AVAaa_vcf_index})==0);
	my $AVAtest_DD_expressed=1;
	$AVAtest_DD_expressed=0 if (scalar(keys %{$AVAdd_vcf_index})==0);
	my $AVAtest_bb_expressed=0;
#Format: %AVAgenome2allele=(chr => (pos => ('A' => allele, 'B' => Allele, 'D' => allele)))
	my %AVAgenome2allele=();
#Format: %AVAallele2genome=(chr => (pos => (allele => 'A', allele => 'B', allele => 'D')))
	my %AVAallele2genome=();
	foreach my $AVAchrom (keys %{$AVAaabbdd_vcfindex}) {
		foreach my $AVApos (sort {$a<=>$b} keys %{${$AVAaabbdd_vcfindex}{$AVAchrom}}) {
			print "SUB(AssignVariationAllele)Test: Chr: $AVAchrom\tPos: $AVApos\t Arr: @{${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}\n" if ($debug or $verbose);
			my ($AVAaabbdd_genohash_index, $AVAaabbdd_genoarr_index)=&ExtractAllele(${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[2], 1);
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
				if (exists ${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}) {
					my ($AVAdd_genohash_index, $AVAdd_genoarr_index)=&ExtractAllele(${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[2], 1);
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
				if (exists ${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}) {
					my ($AVAaa_genohash_index, $AVAaa_genoarr_index)=&ExtractAllele(${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[2], 1);
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
				if ($AVAtest_AA_expressed==1 and $AVAaa_allele=~/^\d$/ and exists ${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}) {
					my ($AVAaabb_genohash_index, $AVAaabb_genoarr_index)=&ExtractAllele(${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[2], 1);
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
					print "chr: $AVAchrom\tPos: $AVApos\tRef: ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[0]\tVar: ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[1]\tAABBDDgeno ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[2]\tAABBgeno: ";
					if (exists ${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}) {
						print ${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tAA: ";
					if (exists ${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}) {
						print ${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tDD: ";
					if (exists ${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}) {
						print ${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[2];
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






###group SNPs based on read length
###&GroupVcf($ReadSam_obj, AABBDD_Vcf_obj, AABB_vcf_obj, AA_vcf_obj, DD_vcf_obj)
###Global: $debug, $verbose
###Dependancy:&ReadVcf, &ReadVariantType, $min_mapq, $bam_genome_tag
sub GroupVcf {
	my ($GVaabbdd_samobj, $GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj)=@_;
	my $GVbb_express=0;
	(my $GVfixAlleles_index, $GVbb_express)=&AssignVariationAllele($GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj);

#Format: %GVfixAlleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %GVfixAlleles_rnaseq=%{$GVfixAlleles_index};
#Format: %GVshare_readIDs = (chr => (readid1 => 1, readid2 => 1)), and then exclude those IDs assigned to alleles
#	my %GVshare_readIDs=();
#Format: %GVreadid_at_pos=(chr => (pos => @readids)); All reads at each position
#	my %GVreadid_at_pos=();
#Format: %GVreadid_by_allele=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
	my %GVreadid_by_allele=();
#Format: %GVreadfragments=(readid => (chrom => (pos => (geno1 => quality, geno2 => quality)));
	my %GVreadfragments=();
#Formar: %GVreadidpairs_by_allele=(chr => (pos => (allele1 => (readID => (1/2 => qual), allele2 => (readID => (1/2 => qual))))))
	my %GVreadidpairs_by_allele=();
#Format: %GVreadid2genome=(readID=> ('A' =>1, 'B' =>1, 'C'=1))
	my %GVreadid2genome=();
#Format: %GVreadid2chrmosome=(readID => ('1AL' => 1; 3B => 1, ...)
#	my %GVreadid2chrmosome=();
#Formar: %GVreadidsum=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
#								'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#								'shared' => 1/0,
#								chr => (pos => (allele => quality))
#							)
#					);
	my %GVreadidsum=();

	foreach my $GVchrom (@cluster_seqids) {

###COMMENT: retrieve all alignments mapped to this reference sequence
		my @GVchr_alignments=$GVaabbdd_samobj->get_features_by_location(-seq_id => "$GVchrom");
		foreach my $GVind_align (@GVchr_alignments) {
			next if ($GVind_align->flag & 0x0004);# makesure mapped
			next unless ($GVind_align->qual >= $min_mapq);#filter mapping quality
#			${$GVshare_readIDs{$GVchrom}}{$GVind_align->name}++;
			${$GVreadidsum{$GVind_align->name}}{'mapped'}++;
			${$GVreadidsum{$GVind_align->name}}{'shared'}++;
			my $GVtag_chrom=$GVind_align->get_tag_values("$bam_chromo_tag");
			while ($GVtag_chrom=~/(\d[ABDLS]{1,2})/g) {
#				${${$GVreadid2chrmosome{$GVind_align->name}}}{$1}++;
				${${$GVreadidsum{$GVind_align->name}}{'chr'}}{$1}++;
			}
			
			my $GVgenome_tagvalue=$GVind_align->get_tag_values("$bam_genome_tag");
			if ($GVgenome_tagvalue=~m/A/i) {
				${$GVreadid2genome{$GVind_align->name}}{'A'}=1;
				${${$GVreadidsum{$GVind_align->name}}{'gen'}}{'A'}++;
			}
			if ($GVgenome_tagvalue=~m/B/i) {
				${$GVreadid2genome{$GVind_align->name}}{'B'}=1;
				${${$GVreadidsum{$GVind_align->name}}{'gen'}}{'B'}++;
			}
			if ($GVgenome_tagvalue=~m/D/i) {
				${$GVreadid2genome{$GVind_align->name}}{'D'}=1;
				${${$GVreadidsum{$GVind_align->name}}{'gen'}}{'D'}++;
			}
		}
		my @GVpositions=sort {$a<=>$b} (keys %{${$GVaabbdd_vcfobj}{$GVchrom}});
		print "SUB(GroupVcf)Test: Reference2: $GVchrom\n" if ($debug or $verbose);
		print "SUB(GroupVcf)Test: Number of variations on $GVchrom: ", scalar(@GVpositions), "\n" if ($debug or $verbose);
###COMMENT: retrieve alignment at each position, uniform read ID assigned to each subgenome
		my %GVpos=();
		foreach my $GVind_pos (@GVpositions) {
			print "SUB(GroupVcf)Test:\t".$GVchrom, "\t", "Pos: $GVind_pos\tRef:${${${$GVaabbdd_vcfobj}{$GVchrom}}{$GVind_pos}}[0]\tVar:${${${$GVaabbdd_vcfobj}{$GVchrom}}{$GVind_pos}}[1]\tGen:${${${$GVaabbdd_vcfobj}{$GVchrom}}{$GVind_pos}}[2]", "\n" if ($debug or $verbose);
			@{$GVpos{$GVind_pos}}= $GVaabbdd_samobj->get_features_by_location(-seq_id => "$GVchrom", -start => "$GVind_pos", -end => "$GVind_pos");
			my $GVaa_allele='?'; my $GVbb_allele='?';my $GVdd_allele='?';
			$GVaa_allele=${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'A'} if (exists ${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'A'});
			$GVbb_allele=${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'B'} if (exists ${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'B'});
			$GVdd_allele=${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'C'} if (exists ${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'C'});
			my %GVallele_assign=();
#Format: %GVallele_assign=('A' => (allele1 =>1, allele2 =>1), 'B' => (allele1 =>1, allele2 =>1), 'D' => (allele1 =>1, allele2 =>1))
			foreach my $GVind_align (@{$GVpos{$GVind_pos}}) {
				next if ($GVind_align->flag & 0x0004);
				if ($GVind_align->qual >= $min_mapq) {
					#$GVread_group: this read contains which allele genotype (0(ref)/1/2/3/.)
					my ($GVread_group, $GVcapture_qual)=&ReadVariantType($GVind_pos, ${${${$GVaabbdd_vcfobj}{$GVchrom}}{$GVind_pos}}[0], ${${${$GVaabbdd_vcfobj}{$GVchrom}}{$GVind_pos}}[1], ${${${$GVaabbdd_vcfobj}{$GVchrom}}{$GVind_pos}}[2], $GVind_align, 2);
					
					

					if ($GVread_group=~/^\d{1,2}$/) {
#						push (@{${$GVreadid_at_pos{$GVchrom}{$GVind_pos}}}, $GVind_align->name);
###For HapCompass fragments  file
						if (! defined $GVcapture_qual or length($GVcapture_qual) !=1) {
							$GVcapture_qual='I'; ### The second highest phred+33 score
						}
						elsif (ord($GVcapture_qual)<33 or ord ($GVcapture_qual) >74) {
							$GVcapture_qual='I'; ### The second highest phred+33 score
						}
						if (exists ${${${$GVreadfragments{$GVind_align->name}}{$GVchrom}}{$GVind_pos}}{$GVread_group}) {
							if (ord($GVcapture_qual) > ord(${${${$GVreadfragments{$GVind_align->name}}{$GVchrom}}{$GVind_pos}}{$GVread_group})) {
								${${${$GVreadfragments{$GVind_align->name}}{$GVchrom}}{$GVind_pos}}{$GVread_group}=$GVcapture_qual;
							}
						}
						else {
							${${${$GVreadfragments{$GVind_align->name}}{$GVchrom}}{$GVind_pos}}{$GVread_group}=$GVcapture_qual;
						}
						push (@{${${$GVreadid_by_allele{$GVchrom}}{$GVind_pos}}{$GVread_group}}, $GVind_align->name);
=old
						if ($GVind_align->flag & 0x0001) {
							if ($GVind_align->flag & 0x0040) {
								${${${${$GVreadidpairs_by_allele{$GVchrom}}{$GVind_pos}}{$GVread_group}}{$GVind_align->name}}{1}=$GVcapture_qual;
							}
							elsif ($GVind_align->flag & 0x0080) {
								${${${${$GVreadidpairs_by_allele{$GVchrom}}{$GVind_pos}}{$GVread_group}}{$GVind_align->name}}{2}=$GVcapture_qual;
							}
							else {
								print STDERR "SUB(GroupVcf)Error: Unknown FLAG pair 1/2 (FLAG: ".$GVind_align->flag.")\n" if ($debug or $verbose);
								return 1;
							}
						}
						else {
							${${${${$GVreadidpairs_by_allele{$GVchrom}}{$GVind_pos}}{$GVread_group}}{$GVind_align->name}}{1}=$GVcapture_qual;
						}
=cut
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
#				delete ${$GVshare_readIDs{$GVchrom}}{$GVind_align->name} if (exists ${$GVshare_readIDs{$GVchrom}}{$GVind_align->name});
				${$GVreadidsum{$GVind_align->name}}{'shared'}=0;
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
#		print "SUB(GroupVcf)Test: number of share reads ids on $GVchrom: ".scalar(keys %{$GVshare_readIDs{$GVchrom}})."\n" if ($debug or $verbose);
	}
#	return (0, \%GVfixAlleles_rnaseq, $GVfixAlleles_index, \%GVreadid_by_allele, \%GVreadid2chrmosome, \%GVshare_readIDs, $GVbb_express); ### For HapTree ###
	return (0, \%GVfixAlleles_rnaseq, $GVfixAlleles_index, \%GVreadfragments, \%GVreadidsum, \%GVreadid_by_allele, $GVbb_express); ### For HapCompass 20150621 ###
	
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
			$FIVaa_alleles=${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'} if (exists ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'});
			$FIVbb_alleles=${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'} if (exists ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'});
			$FIVdd_alleles=${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'} if (exists ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'});
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
				push (@FIVarr_geno_fix, ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'}) if ($FIVout_aa==1);
				push (@FIVarr_geno_fix, ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'}) if ($FIVout_bb==1);
				push (@FIVarr_geno_fix, ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'}) if ($FIVout_dd==1);
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
				${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'}=shift @FIVallele_left if ($FIVmissing_aa ==1);
				${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'}=shift @FIVallele_left if ($FIVmissing_bb ==1);
				${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'}=shift @FIVallele_left if ($FIVmissing_dd ==1);
				push (@FIVarr_geno_fix, ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'}) if ($FIVout_aa==1);
				push (@FIVarr_geno_fix, ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'}) if ($FIVout_bb==1);
				push (@FIVarr_geno_fix, ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'}) if ($FIVout_dd==1);
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
					elsif (exists ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'}) {
						$FIVallele_aa=${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'A'};
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
					elsif (exists ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'}) {
						$FIVallele_bb=${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'B'};
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
					elsif (exists ${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'}) {
						$FIVallele_dd=${${${$FIVfixallele_hashindex}{$FIVchrom}}{$FIVpos}}{'D'};
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
	return ($FIVtest_cmd, $FIVneed_runcompass, \%FIVfixed_geno, \%{$FIVfixallele_hashindex});
}




### convert bam files into fastq
### Bam2FastQ2 ($bamin, 'ABD', %readsum, %allelicreads, outfqprefix_path, [map_code], [MAPQ_code], [RGaware], [path_samtools])
### Global:
### Dependency: $totalreads_AABBDD $reads_per_rg{$_}
### Note: map_code (0=all, 1=mapped, 2=unmapped)
### Note: [RGaware] 0=false; 1=true
sub Bam2FastQ2 {
	my ($BFQbamin, $BFQgeno, $BFQreadsum, $BFQallelicreads, $BFQoutfqprefix, $BFQmapcode, $BFQmapq, $BFQrg_aware, $BFQpath_samtools)=@_;
#Formar: %$BFQreadsum=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
#								'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#								'shared' => 1/0,
#								$chrom => ($pos => (allele => quality))
#								'reads' => (1 => (seq1=>sequence, qual1=> quality), 2=> (seq1=>sequence, qual1=> quality), U => ...)
#							)
#					);
	my $BFQsubinfo='SUB(Bam2FastQ2)';
	$BFQpath_samtools='samtools' unless (defined $BFQpath_samtools);
	$BFQmapcode=0 unless (defined $BFQmapcode);
	$BFQmapq=0 unless (defined $BFQmapq);
	$BFQrg_aware=0 unless (defined $BFQrg_aware);
	my %excludedreads=();
	my $BamKit_failure=0;
	my $BamKit_success=1;
##Format: %returnhash=(A => "A.fastq", B => b.fastq, D => D.fastq);
	my %returnhash=();
##Format: %BFQfpkms=(chr => ('chr1' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#											RG2 => ...)
	my %BFQfpkms=();
##Format: %BFQreadpair=(chr => (RG1 => (readid => 1 => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
#													2 => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
#													'U' => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
#								RG2...);
	my %BFQreadpair=();
	my %BFQgenos=();
	$BFQoutfqprefix='./' unless (defined $BFQoutfqprefix);
	my %BFQreflen=();
#	print "${BFQsubinfo}Test: \n\tBAM: $BFQbamin\n\tGeno: $BFQgeno\n\tFqPrefix: $BFQoutfqprefix\n\tMapcode: $BFQmapcode\n\tMapQ: $BFQmapq\n\tReadGroupAware: $BFQrg_aware\n\tSAMtools: $BFQpath_samtools\n";### For test ###
	
	unless (defined $BFQbamin and -s $BFQbamin) {
		print STDERR "${BFQsubinfo}Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $BFQgeno and $BFQgeno =~/^[ABD]+$/) {
		print STDERR "${BFQsubinfo}Error: unknown geno code\n";
		return $BamKit_failure;
	}
	$BFQgenos{'A'}++ if (defined $BFQgeno and $BFQgeno=~/A/);
	$BFQgenos{'B'}++ if (defined $BFQgeno and $BFQgeno=~/B/);
	$BFQgenos{'D'}++ if (defined $BFQgeno and $BFQgeno=~/D/);
#	print "${BFQsubinfo}Test: \t"; map {print "Geno: $_\t"} (keys %BFQgenos); print "\n";### For test ###
	close BAMIN if (defined fileno(BAMIN));
	unless (open(BAMIN, "$BFQpath_samtools view -h $BFQbamin | ")) {
		print STDERR "${BFQsubinfo}Error: open BAM: $BFQbamin \n";
		return $BamKit_failure;
	}
	close BAMOUT if (defined fileno(BAMOUT));
	unless (open (BAMOUT, " | samtools view -bhS - > $BFQbamin.excluded.bam")) {
		print STDERR "${BFQsubinfo}Error: write BAM: $BFQbamin.excluded.bam \n";
		return $BamKit_failure;
	}
	my $BFQnumline=0;
	while (my $BFQline1=<BAMIN>) {
		$BFQnumline++;
		chomp $BFQline1;
#Get reference seq length
		if ($BFQline1=~/^\@/) {
			if ($BFQline1=~/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/) {
				$BFQreflen{$1}=$2;
			}
			print BAMOUT $BFQline1."\n";
			next;
		}
		my @BFQarr=split(/\t/, $BFQline1);
#Check column number
		if (scalar(@BFQarr)<11) {
			print STDERR "${BFQsubinfo}Warnings: col<11 at line $BFQnumline (Readid: $BFQarr[0]) in BAM $BFQbamin\n";
			next;
		}
#check if exists in hash
		my $BFQis_shared=0;
#		print "${BFQsubinfo}Test: readname: $BFQarr[0]\n"; ### For test ###
		if (exists ${$BFQreadsum}{$BFQarr[0]} and exists ${${$BFQreadsum}{$BFQarr[0]}}{'shared'} and ${${$BFQreadsum}{$BFQarr[0]}}{'shared'} > 0) {
			$BFQis_shared=1;
		}
		my $BFQis_allelic=0;
		my @BFQthisread_gens=();
		if (exists ${$BFQallelicreads}{$BFQarr[0]}) {
			$BFQis_allelic=1;
			my @BFQtemp_gens=();
			if (exists ${${$BFQallelicreads}{$BFQarr[0]}}{'gen'}) {
				@BFQtemp_gens=keys %{${${$BFQallelicreads}{$BFQarr[0]}}{'gen'}};
#				print "${BFQsubinfo}Test: Geno: @BFQtemp_gens Read: $BFQarr[0]\n";### For test ###
			}
			else {
				print STDERR "${BFQsubinfo}Warnings: $BFQarr[0] no 'gen' key in BFQallelicreads\n";
				$excludedreads{$BFQarr[0]}++;
				print BAMOUT $BFQline1."\n";
				next;
			}
			foreach (@BFQtemp_gens) {
				push (@BFQthisread_gens, $_) if (exists $BFQgenos{$_});
			}
			if (scalar(@BFQthisread_gens)<1 or scalar(@BFQthisread_gens)>3) {
#				print STDERR "${BFQsubinfo}Warnings: no subgenome assignments for read $BFQarr[0] at Chrom: Pos $BFQarr[2]:$BFQarr[3]\n";### For test ###
				$excludedreads{$BFQarr[0]}++;
				print BAMOUT $BFQline1."\n";
				next;
			}
		}
		
		if ($BFQis_shared==0 and $BFQis_allelic==0) {
#			print STDERR "${BFQsubinfo}Warnings: read ($BFQarr[0]) is neither shared or allelic at Chrom:Pos $BFQarr[2]:$BFQarr[3]\n";
			print BAMOUT $BFQline1."\n";
			$excludedreads{$BFQarr[0]}++;
			next;
		}
		elsif ($BFQis_shared==1 and $BFQis_allelic==1) {
			print STDERR "${BFQsubinfo}Warnings: read ($BFQarr[0]) is both shared and allelic at Chrom:Pos $BFQarr[2]:$BFQarr[3]\n";
			print BAMOUT $BFQline1."\n";
			$excludedreads{$BFQarr[0]}++;
			next;
		}
#check if mapped
		if ($BFQmapcode==1) {
			if ($BFQarr[1] & 0x0004) {
				print BAMOUT $BFQline1."\n";
				next;
			}
		}
		elsif ($BFQmapcode==2) {
			unless ($BFQarr[1] & 0x0004) {
				print BAMOUT $BFQline1."\n";
				next;
			}
		}
##check if MAPQ threshold
		next unless (defined $BFQarr[4] and $BFQarr[4]>=$BFQmapq);
##Get readgroup info
		my $BFQreadgroup='undefRG';
		if ($BFQrg_aware==1) {
			if ($BFQline1=~/\tRG:Z:(\S+)/) {
				$BFQreadgroup=$1;
			}
		}
##check if mapped to reverse strand
		my $BFQread_id=$BFQarr[0];
		my $BFQreadseq=$BFQarr[9];
		my $BFQreadqual=$BFQarr[10];
		if ($BFQarr[1] & 0x0010) {###0x0010 = reverse strand
			$BFQreadseq=reverse ($BFQreadseq);
			$BFQreadseq=~tr/ATCGatcg/TAGCtagc/;
			$BFQreadqual=reverse ($BFQreadqual);
		}
##check if one of pair and check R1 or R2
		if ($BFQarr[1] & 0x0001) {
			if ($BFQarr[1] & 0x0040) {###read1
				if ($BFQis_shared==1) {
					@{${${${$BFQreadsum}{$BFQread_id}}{'reads'}}{1}}=($BFQreadseq, $BFQreadqual);
					${${${${$BFQreadpair{$BFQarr[2]}}{$BFQreadgroup}}{$BFQread_id}}{1}}{'shared'}++;
				}
				elsif ($BFQis_allelic==1) {
					@{${${${$BFQallelicreads}{$BFQread_id}}{'reads'}}{1}}=($BFQreadseq, $BFQreadqual);
					foreach (@BFQthisread_gens) {
						${${${${$BFQreadpair{$BFQarr[2]}}{$BFQreadgroup}}{$BFQread_id}}{1}}{$_}++;
					}
				}
			}
			elsif ($BFQarr[1] & 0x0080) {###read2
				if ($BFQis_shared==1) {
					@{${${${$BFQreadsum}{$BFQread_id}}{'reads'}}{2}}=($BFQreadseq, $BFQreadqual);
					${${${${$BFQreadpair{$BFQarr[2]}}{$BFQreadgroup}}{$BFQread_id}}{2}}{'shared'}++;
				}
				elsif ($BFQis_allelic==1) {
					@{${${${$BFQallelicreads}{$BFQread_id}}{'reads'}}{2}}=($BFQreadseq, $BFQreadqual);
					foreach (@BFQthisread_gens) {
						${${${${$BFQreadpair{$BFQarr[2]}}{$BFQreadgroup}}{$BFQread_id}}{2}}{$_}++;
					}
				}
			}
			else {
				print STDERR "${BFQsubinfo}Warnings: unknown R1 or R2 (FLAG: $BFQarr[1]) at line $BFQnumline (Readid: $BFQarr[0]) in BAM $BFQbamin\n";
				next;

			}
		}
		else {####sing-end
			if ($BFQis_shared==1) {
				@{${${${$BFQreadsum}{$BFQread_id}}{'reads'}}{'U'}}=($BFQreadseq, $BFQreadqual);
#				${${${${$BFQreadpair{$BFQarr[2]}}{$BFQreadgroup}}{$BFQread_id}}{'U'}}{'shared'}++; ### Ignore single read for FPKM
			}
			elsif ($BFQis_allelic==1) {
				@{${${${$BFQallelicreads}{$BFQread_id}}{'reads'}}{'U'}}=($BFQreadseq, $BFQreadqual);
#				foreach (@BFQthisread_gens) {
#					${${${${$BFQreadpair{$BFQarr[2]}}{$BFQreadgroup}}{$BFQread_id}}{'U'}}{$_}++; ### Ignore single read for FPKM
#				}
			}
		}
	}
	close BAMIN;
	close BAMOUT;
	print "${BFQsubinfo}Warnings: excluded reads num: ".scalar(keys %excludedreads)."\n";

###Output fastq (shared)
	my $BFQshared_fastq=$BFQoutfqprefix.'.shared.fastq';
	unlink $BFQshared_fastq if (-e $BFQshared_fastq);
	close SHAREOUT if (defined fileno(FQOUT));
	unless (open (SHAREOUT, ">$BFQshared_fastq")) {
		print STDERR "${BFQsubinfo}Error: can not write Fastq: $BFQshared_fastq\n";
		return $BamKit_failure;
	}
	foreach my $BFQreadid_shared (keys %{$BFQreadsum}) {
		if (${${$BFQreadsum}{$BFQreadid_shared}}{'shared'} and defined ${${$BFQreadsum}{$BFQreadid_shared}}{'shared'} and ${${$BFQreadsum}{$BFQreadid_shared}}{'shared'}>0 and exists ${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}) {
			if (exists ${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{1}) {
				my ($BFQseq, $BFQqual)=@{${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{1}};
				print SHAREOUT '@'.$BFQreadid_shared."/1\n$BFQseq\n+\n$BFQqual\n";
			}
			if (exists ${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{2}) {
				my ($BFQseq, $BFQqual)=@{${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{2}};
				print SHAREOUT '@'.$BFQreadid_shared."/2\n$BFQseq\n+\n$BFQqual\n";
			}
			if (exists ${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{'U'}) {
				my ($BFQseq, $BFQqual)=@{${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{'U'}};
				print SHAREOUT '@'.$BFQreadid_shared."\n$BFQseq\n+\n$BFQqual\n";
			}
		}
	}
	close SHAREOUT;
	foreach my $BFQind_gen (keys %BFQgenos) {
		my $BFQname=$BFQoutfqprefix.'.'.$BFQind_gen.'.fastq';
		unlink $BFQname if (-e $BFQname);
		if (-s $BFQshared_fastq) {
			unless (CopyFile($BFQshared_fastq, $BFQname)) {
				print STDERR "${BFQsubinfo}Error: copy file failed: $BFQname\n";
			}
		}
		close FQOUT if (defined fileno(FQOUT));
		unless (open (FQOUT, ">>$BFQname")) {
			print STDERR "${BFQsubinfo}Error: can not write Fastq: $BFQname\n";
			return $BamKit_failure;
		}
		foreach my $BFQreadid_allelic (keys %{$BFQallelicreads}) {
			if (exists ${${$BFQallelicreads}{readid_allelic}}{'reads'}) {
				if (exists ${${${$BFQallelicreads}{readid_allelic}}{'reads'}}{1}) {
					my ($BFQseq, $BFQqual)=@{${${${$BFQallelicreads}{readid_allelic}}{'reads'}}{1}};
					print FQOUT '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
				}
				if (exists ${${${$BFQallelicreads}{readid_allelic}}{'reads'}}{2}) {
					my ($BFQseq, $BFQqual)=@{${${${$BFQallelicreads}{readid_allelic}}{'reads'}}{2}};
					print FQOUT '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
				}
				if (exists ${${${$BFQallelicreads}{readid_allelic}}{'reads'}}{'U'}) {
					my ($BFQseq, $BFQqual)=@{${${${$BFQallelicreads}{readid_allelic}}{'reads'}}{'U'}};
					print FQOUT '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
				}
			}
		}
		close FQOUT;
		if ( -s $BFQname) {
			$returnhash{$BFQind_gen}=$BFQname;
		}
		elsif (-e $BFQname) {
			$returnhash{$BFQind_gen}=0;
		}
	}
### Paired reads Counts
##Format: %BFQreadpair=(chr => (RG1 => (readid => 1 => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
##Format: %BFQreadcounts=(chr => ('chr' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
	my %BFQreadcounts=();
	foreach my $BFQchrom (keys %BFQreadpair) {
		foreach my $BFQrg (keys %{$BFQreadpair{$BFQchrom}}) {
			foreach my $BFQreadid (keys %{${$BFQreadpair{$BFQchrom}}{$BFQrg}}) {
				if (exists ${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1} and exists ${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{2}) {
					${$BFQreadcounts{$BFQchrom}}{'chr'}++;
					${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQrg}}{'global'}++;
					if (exists ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{'shared'} and ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{'shared'} >0) {###shared
						foreach (keys %BFQgenos) {
							${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQrg}}{$_}++;
							${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_}++;
						}
					}
					else {##allelic
						foreach (keys %BFQgenos) {
							if (exists ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{$_} and ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{$_} >0){
								${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQrg}}{$_}++;
								${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_}++;
							}
						}
					}
				}
			}
		}
	}
	
### Calculate 
	foreach my $BFQchrom (keys %BFQreadcounts) {
		if (exists ${$BFQreadcounts{$BFQchrom}}{'chr'}) {
			if (exists $BFQreflen{$BFQchrom}and $BFQreflen{$BFQchrom}>0) {
				${$BFQfpkms{$BFQchrom}}{'chr'}=${$BFQreadcounts{$BFQchrom}}{'chr'} * 10^9 / ($BFQreflen{$BFQchrom} *$totalreads_AABBDD);
			}
			else {
				print STDERR "${BFQsubinfo}Error: can not get ref ($BFQchrom) length from header of BAM: $BFQbamin\n";
				${$BFQfpkms{$BFQchrom}}{'chr'}='undef';
			}
		}
		if (exists ${$BFQreadcounts{$BFQchrom}}{'allelic'}) {
			foreach ('A', 'B', 'D') {
				if (exists ${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_}) {
					${${$BFQfpkms{$BFQchrom}}{'allelic'}}{$_}= ${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_} * 10^9 / ($BFQreflen{$BFQchrom} *$totalreads_AABBDD);
				}
				else {
					${${$BFQfpkms{$BFQchrom}}{'allelic'}}{$_}=0;
				}
			}
		}
		if (exists ${$BFQreadcounts{$BFQchrom}}{'read_group'}) {
			foreach my $BFQindrg (keys %reads_per_rg) {
				if (exists ${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQindrg}) {
					foreach ('global', 'A', 'B', 'D') {
						if (exists ${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQindrg}}{$_}) {
							${${${$BFQfpkms{$BFQchrom}}{'read_group'}}{$BFQindrg}}{$_}=${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQindrg}}{$_} * 10^9 / ($BFQreflen{$BFQchrom} * $reads_per_rg{$BFQindrg});
						}
						else {
							${${${$BFQfpkms{$BFQchrom}}{'read_group'}}{$BFQindrg}}{$_}=0;
						}
					}
				}
				else {
					foreach ('global', 'A', 'B', 'D') {
						${${${$BFQfpkms{$BFQchrom}}{'read_group'}}{$BFQindrg}}{$_}=0;
					}
				}
			}
		}
	}
	%BFQreadcounts=();
	return ($BamKit_success, \%returnhash, \%BFQfpkms);
}
