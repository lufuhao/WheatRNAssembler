#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use Bio::DB::Sam;
use Scalar::Util;
use File::Path qw/remove_tree rmtree/;
use File::Copy qw/move/;
use FuhaoPerl5Lib::BamKit qw/ReduceReadNameLength SamCleanHeader SortBam IndexBam CalcFPKM ReadSam/;
use FuhaoPerl5Lib::VcfKit qw/ExtractVcf ReadVcf RunFreebayes ReadVariantType HapcompassVcf HcFragments RunHapCompass GroupFragments CorrectAlleles ReadHcOut/;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::FileKit;
use FuhaoPerl5Lib::FastaKit qw/CdbFasta CdbYank IndexFasta RenameFasta RunFqTrinity/;
use FuhaoPerl5Lib::MiscKit qw/MaxLength/;
use Storable qw/dclone/;
use Data::Dumper qw/Dumper/;
#no autovivification;
use constant USAGE=><<EOH;


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

EOH



###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $numthreads, $debug, $ver, $cleanlog);
#Global
my ($reference, $file_cluster, $list);
my ($file_bam_aabbdd);#, $file_bam_aabb, $file_bam_aa, $file_bam_dd);###filelist
my ($file_vcf_aabbdd, $file_fpkm, $file_alleles);
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
#Trinity
my ($path_trinity, $num_cpus);
#CDHIT
my ($path_cdhitest);
#Bowtie2
#my ($pathbowtie2build, $pathbowtie2p);
my $fq_qual_score;
#LOG
my ($cluster_log, $fpkm_log, $allele_log, $geno_log, $clusterfpkm);###LOG
my ($totalreads_AABBDD, $totalreads_AABB, $totalreads_AA, $totalreads_DD);

GetOptions(
	"help|h!" => \$help,
	"reference|r:s" => \$reference,
	"cluster|c:s" => \$file_cluster,
	"list:s" => \$list,
	"bam1:s" => \$file_bam_aabbdd,
#	"bam2:s" => \$file_bam_aabb, 
#	"bam3:s" => \$file_bam_aa, 
#	"bam4:s" => \$file_bam_dd,
	"vcf1:s" => \$file_vcf_aabbdd, 
	"fpkm:s" => \$file_fpkm,
	"allele:s" => \$file_alleles,
	"taggenome:s" => \$bam_genome_tag,
	"tagchrom:s" => \$bam_chromo_tag,
	"numreads:s" => \$file_totalreads,
	"vardelimiter:s" => \$geno_delimiter,
	"cdbfastapath:s" => \$path_cdbfasta,
	"cdbyankpath:s" => \$path_cdbyank,
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
	"vcfmergepath:s" => \$path_vcfmerge,
	"vcfconcatpath:s" => \$path_vcfconcat,
	"vcfsortpath:s" => \$path_vcfsort,
	"trinitypath:s" => \$path_trinity,
#	"cdhitestpath:s" => \$path_cdhitest,
	"logcluster:s" => \$cluster_log,
	"logfpkm:s" => \$fpkm_log,
	"logcfpkm:s" => \$clusterfpkm,
	"logallele:s" => \$allele_log,
	"loggeno:s" => \$geno_log,
	"clean!" => \$cleanlog,
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
#Trinity
$path_trinity='Trinity' unless (defined $path_trinity);
my $trinity_addcmd='--max_memory '.$javamem.' --run_as_paired --CPU '.$num_cpus.' --group_pairs_distance '.$maxinsert.' --full_cleanup --min_kmer_cov 3 --min_glue 3';
#CDHITEST
$path_cdhitest='cd-hit-est' unless (defined $path_cdhitest);
my $cdhit_addcmd=' -c 1.00 -n 10 -T 0 -r 1 -d 0 -M 30000 ';
#LOG
$cluster_log='0.cluster.log' unless (defined $cluster_log);
$fpkm_log='0.reference.fpkm.log' unless (defined $fpkm_log);
$allele_log='0.allele.log' unless (defined $allele_log);
$geno_log='0.geno.log' unless (defined $geno_log);
$clusterfpkm='0.cluster.FPKM.log' unless (defined $clusterfpkm);
$cleanlog=0 unless (defined $cleanlog);


### input and output ################################################
#Check important file inpput
print "Checking input and output\n";### For test ###
die "Error: undefined reference file\n" unless (defined $reference);
($reference)=AddFilePath($reference);
die "Error: not-found reference file\n" unless (-s $reference);
die "Please specify cluster file\n" unless (defined $file_cluster and -s $file_cluster);
die "Please specify the file of bam list for AABBDD\n" unless (defined $file_bam_aabbdd and -s $file_bam_aabbdd);
#die "Please specify the file of bam list for AABB\n" unless (defined $file_bam_aabb and -s $file_bam_aabb);
#die "Please specify the file of bam list for AA\n" unless (defined $file_bam_aa and -s $file_bam_aa);
#die "Please specify the file of bam list for DD\n" unless (defined $file_bam_dd and -s $file_bam_dd);
die "Please specify the file for number of reads\n" unless (defined $file_totalreads and -s $file_totalreads);
die "Please specify the file for alleles\n" unless (defined $file_alleles and -s $file_alleles);
die "Please specify the file for FPKM\n" unless (defined $file_fpkm and -s $file_fpkm);

#Remove last-run cluster log 
unlink $cluster_log if (-e $cluster_log);
unlink $fpkm_log if (-e $fpkm_log);

#read AABBDD bam files
my $bamAABBDDfiles=' '; 
my $bamAABBfiles=' ';
my $bamAAfiles=' ';
my $bamDDfiles=' ';
my $numbams=0;
my @bam_AABBDD=();
if ($file_bam_aabbdd=~/\.bam$/i) {
	push (@bam_AABBDD, $file_bam_aabbdd);
}
else {
	close BAMAABBDD if (defined fileno(BAMAABBDD));
	open (BAMAABBDD, "<$file_bam_aabbdd") || die "InputOutputError: can not open AABBDD $file_bam_aabbdd\n";
	while (my $line1=<BAMAABBDD>) {
		chomp $line1;
		if (defined $line1 and $line1 =~/^\S+\.bam$/i and -s $line1) {
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
$bamAABBDDfiles=join(',', @bam_AABBDD);
$numbams=scalar(@bam_AABBDD);
undef @bam_AABBDD;

###format: %expressfpkm{chr}=(abd_fpkm, ab_fpkm, a_fpkm, d_fpkm)
my %expressfpkm=();
close FPKMIN if (defined fileno(FPKMIN));
unless (open(FPKMIN, "< $file_fpkm")) {
	die "Error: not read FPKM file\n"
}
while (my $line=<FPKMIN>) {
	chomp $line;
	my @fpkmarr=split(/\t/, $line);
	my $chrom=shift @fpkmarr;
	@{$expressfpkm{$chrom}}=@fpkmarr;
}
close FPKMIN;

my $using_existing_alleles=0;
###format: %ancestralalleles{chr}{$pos}=('A' => 0/1, 'B' => 1/0, 'D' => 1/0)
my %ancestralalleles=();
if (defined $file_alleles and -s $file_alleles) {
	$using_existing_alleles=1;
	close EXISTINGALLELE if (defined fileno(EXISTINGALLELE));
	unless (open (EXISTINGALLELE, "< $file_alleles")) {
		die "Error: can not read fixed allele file\n";
	}
	while (my $line=<EXISTINGALLELE>) {
		chomp $line;
		my @arr=split(/\t/, $line);
		$ancestralalleles{$arr[0]}{$arr[1]}{'A'}=$arr[2];
		$ancestralalleles{$arr[0]}{$arr[1]}{'B'}='?';#$arr[3];
		$ancestralalleles{$arr[0]}{$arr[1]}{'D'}=$arr[4];
	}
	close EXISTINGALLELE;
}


###input vcf
if (defined $file_vcf_aabbdd) {
	if (! -s $file_vcf_aabbdd) {
		die "MainError: can not find AABBDD file: $file_vcf_aabbdd\n";
	}
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
if ($verbose) {
	print "InputOutputReport: VCF inputs: \n";
	if (defined $file_vcf_aabbdd) {print "AABBDD VCF1: $file_vcf_aabbdd\n";} else {print "AABBDD VCF1: null\n";}
}



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





### Main ############################################################
my $RunDir=getcwd;
#Prepare temporary folder
if (! -d "$RunDir/AABBDD") {
	mkdir ("$RunDir/AABBDD", 0766) || die "(Main0)Error: can not create folder $RunDir/AABBDD\n";
}
unlink glob "$RunDir/AABBDD/*";
if (! -d "$RunDir/fasta") {
	mkdir ("$RunDir/fasta", 0766) || die "(Main0)Error: can not create folder $RunDir/fasta\n";
}
if (! -d "$RunDir/kmergenie") {
	mkdir ("$RunDir/kmergenie", 0766) || die "(Main0)Error: can not create folder $RunDir/kmergenie\n";
}

#if (! -d "$RunDir/AABB") {
#	mkdir ("$RunDir/AABB", 0766) || die "(Main0)Error: can not create folder AABB\n";
#}
#unlink glob "$RunDir/AABB/*";
#if (! -d "$RunDir/AA") {
#	mkdir ("$RunDir/AA", 0766) || die "(Main0)Error: can not create folder AA\n";
#}
#unlink glob "$RunDir/AA/*";
#if (! -d "$RunDir/DD") {
#	mkdir ("$RunDir/DD", 0766) || die "(Main0)Error: can not create folder DD\n";
#}
#unlink glob "$RunDir/DD/*";

#(Main1) index fasta
my ($test_cdbfasta, $fasta_index)=&CdbFasta($reference, $path_cdbfasta);
if ($test_cdbfasta == $cmd_failure_code ) {
	die "(Main1)Error: can not index fasta: $reference\n";
}
print "Fasta Index: $fasta_index\n";
# Load clustered fasta id list, 1 cluster perl line, saparated by space
open (CLUSTER, "<$file_cluster") || die "(Main1)Error: can not open file: $file_cluster\n";
# Output cluster succeeds or not. Format: 
open (CLUSTERLOG, ">$cluster_log") || die "(Main1)Error: can not write file: $cluster_log\n";
# output fpkm calculated
open (REFERFPKM, ">$fpkm_log") || die "(Main1)Error: can not write file: $fpkm_log\n";
# output fpkm for each cluster
open (CLUSTFPKM, ">$clusterfpkm") || die "(Main1)Error: can not write file: $clusterfpkm\n";
# output allele log
open (ALLELELOG, ">$allele_log") || die "(Main1)Error: can not write file: $allele_log\n";
open (GENOLOG, ">$geno_log") || die "(Main1)Error: can not geno log: $geno_log\n";
my $cluster_num=0;###for quick find which line/cluster has problem
my @cluster_seqids=();###
while (my $cluster_line=<CLUSTER>) {
	chomp $cluster_line;
	$cluster_num++;

### Step0: control which lines to reads: --list
	my $stage=0;
#	print "(Main$stage)Step: Line control\n"; ### For test ###
	if ($cluster_num<$cluster_start) {
		next;
	}
	elsif ($cluster_num>=$cluster_start) {
		if (defined $cluster_end and $cluster_num>$cluster_end) {
			next;
		}
	}
	next if ($cluster_line=~/^#/);
	print GENOLOG "\n\n\nCluster$cluster_num: $cluster_line\n";
	print "\n\n\n\n\n##### Prcessing Cluster $cluster_num ###\n$cluster_line\n";
	print STDERR "\n\n\n\n\n##### Prcessing Cluster $cluster_num ###\n$cluster_line\n";
	@cluster_seqids=(); ###Empty this in case of abnormal duplicates
	@cluster_seqids=split(/\s+/, $cluster_line);
##COMMENT: Check if empty line
	if (scalar(@cluster_seqids)<1) {
		print STDERR "(Main2)Warnings: line $cluster_num in $file_cluster ignored as empty\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tNumIDinCluster<1\n";
		next;
	}
	my $fastaids=join(',', @cluster_seqids);
##COMMENT: Empty folder
	unless (chdir $RunDir) {
		print STDERR "(Main$stage)Error: can not chdir to : $RunDir\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tChDirAtBegining\n";
		next;
	}
	unlink glob "$RunDir/AABBDD/*";
	unlink glob "$RunDir/AABB/*";
	unlink glob "$RunDir/AA/*";
	unlink glob "$RunDir/DD/*";
	if ($cleanlog) {
		my @cluster=glob "$RunDir/Clust*";
		foreach my $indfolder (@cluster) {
			print "(Main$stage)Testing: delete PATH: $indfolder\n";
#			sleep(1);
			unlink glob "$indfolder/*";
		}
	}



## Stage1: extract fasta reference
	$stage=1;
	print "(Main$stage)Step: Extract AABBDD fasta\n"; ### For test ###
	my $clusterrefer="$RunDir/Clust$cluster_num/ref.$cluster_num.fa";
	if ( ! -d "$RunDir/Clust$cluster_num") {
		unless (mkdir ("$RunDir/Clust$cluster_num", 0766)) {
			print STDERR "(Main$stage)Error: create folder $RunDir/Clust$cluster_num\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tMkdirCluster\n";
			next;
		}
	}
	unless (-s $clusterrefer) {
		unless (&CdbYank($fasta_index, $clusterrefer, \@cluster_seqids, $path_cdbyank)) {
			print STDERR "(Main$stage)Error: failed to extract ClusterID $cluster_num: @cluster_seqids\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tCdbYank\n";
			next;
		}
	}
	unless (-s "$clusterrefer.fai") {
		unless (&IndexFasta($clusterrefer, $path_samtools)) {
			print STDERR "(Main$stage)Error: failed to index ClusterID $cluster_num: @cluster_seqids\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tIndexFasta\n";
			next;
		}
	}




###Stage2: Extract bam
	$stage=2;
	print "(Main$stage)Step: Extract AABBDD BAMs\n"; ### For test ###
	my $cluter_oribam="$RunDir/Clust$cluster_num/AABBDD.$cluster_num.ori.bam";
	unless (&exec_cmd_return("bam_multiextract.sh -b $bamAABBDDfiles -s $fastaids -n $numbams -o $cluter_oribam")) {
		print STDERR "(Main$stage)Error: AABBDD BAM extraction failed\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tBAMExtract\n";
		next;
	}


###Stage3: Extract VCF
	$stage=3;
	print "(Main$stage)Step: Extract AABBDD VCF\n"; ### For test ###
	my $cluster_orivcf="$RunDir/Clust$cluster_num/AABBDD.$cluster_num.ori.vcf.gz";
	if (defined $file_vcf_aabbdd and -s $file_vcf_aabbdd) {
		unless (&ExtractVcf($file_vcf_aabbdd, $cluster_line, $cluster_orivcf)) {
			print STDERR "(Main$stage)Error: AABBDD VCF extraction failed\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tVcfExtractAABBDD\n";
			next;
		}
	}
	else {
		print STDERR "(Main$stage)Error: NoVcfAABBDD\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tNoVcfAABBDD\n";
		next;
	}
	
	


###Stage4: Check Point: detect polymorphism, continue if no, region and merge if yes
	$stage=4;
	my $total_sites=0;
	my $num_allelic_sites=0;
	my $num_not_allelic=0;
#Format: %{$polymorphicsites}=($chr => ($posit => ++));
	my $polymorphicsites={};
	my ($test_loadvcf, $aabbdd_vcf_obj)=&ReadVcf($cluster_orivcf, 20);
	unless ($test_loadvcf) {
		print STDERR "(Main$stage)Error: loadvcf $cluster_orivcf failed\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tLoadVCF\n";
		next;
	}
	foreach my $chrom (keys %{$aabbdd_vcf_obj}) {
		foreach my $posit (keys %{${$aabbdd_vcf_obj}{$chrom}}) {
			$total_sites++;
			if (exists ${$aabbdd_vcf_obj}{$chrom}{$posit}[2]) {
				my $thisgeno=${$aabbdd_vcf_obj}{$chrom}{$posit}[2];
				if ($thisgeno=~/^\d+\/\d+\/\d+$/) {
					my @arr=split(/\//, $thisgeno);
					my %hash1=();
					map {$hash1{$_}++} @arr;
#					print "(Main$stage)Test: Chr:Pos $chrom:$posit Geno: $thisgeno NumAllele: ".scalar(keys %hash1)."\n"; ### For test ###
					if (scalar(keys %hash1)==1) {
						$num_not_allelic++;
					}
					elsif (scalar(keys %hash1)>1) {
						map {${$polymorphicsites}{$chrom}{$posit}{$_}++} (keys %hash1);
						$num_allelic_sites++;
					}
					else {
						print STDERR "(Main$stage)Error: unknown number of genotypes\n";
						print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tNumGeno\n";
						next;
					}
				}
				else {
					print STDERR "(Main$stage)Error: error genotypes\n";
					print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tErrorGeno\n";
					next;
				}
			}
			else {
				print STDERR "(Main$stage)Error: NoGeno\n";
				print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tNoGeno\n";
				next;
			}
		}
	}
	if (1) {
		print "(Main$stage)Info: Total sites: $total_sites\n";
		print "(Main$stage)Info: Allelic sites: $num_allelic_sites\n";
		print "(Main$stage)Info: NotAllelic sites: $num_not_allelic\n";
	}



##Stage5: Defined important hash or array here
	$stage=5;
	my $assignallele={};
#	%{$polymorphismreadids}=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
	my $polymorphismreadids={};
#Format: %{$fragments}=(readid => (chrom => (pos => (geno1 => quality, geno2 => quality)));
	my $fragments={};
#Format: %{$readidsum}=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
#								'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#								'shared' => ++/0,
#								'mapped' => ++/0,
#								'ref' => chr => ++));
	my $readsum={};
#my $test_aa_expressed=0; my $test_bb_expressed=0; my $test_dd_expressed=0;



###Stage6: Group alleles into subgenome
	$stage=6;
	print "\n(Main$stage)Step: GroupReads\n"; ### For test ###
##ReadSam
	my ($test_readbam, $aabbdd_bam_obj)=&ReadSam($cluter_oribam, $clusterrefer, 1);
	unless ($test_readbam) {
		print STDERR "(Main$stage)Error: Reading AABBDD BAM failed\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\treadsamAABBDD\n";
		next;
	}
##GroupReads
	(my $test_groupreads, $fragments, $readsum)=&GroupReads($aabbdd_bam_obj, $aabbdd_vcf_obj);
	if ($test_groupreads) {
		print STDERR "(Main$stage)Error: GroupReads SUB failed\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tGroupReads\n";
		next;
	}
	if (0) {
		print "(Main$stage)Test: \$readsum\n";
		print Dumper $readsum;
		print "\n";
	}
	if (0) { ### test $fragments ###
		print "(Main$stage)Test: \$fragments\n";
		print Dumper $fragments;### For test ###
		print "\n";
	}
	$aabbdd_bam_obj={};



###Stage6: checkpoint
	$stage=6;
	my $fragmentfile="$RunDir/Clust$cluster_num/AABBDD.$cluster_num.fragments";
##Format: %{$corrected_fragments}=($readid => chrom => $partstart => $pos => $allele => $qual);
##Format: %{$corrected_sites}=($chr => $pos => $allele => ++)
##Format: %{}
	my ($test_hc, $corrected_fragments, $corrected_sites, $vcfreferencer)=HcFragments($fragments, $fragmentfile, $readsum);
	unless ($test_hc) {
		print STDERR "(Main$stage)Error: prepare HapCompass Fragments error\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tFragments\n";
		next;
	}
	print "(Main$stage)Info: Number of sites to be phased: $num_allelic_sites => ";
	my $corrected_num_allelic_sites=0;
	foreach (keys %{$corrected_sites}) {
		$corrected_num_allelic_sites+=scalar(keys %{${$corrected_sites}{$_}});
	}
#	print "(Main$stage)Test: $corrected_num_allelic_sites\n"; ### for test ###
	print "(Main$stage)Info: Fragments cleaning: ".scalar(keys %{$fragments})." => ".scalar(keys %{$corrected_fragments})."\n";
	
#	unless (GroupFragments($corrected_fragments)) {
#		print STDERR "(Main$stage)Error: GroupFragments running";
#		next;
#	}
#	print "(Main$stage)Test: \%{\$corrected_sites}\n"; print Dumper $corrected_sites; print "\n";



##Stage7: decide ploidy here
	$stage=7;
	print "\n(Main$stage)Step: Ploidy\n"; ### For test ###
	my $test_run_phase=0;
	my $assembly_desc='';
	my $ploidy=1;
	my %genomecount=();
	my %chromcount=();
	my $final_geno={};
	my $totalmappedreads=0;
	my $total_notag=0;
	my $total_tagged=0;
	my @final_ploidy=();
	my $inferred_genoploidy='';
	my $inferred_chrom='';
	my $genomeinfer_from_reads=0;
	my $unknown_genome=0;
	my %chrom2genome=();
	my $guessgenome='';
	my $num_excluded_reads=0;
	my $num_shared_reads=0;
	my ($test_aa_expressed, $test_bb_expressed, $test_dd_expressed)=(0, 0, 0); 
	foreach my $ind_read (keys %{$readsum}) {
		if (exists ${$readsum}{$ind_read}{'excluded'} and ${$readsum}{$ind_read}{'excluded'}>0) {
			$num_excluded_reads++;
			${$readsum}{$ind_read}{'shared'}=0;
			next;
		}
		if (exists ${$readsum}{$ind_read}{'shared'} and ${$readsum}{$ind_read}{'shared'}>0) {
			$num_shared_reads++;
		}
		if (exists ${${$readsum}{$ind_read}}{'gen'}) {
			foreach my $xgen (keys %{${${$readsum}{$ind_read}}{'gen'}}) {
				$genomecount{$xgen}++;
				if (exists ${$readsum}{$ind_read}{'ref'}) {
					foreach my $xref (keys %{${$readsum}{$ind_read}{'ref'}}) {
						$chrom2genome{$xref}{$xgen}++;
					}
				}
			}
			$total_tagged++;
		}
		else {
			$total_notag++;
		}
		
		if (exists ${${$readsum}{$ind_read}}{'chr'}) {
			foreach (keys %{${${$readsum}{$ind_read}}{'chr'}}) {
				$chromcount{$_}++;
			}
		}
		if (exists ${${$readsum}{$ind_read}}{'mapped'}) {
			$totalmappedreads++;
		}
	}
	unless ($totalmappedreads>0 and ($num_shared_reads+scalar(keys %{$corrected_fragments}))>10) {### minimum 10 reads
		print STDERR "(Main$stage)Error: Not enough reads\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tNotEnoughReads\n";
		next;
	}
	if (($total_tagged/$totalmappedreads) >=0.4) {###Decide ploidy if 40% reads assigned
		$genomeinfer_from_reads=1;
		foreach (keys %genomecount) {
			$inferred_genoploidy.=$_ if (($genomecount{$_}/$total_tagged) >=0.1);
			push (@final_ploidy, $_);
		}
		print "(Main$stage)Info: Inferred geno: $inferred_genoploidy\n";
		unless ($inferred_genoploidy =~/^[ABD]+$/) {
			print STDERR "(Main$stage)Error: unwanted \$inferred_genoploidy: $inferred_genoploidy\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tinferred_genoploidy\n";
			next;
		}
		$ploidy=length($inferred_genoploidy);
		my $num_chr_assigned=0;
		map {my $gen=substr($_, 1, 1); if ($inferred_genoploidy=~/$gen/) {$num_chr_assigned++;}} (keys %chromcount);
		print "(Main$stage)Info: Num_InferChrom: $num_chr_assigned\n";###test###
		if ($num_chr_assigned>=1) {
			foreach (keys %chromcount) {
				my $gen=substr($_, 1, 1); 
				next unless ($inferred_genoploidy=~/$gen/);
				if ($chromcount{$_} >= ($totalmappedreads/$num_chr_assigned)) {
					$inferred_chrom.="$_($chromcount{$_})";
					print "(Main$stage)Info: ChromInfer: $_\t$chromcount{$_}\n";
				}
			}
			print "(Main$stage)Info: Inferred geno: $inferred_chrom\n";
		}
	}
	else {### Decide ploidy from diploid plants
		$genomeinfer_from_reads=0;
		###format: %expressfpkm{chr}=(abd_fpkm, ab_fpkm, a_fpkm, d_fpkm)
		foreach (@cluster_seqids) {
			if (exists $expressfpkm{$_}) {
				$test_aa_expressed=1 if (defined $expressfpkm{$_}[2] and $expressfpkm{$_}[2]>$minimum_fpkm);
				$test_dd_expressed=1 if (defined $expressfpkm{$_}[3] and $expressfpkm{$_}[3]>$minimum_fpkm);
				$test_bb_expressed=1 if (exists $genomecount{'B'} and ($genomecount{'B'}/$total_tagged) >=0.05);
			}
		}
		if ($corrected_num_allelic_sites==0) {
			$ploidy=1;
			$test_run_phase=0;
			my $maxcount=0;
			foreach (keys %genomecount) {
				if ($genomecount{$_}> $maxcount) {
					$guessgenome=$_;
					$maxcount=$genomecount{$_};
				}
				elsif ($genomecount{$_}== $maxcount) {
					$guessgenome.=$_;
				}
			}
			$inferred_genoploidy=$guessgenome;
		}
		elsif ($corrected_num_allelic_sites>0) {
			if ($corrected_num_allelic_sites==1) {
				$test_run_phase=0;
				$unknown_genome=1;
			}
			else {
				$test_run_phase=1;
			}
			$ploidy=$test_aa_expressed+$test_dd_expressed+$test_bb_expressed;
			$inferred_genoploidy.='A' if ($test_aa_expressed>0);
			$inferred_genoploidy.='B' if ($test_bb_expressed>0);
			$inferred_genoploidy.='D' if ($test_dd_expressed>0);
			$inferred_genoploidy='ABD' if ($ploidy==0 or $inferred_genoploidy eq '');
		}
	}###PLOIDY
	unless ($inferred_genoploidy=~/^[ABD]+$/) {
		print STDERR "(Main$stage)Error: genoploidy inferring error: $inferred_genoploidy\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tinferred_genoploidy\n";
		next;
	}
	unless ($ploidy>=0 and $ploidy<4) {
		print STDERR "(Main$stage)Error: ploidy inferring error: $ploidy\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tinferredploidy:$ploidy:$inferred_genoploidy:$genomeinfer_from_reads\n";
		next;
	}
	if (1) {### print out genome assign results ###
		print GENOLOG "##### Genome summary #####\nTotal mapped reads: $totalmappedreads\nTagged: $total_tagged\nNotag: $total_notag\nShared reads: $num_shared_reads\nExcluded reads: $num_excluded_reads\n\n";
		
		map {print GENOLOG "Genome: $_\tCount: $genomecount{$_}\n"} (sort (keys %genomecount)); 
		print GENOLOG "\n";
		map {print GENOLOG "Chrom: $_\tCount: $chromcount{$_}\n"} (sort (keys %chromcount));
		print GENOLOG "##### Genome summary #####\n";
		print GENOLOG "(Main$stage)Info: Ploidy: $ploidy\nAA: $test_aa_expressed\nBB: $test_bb_expressed\nDD: $test_dd_expressed\nPloid string: $inferred_genoploidy\nInferredFromGenome: $genomeinfer_from_reads\n";
	}
	if (1) {### For test ###
		print "##### Genome summary #####\nTotal mapped reads: $totalmappedreads\nTagged: $total_tagged\nNotag: $total_notag\nShared reads: $num_shared_reads\nExcluded reads: $num_excluded_reads\n\n";
		
		map {print "Genome: $_\tCount: $genomecount{$_}\n"} (sort (keys %genomecount)); 
		print "\n";
		map {print "Chrom: $_\tCount: $chromcount{$_}\n"} (sort (keys %chromcount));
		print "##### Genome summary #####\n";
		print "(Main$stage)Info: Ploidy: $ploidy\nAA: $test_aa_expressed\nBB: $test_bb_expressed\nDD: $test_dd_expressed\nPloid string: $inferred_genoploidy\n\nInferredFromGenome: $genomeinfer_from_reads";
	}
#	if ($ploidy!=1 or $corrected_num_allelic_sites!=0 or scalar(keys %{$corrected_fragments}) != 0) {
#		print STDERR "(Main$stage)Error: Multiploidy\n";
#		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tMultiploidy:\n";
#		next;
#	}
	
	
	
###Stage8: correct genome alleles using diploid and genome assignemnts
	$stage=8;
	print "\n(Main$stage)Step: Correct Allele\n"; ### For test ###
	my $comprimised=0;
##Format: %CAgenodepth=(chr => pos => allele => depth)
	my $alleledepth={};
	if ($ploidy==1) {
		if ($corrected_num_allelic_sites>0) {
			print "(Main$stage)Info: Comprimised cluster: ploidy=1 but have $corrected_num_allelic_sites allelic sites\n";
			foreach (keys %{$corrected_fragments}) {
				${$readsum}{$_}{'shared'}=1;
			}
			$comprimised=1;
			$corrected_fragments={};
			$corrected_num_allelic_sites=0;
		}
		$test_run_phase=0;
	}
	else {
		(my $test_CAcode, $assignallele, $alleledepth)=CorrectAlleles($inferred_genoploidy, $corrected_fragments, $readsum, \%ancestralalleles);
		unless ($test_CAcode) {
			print STDERR "(Main$stage)Error: Correct Allele\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tCorrectAllele\n";
			next;
		}
	}
#	print "(Main$stage)Test: assignallele\n"; print Dumper $assignallele; print "\n"; ### Test asignallele ###
#	print "(Main$stage)Test: Depth\n"; print Dumper $alleledepth; print "\n"; ### Test depth ###



###Stage9: Fill in genotypes accroding to the corrected alleles
	$stage=9;
	print "\n(Main$stage)Step: Fillin\n"; ### For test ###
	
	my $fixedAllele_3fillin={};
##Format: $fixedGeno_4fillin{chr}->{pos}='0/1/0'; ('A/B/D')
	my $fixedGeno_4fillin={};
	
	if ($ploidy>1) {
		print "(Main$stage)Info: run_phasing before fillin: $test_run_phase\n";
#		print "(Main$stage)Test: Genome assignment ###\n"; print Dumper $assignallele; ### for test ###
		(my $test_fiv, $test_run_phase, my $fixedgeno_hashindex, my $fillinAlleles_hashindex)=&FillInVariations($assignallele,$aabbdd_vcf_obj, $inferred_genoploidy, $corrected_sites, $alleledepth); ###checl Fillin $polymorphic
#		print "(Main$stage)Test: Genome assignment ###\n"; print Dumper $assignallele; ### for test ###
		if ($test_fiv==0) {###Fillin succeeds
			$fixedAllele_3fillin=$fillinAlleles_hashindex;
			$fixedGeno_4fillin=$fixedgeno_hashindex;
		}
		elsif ($test_fiv==1) {
			print STDERR "(Main$stage)Error: FillIn1 SUB failed\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tFillin1\n";
			next;
		}
		elsif ($test_fiv==2) {
			print STDERR "(Main$stage)Error: FillIn2 SUB failed\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tFillin2\n";
			next;
		}
		elsif ($test_fiv==1) {
			print STDERR "(Main$stage)Error: FillIn3 SUB failed\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tFillin3\n";
			next;
		}
		else {
			print STDERR "(Main$stage)Error: FillIn SUB failed\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tFillin4\n";
			next;
		}
#		print "(Main$stage)Test: ### Genome assigned ###\n"; print Dumper $assignallele; print "\n"; ### For test %{$assignallele}###
#		print "(Main$stage)Test: ### Genoypes Fillin ###\n"; print Dumper $fixedAllele_3fillin; print "\n"; ### For test fixedAllele_3fillin###
		$final_geno=$fixedAllele_3fillin;
		print "(Main$stage)Info: run_phasing after fillin: $test_run_phase\n"; ### For test ###
	}
	
	if (0) {### Test $final_geno ###
		foreach my $chrom (keys %{$final_geno}) {
			print "(Main$stage)Test: Chrom: $chrom\n";
			foreach my $pos (keys %{${$assignallele}{$chrom}}) {
				print "(Main$stage)Test: Chrom: $chrom Pos: $pos\n";
			}
		}
	}
	if (0) {
		foreach my $chrom (keys %{$aabbdd_vcf_obj}) {
			print "(Main$stage)Test: AlleleTesting: chromosome $chrom\n" if ($debug); ### For test ###
			foreach my $posit (sort {$a <=> $b } keys %{${$aabbdd_vcf_obj}{$chrom}}) {
				print "(Main$stage)Test: AlleleTesting: positions $posit\n" if ($debug);### For test ###
				next unless (exists ${$corrected_sites}{$chrom} and exists ${$corrected_sites}{$chrom}{$posit});
				my $printline=$cluster_num."\t".$chrom."\t".$posit."\t".${${${$aabbdd_vcf_obj}{$chrom}}{$posit}}[2]."\t";
				#%ancestralalleles
				$printline.="Diploid: \t";
				if (exists $ancestralalleles{$chrom} and exists $ancestralalleles{$chrom}{$posit}) {
					$printline.=(exists $ancestralalleles{$chrom}{$posit}{'A'}) ? $ancestralalleles{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists $ancestralalleles{$chrom}{$posit}{'B'}) ? $ancestralalleles{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists $ancestralalleles{$chrom}{$posit}{'D'}) ? $ancestralalleles{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$assignallele}
				$printline.="Genome: \t";
				if (exists ${$assignallele}{$chrom} and exists ${$assignallele}{$chrom}{$posit}) {
					$printline.=(exists ${$assignallele}{$chrom}{$posit}{'A'}) ? ${$assignallele}{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$assignallele}{$chrom}{$posit}{'B'}) ? ${$assignallele}{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$assignallele}{$chrom}{$posit}{'D'}) ? ${$assignallele}{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$fixedAllele_3fillin}
				$printline.="Fillin: \t";
				if (exists ${$fixedAllele_3fillin}{$chrom} and exists ${$fixedAllele_3fillin}{$chrom}{$posit}) {
					$printline.=(exists ${$fixedAllele_3fillin}{$chrom}{$posit}{'A'}) ? ${$fixedAllele_3fillin}{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$fixedAllele_3fillin}{$chrom}{$posit}{'B'}) ? ${$fixedAllele_3fillin}{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$fixedAllele_3fillin}{$chrom}{$posit}{'D'}) ? ${$fixedAllele_3fillin}{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$final_geno}
				$printline.="Final: \t";
				if (exists ${$final_geno}{$chrom} and exists ${$final_geno}{$chrom}{$posit}) {
					$printline.=(exists ${$final_geno}{$chrom}{$posit}{'A'}) ? ${$final_geno}{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$final_geno}{$chrom}{$posit}{'B'}) ? ${$final_geno}{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$final_geno}{$chrom}{$posit}{'D'}) ? ${$final_geno}{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$fixedGeno_4fillin}
				$printline.="Geno: \t";
				if (exists ${$fixedGeno_4fillin}{$chrom}) {
					my $fillin_geno=(exists ${$fixedGeno_4fillin}{$chrom}{$posit}) ? ${$fixedGeno_4fillin}{$chrom}{$posit} : 'N';
					$printline.=$fillin_geno."\t";
				}
				else {
					$printline.='N'."\t";
				}
				#$test_run_phase
				$printline.=$test_run_phase."\n";
				print ALLELELOG $printline;
			}
		}
	}


##COMMENT: Main($stage) SNP phasing
	$stage=10;
	print "\n(Main$stage)Step: phasing\n"; ### For test ###
	if ($test_run_phase) {
		my $hapcompassvcf="$RunDir/Clust$cluster_num/AABBDD.$cluster_num.hapcompass.vcf";
		my $phasedsolution="$RunDir/Clust$cluster_num/AABBDD.$cluster_num.phased";
		unless (HapcompassVcf($cluster_orivcf, $hapcompassvcf, $fixedGeno_4fillin, $corrected_sites, $vcfreferencer)) {
			print STDERR "(Main$stage)Error: prepare HapCompass VCF error\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tHapcompassVcf\n";
			next;
		}
		unless (RunHapCompass(0, $fragmentfile, $hapcompassvcf, $ploidy, "$fragmentfile.hapcompass", $phasedsolution, ' ', $path_java, $path_hapcompassjar)) {
			print STDERR "(Main$stage)Error: HapCompass running error\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tHapCompass\n";
			next;
		}
		my $phased_geno=dclone($fixedAllele_3fillin);
		my $excluded_chroms={};
		unless (ReadHcOut($phasedsolution, $inferred_genoploidy, $phased_geno, $excluded_chroms)) {
			print STDERR "(Main$stage)Error: HapCompass read hapcompass error\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tHapCompassOut\n";
			next;
		}
		$final_geno=$phased_geno;
		if (scalar(keys %{$excluded_chroms}) > 0) {
			my @temp_excluded=keys %{$excluded_chroms};
			print CLUSTERLOG $cluster_num."\t?\t1\t$stage\tExcluded: @temp_excluded\n";
			foreach (@temp_excluded) {
				delete ${$corrected_sites}{$_} if (exists ${$corrected_sites}{$_});
				foreach (keys %{$corrected_sites}) {
					$corrected_num_allelic_sites+=scalar(keys %{${$corrected_sites}{$_}});
				}
#				print "(Main$stage)Test: $corrected_num_allelic_sites\n";
			}
			foreach my $readid (keys %{$readsum}) {
				my $test_excluded=1;
				foreach my $excludchr (@temp_excluded) {
					unless (exists ${$readsum}{$readid}{'ref'} and exists ${$readsum}{$readid}{'ref'}{$excludchr}) {
						$test_excluded=0;
					}
				}
				if ($test_excluded==1) {
					${$readsum}{$readid}{'excluded'}++;
					delete ${$corrected_fragments}{$readid} if (exists ${$corrected_fragments}{$readid});
				}
				if (exists ${$corrected_fragments}{$readid}) {
					foreach (@temp_excluded) {
						if (exists ${$corrected_fragments}{$readid}{$_}) {
							delete ${$corrected_fragments}{$readid}{$_};
							unless (scalar(keys %{${$corrected_fragments}{$readid}})>0) {
								delete ${$corrected_fragments}{$readid};
							}
						}
					}
				}
			}
			
		}
	}
	else {
		print "(Main$stage)Info: no need to run HapCompass\n";
		print GENOLOG "(Main$stage)Info: no need to run HapCompass\n";
	}
	if (1) {
		foreach my $chrom (keys %{$aabbdd_vcf_obj}) {
			print "(Main$stage)Test: AlleleTesting: chromosome $chrom\n" if ($debug); ### For test ###
			foreach my $posit (sort {$a <=> $b } keys %{${$aabbdd_vcf_obj}{$chrom}}) {
				print "(Main$stage)Test: AlleleTesting: positions $posit\n" if ($debug);### For test ###
				next unless (exists ${$corrected_sites}{$chrom} and exists ${$corrected_sites}{$chrom}{$posit});
				my $printline=$cluster_num."\t".$chrom."\t".$posit."\t".${${${$aabbdd_vcf_obj}{$chrom}}{$posit}}[2]."\t";
				#%ancestralalleles
				$printline.="Diploid: \t";
				if (exists $ancestralalleles{$chrom} and exists ${$ancestralalleles{$chrom}}{$posit}) {
					$printline.=(exists $ancestralalleles{$chrom}{$posit}{'A'} and $ancestralalleles{$chrom}{$posit}{'A'}=~/^\d+$/) ? $ancestralalleles{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists $ancestralalleles{$chrom}{$posit}{'B'} and $ancestralalleles{$chrom}{$posit}{'B'}=~/^\d+$/) ? $ancestralalleles{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists $ancestralalleles{$chrom}{$posit}{'D'} and $ancestralalleles{$chrom}{$posit}{'D'}=~/^\d+$/) ? $ancestralalleles{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$assignallele}
				$printline.="Genome: \t";
				if (exists ${$assignallele}{$chrom} and exists ${${$assignallele}{$chrom}}{$posit}) {
					$printline.=(exists ${$assignallele}{$chrom}{$posit}{'A'}) ? ${$assignallele}{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$assignallele}{$chrom}{$posit}{'B'}) ? ${$assignallele}{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$assignallele}{$chrom}{$posit}{'D'}) ? ${$assignallele}{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$fixedAllele_3fillin}
				$printline.="Fillin: \t";
				if (exists ${$fixedAllele_3fillin}{$chrom} and exists ${${$fixedAllele_3fillin}{$chrom}}{$posit}) {
					$printline.=(exists ${$fixedAllele_3fillin}{$chrom}{$posit}{'A'}) ? ${$fixedAllele_3fillin}{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$fixedAllele_3fillin}{$chrom}{$posit}{'B'}) ? ${$fixedAllele_3fillin}{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$fixedAllele_3fillin}{$chrom}{$posit}{'D'}) ? ${$fixedAllele_3fillin}{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#%{$fixedGeno_4fillin}
				$printline.="Geno: \t";
				if (exists ${$fixedGeno_4fillin}{$chrom}) {
					my $fillin_geno=(exists ${$fixedGeno_4fillin}{$chrom}{$posit}) ? ${$fixedGeno_4fillin}{$chrom}{$posit} : 'N';
					$printline.=$fillin_geno."\t";
				}
				else {
					$printline.='N'."\t";
				}
				#%{$final_geno}
				$printline.="Final: \t";
				if (exists ${$final_geno}{$chrom} and exists ${$final_geno}{$chrom}{$posit}) {
					$printline.=(exists ${$final_geno}{$chrom}{$posit}{'A'}) ? ${$final_geno}{$chrom}{$posit}{'A'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$final_geno}{$chrom}{$posit}{'B'}) ? ${$final_geno}{$chrom}{$posit}{'B'} : 'N';
					$printline.="\t";
					$printline.=(exists ${$final_geno}{$chrom}{$posit}{'D'}) ? ${$final_geno}{$chrom}{$posit}{'D'} : 'N';
					$printline.="\t";
				}
				else {
					$printline.='N'."\t".'N'."\t".'N'."\t";
				}
				#$test_run_phase
				$printline.=$test_run_phase."\n";
				print ALLELELOG $printline;
			}
		}
	}

	if (0) {### Test %{$final_geno} ###
		foreach my $chrom (keys %{$final_geno}) {
			print "(Main$stage)Test: Chrom: $chrom\n";
			foreach my $pos (keys %{${$assignallele}{$chrom}}) {
				print "(Main$stage)Test: Chrom: $chrom Pos: $pos\n";
				if (defined ${${$assignallele}{$chrom}}{$pos}){
					print "NOT defined\n";
				}
				else {
					print "Defined\n";
				}
			}
		}
	}



###Stage11: Get phlymorphic readids
	$stage=11;
	print "\n(Main$stage)Step: Get allelic reads\n"; ### For test ###
##Format: %{$final_geno}=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
##Format: %{$corrected_fragments}=($readid => chrom => $partstart => $pos => $allele => $qual);
#	%polymorphismreadids=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++),
	if ($corrected_num_allelic_sites>0) {
#		print "(Main$stage)Test: \$corrected_fragments\n"; print Dumper $corrected_fragments; print "\n"; ### For test ###
#		print "(Main$stage)Test: \$final_geno\n"; print Dumper $final_geno; print "\n"; ### For test ###
		(my $test_GPRcode, $polymorphismreadids)=&GetPolymorphicReads($corrected_fragments, $final_geno);
	}
	if (1) {
		print "(Main$stage)Test: Total mapped reads: ".$totalmappedreads."\n";
		print "(Main$stage)Test: Total allelic reads: ".scalar(keys %{$polymorphismreadids})."\n";
		print GENOLOG "(Main$stage)Test: Total mapped reads: ".$totalmappedreads."\n";
		print GENOLOG "(Main$stage)Test: Total allelic reads: ".scalar(keys %{$polymorphismreadids})."\n";
	}



###Stage12: Bam2fastq2
	$stage=12;
	print "\n(Main$stage)Step: bam2fastq\n"; ### For test ###
	my ($test_bam2fastq, $fastqfiles, $rgfpkms, $clusterfpkm)=&Bam2FastQ2($cluter_oribam, $inferred_genoploidy, $readsum, $polymorphismreadids, "$RunDir/Clust$cluster_num/reads", 1, 0, 1, $path_samtools);
	if (! $test_bam2fastq) {
		print STDERR "(Main$stage)Error: Bam2FastQ2\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tBam2FastQ2\n";
		next;
	}
##Format: %{$rgfpkms}=(chr => ('chr' => FPKM
#						'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#						'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#											RG2 => ...)
	if (1) {
#		print Dumper $rgfpkms; ### For test ###
		foreach my $idvchrom (keys %{$rgfpkms}) {
			if (exists ${$rgfpkms}{$idvchrom}{'chr'}) {
				print REFERFPKM "Clust$cluster_num\tChr\t".${$rgfpkms}{$idvchrom}{'chr'}."\n";
			}
			if (exists ${$rgfpkms}{$idvchrom}{'allelic'}) {
				foreach (sort (keys %{${$rgfpkms}{$idvchrom}{'allelic'}})) {
					print REFERFPKM "Clust$cluster_num\tAllelic\t$_\t".${$rgfpkms}{$idvchrom}{'allelic'}{$_}."\n";
				}
			}
			if (exists ${$rgfpkms}{$idvchrom}{'read_group'}) {
				foreach my $ind_rg (sort (keys %{${$rgfpkms}{$idvchrom}{'read_group'}})) {
					print REFERFPKM "Clust$cluster_num\tReadGroup\t$ind_rg";
					foreach my $indg (sort (keys %{${$rgfpkms}{$idvchrom}{'read_group'}{$ind_rg}})) {
						print REFERFPKM "\t$indg\t".${$rgfpkms}{$idvchrom}{'read_group'}{$ind_rg}{$indg};
					}
					print REFERFPKM "\n";
				}
			}
		}
	}
##Format $clusterfpkm=('genomesize' => size, 
#							'cluster' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#											RG2 => ...)
	if (1) {
		if (exists ${$clusterfpkm}{'genomesize'}){
			print CLUSTFPKM "Clust$cluster_num\tCLUST\tGenomeSize:".${$clusterfpkm}{'genomesize'}."\n";
		}
		else {
			print CLUSTFPKM "Clust$cluster_num\tCLUST\tGenomeSize:N\n";
		}
		if (exists ${$clusterfpkm}{'cluster'}){
			print CLUSTFPKM "Clust$cluster_num\tClustGlobal\t".${$clusterfpkm}{'cluster'}."\n";
		}
		else {
			print CLUSTFPKM "Clust$cluster_num\tClustGlobal\tN\n";
		}
		if (exists ${$clusterfpkm}{'allelic'}) {
			foreach (sort (keys %{${$clusterfpkm}{'allelic'}})) {
				print CLUSTFPKM "Clust$cluster_num\tAllelic\t$_\t".${$clusterfpkm}{'allelic'}{$_}."\n";
			}
		}
		if (exists ${$clusterfpkm}{'read_group'}) {
			foreach my $ind_rg (sort (keys %{${$clusterfpkm}{'read_group'}})) {
				print CLUSTFPKM "Clust$cluster_num\tReadGroup\t$ind_rg";
				foreach my $indg (sort (keys %{${$clusterfpkm}{'read_group'}{$ind_rg}})) {
					print CLUSTFPKM "\t$indg\t".${$clusterfpkm}{'read_group'}{$ind_rg}{$indg};
				}
				print CLUSTFPKM "\n";
			}
		}
	}


### (Main13) Trinity
	$stage=13;
	print "\n(Main$stage)Step: Trinity\n"; ### For test ###
	my %fasta_assembled=();
	my $test_trinity=0;
	foreach (keys %{$fastqfiles}) {
		unless (defined ${$fastqfiles}{$_} and -s ${$fastqfiles}{$_}) {
			print STDERR "(Main$stage)Error: Geno: $_ fastq file empty or not exists: ${$fastqfiles}{$_}\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tTrinity-$_\n";
			$test_trinity=1;
			next;
		}
		unless (RunFqTrinity(${$fastqfiles}{$_}, "$RunDir/Clust$cluster_num/Clust$cluster_num.trinity.$_.fasta", $trinity_addcmd, $path_trinity)) {
			print STDERR "(Main$stage)Error: Trinity-$_\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tTrinity-$_\n";
			$test_trinity=1;
			next;
		}
		$fasta_assembled{$_}="$RunDir/Clust$cluster_num/Clust$cluster_num.trinity.$_.fasta";
	}
	if ($test_trinity==1) {
		print STDERR "(Main$stage)Error: Trinity may partially failed\n";
		next;
	}
	unless (chdir $RunDir) {
		print STDERR "(Main$stage)Error: can not chdir to : $RunDir\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tChDirAfterTrinity\n";
		next;
	}


###Stage14: rename fasta
	$stage=14;
	print "\n(Main$stage)Step: FastaRenamer\n"; ### For test ###
	foreach (keys %fasta_assembled) {
		unless (RenameFasta($fasta_assembled{$_}, "$fasta_assembled{$_}.rename.fa", "Cluster${cluster_num}_${_}_", 9, "Geno=$inferred_genoploidy Chrom=$inferred_chrom")) {
			print STDERR "(Main$stage)Error: FastaRenamer: $_\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tFastaRenamer\n";
			next;
		}
		unless (MoveFile("$fasta_assembled{$_}.rename.fa", "$RunDir/fasta/")) {
			print STDERR "(Main$stage)Error: MoveFile: $_\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tMoveFile\n";
			next;
		}
	}
### merge fasta
=notsure
### (Main14) CD-HIT
	$stage=13;
	my %fasta_cdhit=();
	foreach (keys %fasta_assembled) {
		unless (&CdHitEst($fasta_assembled{$_}, "$fasta_assembled{$_}.cdhit.fasta", $cdhit_addcmd, $path_cdhitest)) {
			print STDERR "(Main$stage)Error: cdhit: $_\n";
			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tCdHitEst\n";
			next;
		}
		else {
			push (@fasta_cdhit, "$_.cdhit.fasta");
			unlink $_ if (-e $_);
		}
	}
=cut
###Stage15: Cleaning
	$stage=15;
	print "\n(Main$stage)Step: cleaning\n"; ### For test ###
	unless (chdir $RunDir) {
		print STDERR "(Main$stage)Error: can not chdir to : $RunDir\n";
		print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tChDirAtLAST\n";
		next;
	}
	if ($cleanlog) {
#		sleep(2);
		unlink glob "$RunDir/Clust$cluster_num/*";
#		unless (rmtree("$RunDir/Clust$cluster_num")) {
#			print STDERR "(Main$stage)Error: DeletePATH: $RunDir/Clust$cluster_num\n";
#			print CLUSTERLOG $cluster_num."\tFail\t1\t$stage\tDeletePATH\n";
#		}
	}



### Stage16: Ending
	$stage=16;
	print CLUSTERLOG $cluster_num."\tSucceed\t0\t$stage\t0\t$genomeinfer_from_reads\t$inferred_genoploidy\t\n";
}
close CLUSTER;
close CLUSTERLOG;
close REFERFPKM;
close CLUSTFPKM;
close ALLELELOG;
close GENOLOG;


#####################################################################
###                         sub functions                         ###
#####################################################################



### extract batch of sequence alignment from Bams, and index
### &GetBam(bam_arr, sequence_arr, output)
### Global:$cluster_num, $path_samtools
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
		if (scalar(@{$EBfiles_bam_index}) == 1) {
			$EBbamout="$EBoutput.multihead.bam";
		}
		else {
			$EBbamout="$EBoutput.$EBi.bam";
		}
		my $EBcmd="$path_samtools view -b -h -F 4 $EBind_bam $EBsam_seq > $EBbamout";
		if (! &exec_cmd_return($EBcmd)) {
			print STDERR "${EBsubinfo}Error: bam extract $EBsam_seq (Cluster: $cluster_num) from $EBind_bam failed\n";
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
			print STDERR "${EBsubinfo}Error: bam merge (Cluster: $cluster_num) from $EBmerge_bams failed\n";
			return 0;
		}
		else {
			if (-s "$EBoutput.multihead.bam") {
				map {unlink($_)} @EBtemp_bams;###delete temporary files ### For test ###
				
			}
			else{
				print STDERR "${EBsubinfo}Error: bam merge (Cluster: $cluster_num) from $EBmerge_bams not found\n";
				return 0;
			}
		}
	}
##COMMENT: cleanheader
	my $EBtest_cleanheader=&SamCleanHeader("$EBoutput.multihead.bam", "$EBoutput.withRG.bam", "$EBoutput.noRG.bam", $EBseq_obj,  $path_samtools);
	if ($EBtest_cleanheader==0) {
		print STDERR "${EBsubinfo}Error: bam header clean (Cluster: $cluster_num) failed\n";
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
			print STDERR "${EBsubinfo}Error: bam headerclean output (Cluster: $cluster_num) not found\n";
			return 0;
		}
	}
	return 1;
}



### Read FPKM for cluster
### &ReadFpkm($filefpkm)
### Global: $minimum_fpkm
### Dependency:
### Note:
sub ReadFpkm {
	my ($RFfpkmin, $RFseqids_arrindex)=@_;
	
	my $RFsubinfo='SUB(Main::ReadFpkm)';
	my @RFtest_expressed=();
	my $RFtotalnum=scalar(@{$RFseqids_arrindex});
	my %RFseqid=();
	map {$RFseqid{$_}++} @{$RFseqids_arrindex};
	my ($RFfpkm1, $RFfpkm2, $RFfpkm3, $RFfpkm4)=(0, 0, 0, 0);
	local *RFFPKMIN;
	close RFFPKMIN if (defined fileno(RFFPKMIN));
	unless (open(RFFPKMIN, "< $RFfpkmin")) {
		print STDERR "${RFsubinfo}Error: can not open FPKM input\n";
		return 0;
	}
	
	my $RFmatchednum=0;
	while (my $RFline=<RFFPKMIN>) {
		chomp $RFline;
		my @RFarr=split(/\t/, $RFline);
		if (exists $RFseqid{$RFarr[0]}) {
			$RFmatchednum++;
			$RFfpkm1=1 if ($RFarr[1] > $minimum_fpkm);
			$RFfpkm2=1 if ($RFarr[2] > $minimum_fpkm);
			$RFfpkm3=1 if ($RFarr[3] > $minimum_fpkm);
			$RFfpkm4=1 if ($RFarr[4] > $minimum_fpkm);
			if ($RFfpkm1==1 and $RFfpkm2==1 and $RFfpkm3==1 and $RFfpkm4==1) {
				last;
			}
			elsif ($RFmatchednum==$RFtotalnum) {
				last;
			}
		}
	}
	close RFFPKMIN;
	@RFtest_expressed=($RFfpkm1, $RFfpkm2, $RFfpkm3, $RFfpkm4);
	
	return (1, \@RFtest_expressed);
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



### Group SNPs based on read length
### &GroupReads($ReadSam_obj, AABBDD_Vcf_obj)
### Global: $debug, $verbose
### Dependancy:&ReadVariantType, $min_mapq, $bam_genome_tag, $bam_chromo_tag
### Note:
sub GroupReads {
	my ($GRaabbdd_samobj, $GRaabbdd_vcfobj)=@_;

	my $GRsubinfo="SUB(GroupReads)";
#Format: %GRreadid_by_allele=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
#	my %GRreadid_by_allele=();
#Format: %GRreadfragments=(readid => (chrom => (pos => (geno1 => quality, geno2 => quality)));
	my %GRreadfragments=();
#Format: %GRreadid2genome=(readID=> ('A' =>1, 'B' =>1, 'C'=1))
	my %GRreadid2genome=();
#Format: %GRreadidsum=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
#								'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#								'shared' => ++/0,
#								'mapped' => 1/0,
#								'ref' => Reference=> ++
#								'range' => (Reference => ((start, end), (start2-end2)...)
#								'excluded'=> 0/++)
#							)
#					);
	my %GRreadidsum=();
#Format: %GRallele_assign=($chr => ($pos => ('A' => (allele1 =>++, allele2 =>++), 
#											'B' => (allele1 =>++, allele2 =>++), 
#											'D' => (allele1 =>++, allele2 =>++)))
	my %GRallele_assign=();
	my $GRtotalreads_unmapped=0;
	my $GRtotalreads_lowqual=0;

###Retrieve all reads into share group
	foreach my $GRchrom (@cluster_seqids) {
###COMMENT: retrieve all alignments mapped to this reference sequence
		my @GRchr_alignments=$GRaabbdd_samobj->get_features_by_location(-seq_id => "$GRchrom");
		foreach my $GRind_align (@GRchr_alignments) {
			if ($GRind_align->flag & 0x0004) {# makesure mapped
				$GRtotalreads_unmapped++;
				next;
			}
			unless ($GRind_align->qual >= $min_mapq) {#filter mapping quality
				$GRtotalreads_lowqual++;
				next;
			}
			my $GRthisreadid=$GRind_align->name;
			${$GRreadidsum{$GRthisreadid}}{'mapped'}++;
			${$GRreadidsum{$GRthisreadid}}{'shared'}++;
			$GRreadidsum{$GRthisreadid}{'ref'}{$GRchrom}++;
			my @GRstartend=($GRind_align->start, $GRind_align->end);
			unless (defined $GRstartend[0] and $GRstartend[0]=~/^\d+$/ and defined $GRstartend[1] and $GRstartend[1]=~/^\d+$/) {
				print STDERR "${GRsubinfo}Test: undefined (Chr: $GRchrom) range for read: $GRthisreadid\n";
			}
			else {
				push (@{$GRreadidsum{$GRthisreadid}{'range'}{$GRchrom}}, \@GRstartend);
			}
			${$GRreadidsum{$GRthisreadid}}{'excluded'}=0;
#			next if (exists $GRreadidsum{$GRthisreadid});###Reduce mem
			my $GRtag_chrom=$GRind_align->get_tag_values("$bam_chromo_tag");
			while ($GRtag_chrom=~/(\d[ABDLS]{1,2})/g) {
				${${$GRreadidsum{$GRthisreadid}}{'chr'}}{$1}++;
			}
			my $GRgenome_tagvalue=$GRind_align->get_tag_values("$bam_genome_tag");
			if ($GRgenome_tagvalue=~m/A/i) {
				${$GRreadid2genome{$GRthisreadid}}{'A'}=1;
				${${$GRreadidsum{$GRthisreadid}}{'gen'}}{'A'}++;
			}
			if ($GRgenome_tagvalue=~m/B/i) {
				${$GRreadid2genome{$GRthisreadid}}{'B'}=1;
				${${$GRreadidsum{$GRthisreadid}}{'gen'}}{'B'}++;
			}
			if ($GRgenome_tagvalue=~m/D/i) {
				${$GRreadid2genome{$GRthisreadid}}{'D'}=1;
				${${$GRreadidsum{$GRthisreadid}}{'gen'}}{'D'}++;
			}
		}
	}

### read allele for each reads
	foreach my $GRchrom (keys %{$GRaabbdd_vcfobj}) {
		my @GRpositions=sort {$a<=>$b} (keys %{${$GRaabbdd_vcfobj}{$GRchrom}});
		print "${GRsubinfo}Test: Number of variations on $GRchrom: ", scalar(@GRpositions), "\n" if ($debug or $verbose);
		foreach my $GRind_pos (sort {$a <=> $b} @GRpositions) {
			print "${GRsubinfo}Test:\t".$GRchrom, "\t", "Pos: $GRind_pos\tRef:${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[0]\tVar:${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[1]\tGen:${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[2]", "\n" if ($debug or $verbose);
#			my @GRallelearr=split(/\//, ${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[2]);
#			my %GRallelehash=();
#			map {$GRallelehash{$_}++} @GRallelearr;
#			my $GRnumallele=scalar(keys %GRallelehash);
#			if ($GRnumallele<1 or $GRnumallele>3) {
#				print STDERR "${GRsubinfo}Error: Chr:Pos $GRchrom:$GRind_pos: Number of allele <1 or >3\n";
#				return 1;
#			}
			my $GRend_pos=$GRind_pos+length(${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[0])-1;
			my @GRchr_alignments=();
			@GRchr_alignments= $GRaabbdd_samobj->get_features_by_location(-seq_id => "$GRchrom", -start => "$GRind_pos", -end => "$GRend_pos");
			
			foreach my $GRind_align (@GRchr_alignments) {
				my $GRthisreadid=$GRind_align->name;
				next if ($GRind_align->flag & 0x0004);
				${$GRreadidsum{$GRthisreadid}}{'shared'}=0;# if ($GRnumallele>1);
				next unless ($GRind_align->qual >= $min_mapq);
				my $GRstrand='U';
				if ($GRind_align->flag & 0x0040) {
					$GRstrand=1;
				}
				elsif ($GRind_align->flag & 0x0080) {
					$GRstrand=2;
				}
				#$GRalelecode: this read contains which allele genotype (0(ref)/1/2/3/.)
				my ($GRalelecode, $GRcapture_qual)=&ReadVariantType($GRind_pos, ${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[0], ${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[1], ${${${$GRaabbdd_vcfobj}{$GRchrom}}{$GRind_pos}}[2], $GRind_align, 2);
				unless ($GRalelecode =~ /^\d+$/) {
#				unless ($GRalelecode =~ /^\d+$/ and exists $GRallelehash{$GRalelecode}) {
					$GRalelecode='N';
#					$GRcapture_qual='#';
#					next;
				}
#				if ($GRnumallele>1) {
###For HapCompass fragments  file
					if (! defined $GRcapture_qual or length($GRcapture_qual) !=1) {
						$GRcapture_qual='I'; ### The second highest phred+33 score
					}
					elsif (ord($GRcapture_qual)<33 or ord ($GRcapture_qual) >74) {
						$GRcapture_qual='I'; ### The second highest phred+33 score
					}
					if (exists $GRreadfragments{$GRthisreadid} and exists $GRreadfragments{$GRthisreadid}{$GRchrom} and exists $GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos} and exists $GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos}{$GRstrand} and exists $GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos}{$GRstrand}{$GRalelecode}) {
						if (ord($GRcapture_qual) > ord($GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos}{$GRstrand}{$GRalelecode})) {
							$GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos}{$GRstrand}{$GRalelecode}=$GRcapture_qual;
						}
					}
					else {
						$GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos}{$GRstrand}{$GRalelecode}=$GRcapture_qual;
					}
#					push (@{${${$GRreadid_by_allele{$GRchrom}}{$GRind_pos}}{$GRalelecode}}, $GRthisreadid);
#				}
#				else {
#					$GRreadfragments{$GRthisreadid}{$GRchrom}{$GRind_pos}{$GRstrand}{$GRalelecode}=$GRcapture_qual;
#				}
			}
		}
	}

#Merge parts that can overlap
	foreach my $GRreadid (keys %GRreadidsum) {
		next unless (exists $GRreadidsum{$GRreadid}{'range'});
		foreach my $GRchrom (keys %{$GRreadidsum{$GRreadid}{'range'}}) {
			my @GRarr=@{$GRreadidsum{$GRreadid}{'range'}{$GRchrom}};
#			print "${GRsubinfo}Test: range before merge\n"; ### For test ###
#			print Dumper \@GRarr; ### For test ###
			
			my ($GRtestmerge, $GRmergearr)=&MergeRange(\@GRarr);
			unless (defined $GRtestmerge and $GRtestmerge>0) {
				print STDERR "${GRsubinfo}Error: mergerage failed\n";
				return 1;
			}
#			print "${GRsubinfo}Test: range after merge\n";### For test ###
#			print Dumper $GRmergearr;### For test ###
			my %GRmergehash=();
			map {$GRmergehash{$_->[0]}=$_->[1]} @{$GRmergearr};
#			print "${GRsubinfo}Test: hash after merge\n";### For test ###
#			print Dumper \%GRmergehash;### For test ###
			$GRreadidsum{$GRreadid}{'merge'}{$GRchrom}=\%GRmergehash;
		}
	}

	if (1) {
		print "${GRsubinfo}Test: total reads: ", scalar(keys %GRreadidsum), "\n";
		print "${GRsubinfo}Test: total reads unmapped: $GRtotalreads_unmapped\n";
		print "${GRsubinfo}Test: total reads low-quality: $GRtotalreads_lowqual\n";
		print "${GRsubinfo}Test: total allelic reads: ", scalar(keys %GRreadfragments), "\n";
	}
	return (0, \%GRreadfragments, \%GRreadidsum);
}



### intersect Bioranges into intervels
### &MergeRange(\@arr)
### Global:
### Dependency:
### Note:
### Example: [ [1, 14], [2, 10], [40, 60], [30,70], [71, 90]] => [[1,14], [40, 90]]
sub MergeRange {
	my $MRrangeori=shift;
	
	my $MRsubinfo="SUB(MergeRange)";
	my %MRrangehash=();
	my %MRmergehash=();
	my @MRretarr=();
	
	if (scalar(@{$MRrangeori})==0) {
		print STDERR $MRsubinfo, "Error: empty range\n";
		return 0;
	}
	elsif (scalar(@{$MRrangeori})==1) {
#		print $MRsubinfo, "Info: 1 range\n"; ### For test ###
		return (1, $MRrangeori);
	}
	
	foreach my $MRindarr (@{$MRrangeori}) {
		my ($MRstart, $MRend)=@{$MRindarr};
		unless (defined $MRstart and $MRstart=~/^\d+$/ and defined $MRend and $MRend=~/^\d+$/) {
			print STDERR $MRsubinfo, "Error: non number range\n";
			return 0;
		}
		if (exists $MRrangehash{$MRstart}) {
			$MRrangehash{$MRstart}=$MRend if ($MRend> $MRrangehash{$MRstart});
		}
		else {
			$MRrangehash{$MRstart}=$MRend;
		}
	}
#	print $MRsubinfo, "Info: merge\n"; ### For test ###
	my $MRtest_first=0;
	my ($MRmin, $MRmax);
	my $MRprinted=0;
	foreach my $MRindstart (sort {$a <=> $b} keys %MRrangehash) {
		$MRtest_first++;
		if ($MRtest_first==1) {
			$MRmin=$MRindstart;
			$MRmax=$MRrangehash{$MRindstart};
			next;
		}
		else {
			if ($MRindstart>$MRmin and $MRindstart <=($MRmax+1)) {
				if ($MRmax<$MRrangehash{$MRindstart}) {
					$MRmax=$MRrangehash{$MRindstart};
				}
			}
			else {
				my @MRarr2=($MRmin, $MRmax);
				push (@MRretarr, \@MRarr2);
				$MRmin=$MRindstart;
				$MRmax=$MRrangehash{$MRindstart};
			}
		}
		if ($MRtest_first==scalar(keys %MRrangehash)) {
			my @MRarr3=($MRmin, $MRmax);
			push (@MRretarr, \@MRarr3);
		}
	}
#	print Dumper \@MRretarr; ### For test ###
	return (1, \@MRretarr);
}



### FillInVariations
###&FillInVariations()
###Global: $debug
###Dependency:&ExtractAllele
###Note
### Return (0=successful, 1=ExtraAlleles, 2=InsufficientAlleles, 3=Unknown)
sub FillInVariations {
	my ($FIVfixallele_hashindex, $FIVvcf_obj, $FIVgenoploidy, $FIVpolymorphicstes, $FIVdepth)=@_;
##Format: %{$FIVdepth}=(chr => pos => allele => depth)

	my $FIVsubinfo ='SUB(FillInVariations)';
##COMMENT: Test SUB(FillInVariations) successful or not
	my $FIVneed_runcompass=0;###Test need to run next step compass (if there is any unfixed allele) or not (if all alleles are fixed);
#	%{$FIVpolymorphicstes}{chr}{$pos}
##Format: %{$FIVfixallele_hashindex}=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my $FIVfillinalleles=dclone($FIVfixallele_hashindex);
##Format: ${$FIVvcf_obj}{chr}->{pos}=(ref, var, gen/gen2, ...);
##Format: $FIVfixed_geno{chr}->{pos}='0/1/0'; ('A/B/D')
	my %FIVfixed_geno=();
##COMMENT: decide which subgenome to output
	my %FIVgenohash=();
	if ($FIVgenoploidy=~/A/) {
		$FIVgenohash{'A'}++;
#		print $FIVsubinfo, "AA\n"; ### For test ###
	}
	if ($FIVgenoploidy=~/B/) {
		$FIVgenohash{'B'}++;
#		print $FIVsubinfo, "BB\n"; ### For test ###
	}
	if ($FIVgenoploidy=~/D/) {
		$FIVgenohash{'D'}++;
#		print $FIVsubinfo, "DD\n"; ### For test ###
	}
	
	my $FIVploidy=scalar(keys %FIVgenohash);
	if ($FIVploidy<2 or $FIVploidy>3) {
		print $FIVsubinfo, "Error: Total ploidy: $FIVploidy\n";### For test ###
		return 1;
	}
#	print $FIVsubinfo, "Test: Total ploidy: $FIVploidy\n"; ### For test ###
##COMMENT: Fill in unknown alleles
	foreach my $FIVchrom (sort keys %{$FIVpolymorphicstes}) {
		foreach my $FIVpos (sort keys %{${$FIVpolymorphicstes}{$FIVchrom}}) {
			my @FIVallelecode=keys %{${$FIVpolymorphicstes}{$FIVchrom}{$FIVpos}};
### Get read depth after clean
			my %FIValleledep=();
			foreach (@FIVallelecode) {
				unless (exists ${$FIVdepth}{$FIVchrom} and exists ${$FIVdepth}{$FIVchrom}{$FIVpos} and exists ${$FIVdepth}{$FIVchrom}{$FIVpos}{$_} and ${$FIVdepth}{$FIVchrom}{$FIVpos}{$_}>0) {
					print STDERR $FIVsubinfo, "Error: depth==0 at Chr:Pos:Allele $FIVchrom:$FIVpos:$_\n";
					return 1;
				}
				$FIValleledep{$_}=${$FIVdepth}{$FIVchrom}{$FIVpos}{$_};
			}
			unless (scalar(keys %FIValleledep)>1) {
				print STDERR $FIVsubinfo, "Error: allele <1 at Chr:Pos $FIVchrom:$FIVpos\n";
				return 1;
			}
###rephase geno hash
			my @FIVvcfgeno=();
			if (scalar(@FIVallelecode)==3) {
				@FIVvcfgeno=@FIVallelecode;
			}
			elsif (scalar(@FIVallelecode)==2) {
				@FIVvcfgeno=@FIVallelecode;
				my $FIVmaxdepth=0;
				my $maxcount=0;
				my $FIVmaxallelecoce;
				foreach (@FIVallelecode) {
					if ($FIValleledep{$_}>$FIVmaxdepth) {
						$FIVmaxdepth=$FIValleledep{$_};
						$maxcount=1;
						$FIVmaxallelecoce=$_;
					}
					elsif ($FIValleledep{$_}==$FIVmaxdepth) {
						$maxcount++;
					}
				}
				if (defined $FIVmaxallelecoce and $FIVmaxallelecoce=~/^\d+$/) {
					push (@FIVvcfgeno, $FIVmaxallelecoce);
				}
				else {
					if (exists ${$FIVvcf_obj}{$FIVchrom} and exists ${$FIVvcf_obj}{$FIVchrom}{$FIVpos} and defined ${${$FIVvcf_obj}{$FIVchrom}{$FIVpos}}[2] and ${${$FIVvcf_obj}{$FIVchrom}{$FIVpos}}[2]=~/^\d+\/\d+\/\d+$/) {
						my @FIVvcf_thisgeno=split(/\//, ${${$FIVvcf_obj}{$FIVchrom}{$FIVpos}}[2]);
						unless (scalar(@FIVvcf_thisgeno)==3) {
							print STDERR $FIVsubinfo, "Error: genotype not triple in VCF at Chr:Pos:Geno $FIVchrom:$FIVpos:${${$FIVvcf_obj}{$FIVchrom}{$FIVpos}}[2]\n";
							return 1;
						}
						my %FIVvcf_genohash=();
						map {$FIVvcf_genohash{$_}++} @FIVvcf_thisgeno;
						foreach (keys %FIVvcf_genohash) {
							if ($FIVvcf_genohash{$_}>=2) {
								push (@FIVvcfgeno, $_);
							}
						}
					}
					else {
						print STDERR $FIVsubinfo, "Error: error to retrieve genotype in VCF at Chr:Pos $FIVchrom:$FIVpos\n";
						return 1;
					}
				}
			}
			else {
				print STDERR $FIVsubinfo, "Error: unknown fillin Ploidy: $FIVploidy AlleleNum: ".scalar(@FIVallelecode)." Allele @FIVallelecode at Chr:Pos $FIVchrom:$FIVpos\n";
				return 1;
			}
			unless (scalar(@FIVvcfgeno)==3) {
				print STDERR $FIVsubinfo, "Error: ploidy !=3 at Chr:Pos $FIVchrom:$FIVpos, @FIVvcfgeno\n";
				return 1;
			}
###Format: %vcfgenohash=(0 => ++, 1=> ++, 2 => ++);
			my %FIVvcfgenhash=();
			map {$FIVvcfgenhash{$_}++} @FIVvcfgeno;
			
			my $FIVnum_allele_notsure=0;
			my ($FIVaa_alleles,$FIVbb_alleles,$FIVdd_alleles)=('?', '?', '?');
			if (exists ${$FIVfillinalleles}{$FIVchrom} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}) {
				$FIVaa_alleles=${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'} if (exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'});
				$FIVbb_alleles=${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'} if (exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'});
				$FIVdd_alleles=${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'} if (exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'});
			}
			my ($FIVmissing_aa, $FIVmissing_bb, $FIVmissing_dd)=(0, 0, 0);
			my %FIVallele_taken=();
			if (exists $FIVgenohash{'A'}) {
				if ($FIVaa_alleles =~/^\d+$/ and exists $FIVvcfgenhash{$FIVaa_alleles}) {
					$FIVvcfgenhash{$FIVaa_alleles}--;
					$FIVallele_taken{$FIVaa_alleles}++;
				}
				else {
					$FIVnum_allele_notsure++;
					$FIVmissing_aa=1;
				}
			}
			if (exists $FIVgenohash{'B'}) {
				if ($FIVbb_alleles =~/^\d+$/ and exists $FIVvcfgenhash{$FIVbb_alleles}) {
					$FIVvcfgenhash{$FIVbb_alleles}--;
					$FIVallele_taken{$FIVbb_alleles}++;
				}
				else {
					$FIVnum_allele_notsure++;
					$FIVmissing_bb=1;
				}
			}
			if (exists $FIVgenohash{'D'}) {
				if ($FIVdd_alleles =~/^\d+$/ and exists $FIVvcfgenhash{$FIVdd_alleles}) {
					$FIVvcfgenhash{$FIVdd_alleles}--;
					$FIVallele_taken{$FIVdd_alleles}++;
				}
				else {
					$FIVnum_allele_notsure++;
					$FIVmissing_dd=1;
				}
			}
			my @FIVallele_left=();
			foreach (keys %FIVvcfgenhash) {
				push (@FIVallele_left, $_) if ($FIVvcfgenhash{$_}>0);
			}
			if (scalar(@FIVallele_left)<1 and $FIVnum_allele_notsure>0) {
				print STDERR $FIVsubinfo, "Error: no left allele at $FIVchrom:$FIVpos\n";###For test###
				return 2;
			}
			my @FIVarr_geno_fix=();
			if ($FIVnum_allele_notsure==0) {
				push (@FIVarr_geno_fix, ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'}) if (exists $FIVgenohash{'A'});
				push (@FIVarr_geno_fix, ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'}) if (exists $FIVgenohash{'B'});
				push (@FIVarr_geno_fix, ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'}) if (exists $FIVgenohash{'D'});
			}
			elsif ($FIVnum_allele_notsure==1) {
				if (scalar(@FIVallele_left)!=1) {
					my @FIVallele_left_new=();
					foreach (@FIVallele_left) {
						push (@FIVallele_left_new, $_) unless (exists $FIVallele_taken{$_});
					}
					if (scalar(@FIVallele_left_new)!=1) {
						print STDERR $FIVsubinfo, "Error: extra at $FIVchrom: $FIVpos\n";###For test###
						print "AllelesLeft: @FIVallele_left_new\nAlleleNum: $FIVnum_allele_notsure\n";###For test###
						return 1;
					}
					else {
						@FIVallele_left=@FIVallele_left_new;
					}
				}
				${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'}=shift @FIVallele_left if ($FIVmissing_aa ==1);
				${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'}=shift @FIVallele_left if ($FIVmissing_bb ==1);
				${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'}=shift @FIVallele_left if ($FIVmissing_dd ==1);
				push (@FIVarr_geno_fix, ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'}) if (exists $FIVgenohash{'A'});
				push (@FIVarr_geno_fix, ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'}) if (exists $FIVgenohash{'B'});
				push (@FIVarr_geno_fix, ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'}) if (exists $FIVgenohash{'D'});
			}
			elsif ($FIVnum_allele_notsure==2) {
				my $uniqueallele=0;
				if (scalar(@FIVallele_left)==1){###Check if only one allele left
					$FIVallele_left[1]=$FIVallele_left[0];
					$uniqueallele=1;
				}
				if (scalar(@FIVallele_left)!=2) {###check if two allele left, return error if not
					print STDERR $FIVsubinfo, "Error: extra at $FIVchrom: $FIVpos\n";###For test###
					print "AllelesLeft: @FIVallele_left\nAlleleNum: $FIVnum_allele_notsure\n";###For test###
					return 1;
				}
				my ($FIVallele_aa,$FIVallele_bb, $FIVallele_dd)=('?', '?', '?');
				if (exists $FIVgenohash{'A'}) {
					if ($FIVmissing_aa ==1 and scalar(@FIVallele_left)>0) {
						$FIVallele_aa=shift @FIVallele_left;
						if ($uniqueallele) {
							${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'}=$FIVallele_aa;
						}
						else {
							$FIVneed_runcompass=1;
						}
					}
					elsif (exists ${$FIVfillinalleles}{$FIVchrom} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'}) {
						$FIVallele_aa=${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'A'};
					}
					else {
						print STDERR $FIVsubinfo, "Error: AA allele at $FIVchrom: $FIVpos\n";###For test###
						return 3;
					}
					push (@FIVarr_geno_fix, $FIVallele_aa);
				}
				if (exists $FIVgenohash{'B'}) {
					if ($FIVmissing_bb ==1 and scalar(@FIVallele_left)>0) {
						$FIVallele_bb=shift @FIVallele_left;
						if ($uniqueallele) {
							${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'}=$FIVallele_bb;
						}
						else {
							$FIVneed_runcompass=1;
						}
					}
					elsif (exists ${$FIVfillinalleles}{$FIVchrom} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'}) {
						$FIVallele_bb=${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'B'};
					}
					else {
						print STDERR $FIVsubinfo, "Error: BB allele at $FIVchrom: $FIVpos\n";###For test###
						return 3;
					}
					push (@FIVarr_geno_fix, $FIVallele_bb);
				}
				if (exists $FIVgenohash{'D'}) {
					if ($FIVmissing_dd ==1 and scalar(@FIVallele_left)>0) {
						$FIVallele_dd=shift @FIVallele_left;
						if ($uniqueallele) {
							${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'}=$FIVmissing_dd;
						}
						else {
							$FIVneed_runcompass=1;
						}
					}
					elsif (exists ${$FIVfillinalleles}{$FIVchrom} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos} and exists ${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'}) {
						$FIVallele_dd=${$FIVfillinalleles}{$FIVchrom}{$FIVpos}{'D'};
					}
					else {
						print STDERR $FIVsubinfo, "Error: DD allele at $FIVchrom: $FIVpos\n";###For test###
						return 3;
					}
					push (@FIVarr_geno_fix, $FIVallele_dd);
				}
			}
			elsif ($FIVnum_allele_notsure==3) {
				@FIVarr_geno_fix=@FIVvcfgeno;
				$FIVneed_runcompass=1;
			}
			${$FIVfixed_geno{$FIVchrom}}{$FIVpos}=join('/', @FIVarr_geno_fix);
#			print $FIVsubinfo, "Test: Fixed Alleles: ".${$FIVfixed_geno{$FIVchrom}}{$FIVpos}."\t$FIVchrom: $FIVpos\n";# if ($debug);
		}
	}
	return (0, $FIVneed_runcompass, \%FIVfixed_geno, $FIVfillinalleles);
}		



### Get polymorphic reads for each geno
### &GetPolymorphicReads($fragment_corrected, $GPRfinalgeno)
### Global: 
### Dependency: 
### Note:
sub GetPolymorphicReads {
	my ($GPRcorrfrag, $GPRfinalgeno)=@_;
	
	my $GPRsubinfo='SUB(GetPolymorphicReads)';
##Format: %GPRexcludedreads=(readid => ++);
	my %GPRexcludedreads=();
##Format: %GPRcoveredsites=(chr => pos => allele => ++);
	my %GPRcoveredsites=();
##Format: %GPRreadcount=('A' => ++, 'B' => ++, 'D' => ++)
	my %GPRreadcount=('A' => 0, 'B' => 0, 'D' => 0);
	
##Format: %{$GPRfinalgeno}=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
##Format: %{$GPRcorrfrag}=($readid => chrom => $partstart => $pos => $allele => $qual);
##Format: %GPRpolymorphismreadids=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++),
	my %GPRpolymorphismreadids=();
	
	my $errorcode=0;
	READID14: {foreach my $GPRthisreadid (keys %{$GPRcorrfrag}) {
		my %GPRassignreads=();
		my %GPRthisreadsites=();
		foreach my $GPRthischrom (keys %{${$GPRcorrfrag}{$GPRthisreadid}}) {
			foreach my $thisstart (keys %{${$GPRcorrfrag}{$GPRthisreadid}{$GPRthischrom}}) {
				foreach my $GPRthispos (keys %{${$GPRcorrfrag}{$GPRthisreadid}{$GPRthischrom}{$thisstart}}) {
					my @GPRthisallelearr=keys %{${$GPRcorrfrag}{$GPRthisreadid}{$GPRthischrom}{$thisstart}{$GPRthispos}};
					unless (scalar(@GPRthisallelearr)==1 and defined $GPRthisallelearr[0] and $GPRthisallelearr[0]=~/^\d+$/) {
						print STDERR $GPRsubinfo, "Error: invalid allele at: Read:Chrom:Pos:Alleles $GPRthischrom:$thisstart:$GPRthispos:@GPRthisallelearr\n";
						$errorcode++;
						last READID14;
					}
					my $thisallelecode=shift @GPRthisallelearr;
					$GPRthisreadsites{$GPRthischrom}{$GPRthispos}{$thisallelecode}++;
					if (exists ${$GPRfinalgeno}{$GPRthischrom} and exists ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}) {
						if (exists ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'A'} and ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'A'}=~/^\d+$/) {
							if ($thisallelecode eq ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'A'}) {
								$GPRassignreads{'A'}{'include'}++;
							}
							else {
								$GPRassignreads{'A'}{'exclude'}++;
							}
						}
						if (exists ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'B'} and ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'B'}=~/^\d+$/) {
							if ($thisallelecode eq ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'B'}) {
								$GPRassignreads{'B'}{'include'}++;
							}
							else {
								$GPRassignreads{'B'}{'exclude'}++;
							}
						}
						if (exists ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'D'} and ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'D'}=~/^\d+$/) {
							if ($thisallelecode eq ${$GPRfinalgeno}{$GPRthischrom}{$GPRthispos}{'D'}) {
								$GPRassignreads{'D'}{'include'}++;
							}
							else {
								$GPRassignreads{'D'}{'exclude'}++;
							}
						}
					}
					else {
						print STDERR $GPRsubinfo, "Error: not existed Read:Chr:Pos $GPRthisreadid:$GPRthischrom:$GPRthispos in \%finalgeno\n";
						$errorcode++;
						next READID14;
					}
				}
			}
		}
		my @GPRthisread2geno=();
		foreach my $GPRgeno (keys %GPRassignreads) {
			unless (exists $GPRassignreads{$GPRgeno}{'include'}) {
				push (@GPRthisread2geno, $GPRgeno);
			}
		}
		if (scalar(@GPRthisread2geno)==0) {
			$GPRexcludedreads{$GPRthisreadid}++;
		}
		else {
			foreach (@GPRthisread2geno) {
				$GPRpolymorphismreadids{$GPRthisreadid}{'gen'}{$_}++;
				$GPRreadcount{$_}++;
			}
			
			foreach my $GPRindchr (keys %GPRthisreadsites) {
				foreach my $GPRindpos (keys %{$GPRthisreadsites{$GPRindchr}}) {
					foreach my $GPRindallele (keys %{$GPRthisreadsites{$GPRindchr}{$GPRindpos}}) {
						$GPRcoveredsites{$GPRindchr}{$GPRindpos}{$GPRindallele}++;
					}
				}
			}
		}
	}}###READID14
	if ($errorcode) {
		print STDERR $GPRsubinfo, "Error: errorcode==1\n";
		return 0;
	}
	if (0) {
		print $GPRsubinfo, "Test: \%GPRpolymorphismreadids\n";
		print Dumper \%GPRpolymorphismreadids;
		print "\n";
		print $GPRsubinfo, "Test: \%GPRexcludedreads\n";
		print Dumper \%GPRexcludedreads;
		print "\n";
		print $GPRsubinfo, "Test: \%GPRcoveredsites\n";
		print Dumper \%GPRcoveredsites;
		print "\n";
	}
	foreach (sort keys %GPRreadcount) {
		print $GPRsubinfo, "Info: Geno:NumReads $_ $GPRreadcount{$_}\n";
	}
	return (1, \%GPRpolymorphismreadids)
}



### convert bam files into fastq
### Bam2FastQ2 ($bamin, 'ABD', %readsum, %allelicreads, outfqprefix_path, [map_code], [MAPQ_code], [RGaware], [path_samtools])
### Global:
### Dependency: $totalreads_AABBDD $reads_per_rg{$_}
### Note: map_code (0=all, 1=mapped, 2=unmapped)
### Note: [RGaware] 0=false; 1=true
sub Bam2FastQ2 {
	my ($BFQbamin, $BFQgeno, $BFQreadsum, $BFQallelicreads, $BFQoutfqprefix, $BFQmapcode, $BFQmapq, $BFQrg_aware, $BFQpath_samtools)=@_;
	
	local *BAMIN; local *BAMOUT; local *BAMKEEP; local *BAMSHARE; local *BAMALLELIC;
	local *ALLFQ; local *SHAREOUT; local *FQOUT;
#Formar: %$BFQreadsum=(readid => ('gen' => ( 'A' => ++, 'B' => ++, 'D' => ++), 
#								'chr' => ( '1AL' => ++, '!BL' => ++, '2DL' => ++, ...),
#								'shared' => ++/0,
#								$chrom => ($pos => (allele => quality))
#								'reads' => (1 => (seq1=>sequence, qual1=> quality), 2=> (seq1=>sequence, qual1=> quality), U => ...)
#								'excluded' => 0/++
#							)
#					);
	my $BFQsubinfo='SUB(Bam2FastQ2)';
	$BFQpath_samtools='samtools' unless (defined $BFQpath_samtools);
	$BFQmapcode=0 unless (defined $BFQmapcode);
	$BFQmapq=0 unless (defined $BFQmapq);
	$BFQrg_aware=0 unless (defined $BFQrg_aware);
	my %BFQexcludedreads=();
	my $BamKit_failure=0;
	my $BamKit_success=1;
##Format: %returnhash=(A => "A.fastq", B => b.fastq, D => D.fastq);
	my %returnhash=();
##Format: %BFQfpkms=(chr => ('chr1' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#											RG2 => ...)
	my %BFQfpkms=();
##Format %BFQclusterfpkm=('genomesize' => size, 
#							'cluster' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#											RG2 => ...)
#Calculate overal FPKM using estimated Cluster size
	my %BFQclusterfpkm=();
##Format: %BFQreadpair=(chr => (RG1 => (readid => 1 => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
#													2 => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
#													'U' => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
#								RG2...);
	my %BFQreadpair=();
	my %BFQgenos=();
	$BFQoutfqprefix='./' unless (defined $BFQoutfqprefix);
	my %BFQreflen=();
	my $BFQouttestbam=1; ### Output BAMs for debugging
	
#	print "${BFQsubinfo}Test: \n\tBAM: $BFQbamin\n\tGeno: $BFQgeno\n\tFqPrefix: $BFQoutfqprefix\n\tMapcode: $BFQmapcode\n\tMapQ: $BFQmapq\n\tReadGroupAware: $BFQrg_aware\n\tSAMtools: $BFQpath_samtools\n";### For test ###
	
	unless (defined $BFQbamin and -s $BFQbamin) {
		print STDERR "${BFQsubinfo}Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $BFQgeno and $BFQgeno =~/^[ABD]+$/) {
		print STDERR "${BFQsubinfo}Error: unknown geno code: $BFQgeno\n";
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
	
	if ($BFQouttestbam) {
		close BAMKEEP if (defined fileno(BAMKEEP));
		unless (open (BAMKEEP, " | samtools view -bhS - > $BFQbamin.included.bam")) {
			print STDERR "${BFQsubinfo}Error: write BAM: $BFQbamin.included.bam \n";
			return $BamKit_failure;
		}
		close BAMSHARE if (defined fileno(BAMSHARE));
		unless (open(BAMSHARE, " | samtools view -bhS - > $BFQbamin.shared.bam")) {
			print STDERR "${BFQsubinfo}Error: write BAM: $BFQbamin.shared.bam \n";
			return $BamKit_failure;
		}
		close BAMALLELIC if (defined fileno(BAMALLELIC));
		unless (open(BAMALLELIC, " | samtools view -bhS - > $BFQbamin.allelic.bam")) {
			print STDERR "${BFQsubinfo}Error: write BAM: $BFQbamin.allelic.bam \n";
			return $BamKit_failure;
		}
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
			if ($BFQouttestbam) {
				print BAMKEEP $BFQline1."\n";
				print BAMSHARE $BFQline1."\n";
				print BAMALLELIC $BFQline1."\n";
			}
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
		unless (exists ${$BFQreadsum}{$BFQarr[0]}) {
			print BAMOUT $BFQline1."\n";
			next;
		}
		if (exists ${$BFQreadsum}{$BFQarr[0]}{'excluded'} and ${$BFQreadsum}{$BFQarr[0]}{'excluded'}>0) {
			print BAMOUT $BFQline1."\n";
			next;
		}
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
				$BFQexcludedreads{$BFQarr[0]}++;
				print BAMOUT $BFQline1."\n";
				next;
			}
			foreach (@BFQtemp_gens) {
				push (@BFQthisread_gens, $_) if (exists $BFQgenos{$_});
			}
			if (scalar(@BFQthisread_gens)<1 or scalar(@BFQthisread_gens)>3) {
#				print STDERR "${BFQsubinfo}Warnings: no subgenome assignments for read $BFQarr[0] at Chrom: Pos $BFQarr[2]:$BFQarr[3]\n";### For test ###
				$BFQexcludedreads{$BFQarr[0]}++;
				print BAMOUT $BFQline1."\n";
				next;
			}
		}
		
		if ($BFQis_shared==0 and $BFQis_allelic==0) {
#			print STDERR "${BFQsubinfo}Warnings: read ($BFQarr[0]) is neither shared or allelic at Chrom:Pos $BFQarr[2]:$BFQarr[3]\n";
			print BAMOUT $BFQline1."\n";
			$BFQexcludedreads{$BFQarr[0]}++;
			next;
		}
		elsif ($BFQis_shared==1 and $BFQis_allelic==1) {
			print STDERR "${BFQsubinfo}Warnings: read ($BFQarr[0]) is both shared and allelic at Chrom:Pos $BFQarr[2]:$BFQarr[3]\n";
			print BAMOUT $BFQline1."\n";
			$BFQexcludedreads{$BFQarr[0]}++;
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
		if ($BFQouttestbam) {
			if ($BFQis_shared==1) {
				print BAMSHARE $BFQline1."\n";
			}
			elsif ($BFQis_allelic==1) {
				print BAMALLELIC $BFQline1."\n";
			}
			print BAMKEEP $BFQline1."\n";
		}
	}
	close BAMIN;
	close BAMOUT;
	if ($BFQouttestbam) {
		close BAMKEEP;
		close BAMSHARE;
		close BAMALLELIC;
	}
	print "${BFQsubinfo}Warnings: excluded reads num: ".scalar(keys %BFQexcludedreads)."\n";

###Output fastq (shared)
	my $BFQshared_fastq=$BFQoutfqprefix.'.shared.fastq';
	my $BFQallfastq=$BFQoutfqprefix.'.all.fastq';
	unlink $BFQshared_fastq if (-e $BFQshared_fastq);
	
	my $BFQnumshared=0;
	close ALLFQ if (defined fileno(ALLFQ));
	unless (open(ALLFQ, "> $BFQallfastq")) {
		print STDERR "${BFQsubinfo}Error: can not write ALL Fastq: $BFQallfastq\n";
		return $BamKit_failure;
	}
	close SHAREOUT if (defined fileno(FQOUT));
	unless (open (SHAREOUT, ">$BFQshared_fastq")) {
		print STDERR "${BFQsubinfo}Error: can not write SHARED Fastq: $BFQshared_fastq\n";
		return $BamKit_failure;
	}
	foreach my $BFQreadid_shared (keys %{$BFQreadsum}) {
		if (exists ${${$BFQreadsum}{$BFQreadid_shared}}{'shared'} and ${${$BFQreadsum}{$BFQreadid_shared}}{'shared'}>0 and exists ${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}) {
			$BFQnumshared++;
			if (exists ${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{1}) {
				my ($BFQseq, $BFQqual)=@{${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{1}};
				print SHAREOUT '@'.$BFQreadid_shared."/1\n$BFQseq\n+\n$BFQqual\n";
				print ALLFQ '@'.$BFQreadid_shared."/1\n$BFQseq\n+\n$BFQqual\n";
			}
			if (exists ${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{2}) {
				my ($BFQseq, $BFQqual)=@{${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{2}};
				print SHAREOUT '@'.$BFQreadid_shared."/2\n$BFQseq\n+\n$BFQqual\n";
				print ALLFQ '@'.$BFQreadid_shared."/2\n$BFQseq\n+\n$BFQqual\n";
			}
			if (exists ${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{'U'}) {
				my ($BFQseq, $BFQqual)=@{${${${$BFQreadsum}{$BFQreadid_shared}}{'reads'}}{'U'}};
				print SHAREOUT '@'.$BFQreadid_shared."\n$BFQseq\n+\n$BFQqual\n";
				print ALLFQ '@'.$BFQreadid_shared."\n$BFQseq\n+\n$BFQqual\n";
			}
		}
	}
	close SHAREOUT;
###Output Allele FQ
	my %BFQcount_geno2reads=('A' => 0, 'B' => 0, 'D' => 0);
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
			next unless (exists ${$BFQallelicreads}{$BFQreadid_allelic}{'gen'} and exists ${$BFQallelicreads}{$BFQreadid_allelic}{'gen'}{$BFQind_gen} and ${$BFQallelicreads}{$BFQreadid_allelic}{'gen'}{$BFQind_gen}>0);
			$BFQcount_geno2reads{$BFQind_gen}++;
			if (exists ${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}) {
				if (exists ${${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}}{1}) {
					my ($BFQseq, $BFQqual)=@{${${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}}{1}};
					print FQOUT '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
					print ALLFQ '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
				}
				if (exists ${${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}}{2}) {
					my ($BFQseq, $BFQqual)=@{${${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}}{2}};
					print FQOUT '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
					print ALLFQ '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
				}
				if (exists ${${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}}{'U'}) {
					my ($BFQseq, $BFQqual)=@{${${${$BFQallelicreads}{$BFQreadid_allelic}}{'reads'}}{'U'}};
					print FQOUT '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
					print ALLFQ '@'.$BFQreadid_allelic."/1\n$BFQseq\n+\n$BFQqual\n";
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
	close ALLFQ;
	foreach (sort keys %BFQcount_geno2reads) {
		print STDERR $BFQsubinfo, "Info: NumReads(Geno:Shared:Allelic) $_:$BFQnumshared:$BFQcount_geno2reads{$_}\n";
	}
	
### Calculate GetGenomeSize
	my ($BFQtest_cmd_genomesize, $BFQgenomesize)=&GetGenomeSize($BFQallfastq);
	unless ($BFQtest_cmd_genomesize==1 and $BFQgenomesize=~/^\d+$/ and $BFQgenomesize>0) {
		print STDERR "${BFQsubinfo}Error: unknown genome size\n";
		return $BamKit_failure;
	}
	print "${BFQsubinfo}Info: Estimated Cluster size: $BFQgenomesize\n";
	
### Paired reads Counts
##Format: %BFQreadpair=(chr => (RG1 => (readid => 1 => ('A' => 1/0, 'B' =>1/0, 'D' => 1/0, 'shared' => 1/0)
##Format: %BFQreadcounts=(chr => ('chr' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
	my %BFQreadcounts=();
#Format: %BFQoveralcounts=('cluster' => FPKM
#							'read_group' => (RG1 => ('global' => FPKM, 'A' => FPKM, 'B' => FPKM, 'D' => FPKM)
#							'allelic' => ('A' => FPKM, 'B' => FPKM, 'D' => FPKM)
	my %BFQoveralcounts=();
##Format: %BFQreadcounted=(readid => ++)
#test readid if counted or not for %BFQreadcounted
	my %BFQreadcounted=();
	
	foreach my $BFQchrom (keys %BFQreadpair) {
		foreach my $BFQrg (keys %{$BFQreadpair{$BFQchrom}}) {
			foreach my $BFQreadid (keys %{${$BFQreadpair{$BFQchrom}}{$BFQrg}}) {
				if (exists ${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1} and exists ${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{2}) {
					${$BFQreadcounts{$BFQchrom}}{'chr'}++;
					${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQrg}}{'global'}++;
					unless (exists $BFQreadcounted{$BFQreadid}) {
						$BFQoveralcounts{'cluster'}++;
						$BFQoveralcounts{'read_group'}{$BFQrg}{'global'}++;
					}
					if (exists ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{'shared'} and ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{'shared'} >0) {###shared
						foreach (keys %BFQgenos) {
							${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQrg}}{$_}++;
							${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_}++;
							unless (exists $BFQreadcounted{$BFQreadid}) {
								$BFQoveralcounts{'read_group'}{$BFQrg}{$_}++;
								$BFQoveralcounts{'allelic'}{$_}++;
								$BFQreadcounted{$BFQreadid}++;
							}
						}
					}
					else {##allelic
						foreach (keys %BFQgenos) {
							if (exists ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{$_} and ${${${${$BFQreadpair{$BFQchrom}}{$BFQrg}}{$BFQreadid}}{1}}{$_} >0){
								${${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQrg}}{$_}++;
								${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_}++;
								unless (exists $BFQreadcounted{$BFQreadid}) {
									$BFQoveralcounts{'read_group'}{$BFQrg}{$_}++;
									$BFQoveralcounts{'allelic'}{$_}++;
									$BFQreadcounted{$BFQreadid}++;
								}
							}
						}
					}
				}
			}
		}
	}
	
### Calculate FPKMs for each reference sequence
	foreach my $BFQchrom (keys %BFQreadcounts) {
		if (exists ${$BFQreadcounts{$BFQchrom}}{'chr'}) {
			if (exists $BFQreflen{$BFQchrom} and $BFQreflen{$BFQchrom}>0) {
				${$BFQfpkms{$BFQchrom}}{'chr'}=${$BFQreadcounts{$BFQchrom}}{'chr'} * 10^9 / ($BFQreflen{$BFQchrom} * $totalreads_AABBDD);
			}
			else {
				print STDERR "${BFQsubinfo}Error: can not get ref ($BFQchrom) length from header of BAM: $BFQbamin\n";
				${$BFQfpkms{$BFQchrom}}{'chr'}='undef';
			}
		}
		foreach ('A', 'B', 'D') {
			if (exists ${$BFQreadcounts{$BFQchrom}}{'allelic'} and exists ${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_}) {
				${${$BFQfpkms{$BFQchrom}}{'allelic'}}{$_}= ${${$BFQreadcounts{$BFQchrom}}{'allelic'}}{$_} * 10^9 / ($BFQreflen{$BFQchrom} * $totalreads_AABBDD);
			}
			else {
				${${$BFQfpkms{$BFQchrom}}{'allelic'}}{$_}=0;
			}
		}
		foreach my $BFQindrg (keys %reads_per_rg) {
			if (exists ${$BFQreadcounts{$BFQchrom}}{'read_group'} and exists ${${$BFQreadcounts{$BFQchrom}}{'read_group'}}{$BFQindrg}) {
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
	
	
### Calculate overal FPKMs into %BFQclusterfpkm %BFQoveralcounts
	if (defined $BFQgenomesize and $BFQgenomesize>0) {
		$BFQclusterfpkm{'genomesize'}=$BFQgenomesize;
		if (exists $BFQoveralcounts{'cluster'}) {
			$BFQclusterfpkm{'cluster'}=$BFQoveralcounts{'cluster'} * 10^9 / ($BFQgenomesize * $totalreads_AABBDD);
		}
		else {
			$BFQclusterfpkm{'cluster'}=0;
		}
		foreach my $BFQindrg (keys %reads_per_rg) {
			if (exists $BFQoveralcounts{'read_group'} and exists $BFQoveralcounts{'read_group'}{$BFQindrg}) {
				foreach ('global', 'A', 'B', 'D') {
					if (exists $BFQoveralcounts{'read_group'}{$BFQindrg}{$_}) {
						$BFQclusterfpkm{'read_group'}{$BFQindrg}{$_}=$BFQoveralcounts{'read_group'}{$BFQindrg}{$_}  * 10^9 / ($BFQgenomesize * $reads_per_rg{$BFQindrg});
					}
					else {
						$BFQclusterfpkm{'read_group'}{$BFQindrg}{$_}=0;
					}
				}
			}
			else {
				foreach ('global', 'A', 'B', 'D') {
					$BFQclusterfpkm{'read_group'}{$BFQindrg}{$_}=0;
				}
			}
		}
		foreach ('A', 'B', 'D') {
			if (exists $BFQoveralcounts{'allelic'} and exists $BFQoveralcounts{'allelic'}{$_}) {
				$BFQclusterfpkm{'allelic'}{$_}= $BFQoveralcounts{'allelic'}{$_} * 10^9 / ($BFQgenomesize * $totalreads_AABBDD);
			}
			else {
				$BFQclusterfpkm{'allelic'}{$_}=0;
			}
		}
	}
	else {
		$BFQclusterfpkm{'genomesize'}=0;
	}

	%BFQoveralcounts=();
	%BFQreadcounted=();
	%BFQreadcounts=();
	if (0) {###Test cluster FPKM hash %BFQclusterfpkm
		print "${BFQsubinfo}Test: Reference FPKM\n";
		print Dumper \%BFQfpkms;
		print "\n";
	}
	if (0) {###Test cluster FPKM hash %BFQclusterfpkm
		print "${BFQsubinfo}Test: Cluster FPKM\n";
		print Dumper \%BFQclusterfpkm;
		print "\n";
	}
	return ($BamKit_success, \%returnhash, \%BFQfpkms, \%BFQclusterfpkm);
}



### Estimate genome size by kmergenie
### &GetGenomeSize (fastq, addcmd, path_kmergenie)
### Global:
### Dependency:
### Note:
sub GetGenomeSize {
	my ($GGSfastq, $GGSprefix, $GGSaddcmd, $GGSpath_kmergenie)=@_;
	
	my $GGSsubinfo='SUB(GetGenomeSize)';
	$GGSprefix="$RunDir/kmergenie/MyFastq" unless (defined $GGSprefix);
	$GGSpath_kmergenie='kmergenie' unless (defined $GGSpath_kmergenie);
	$GGSaddcmd= ' -l 27 -k 87 -o '.$GGSprefix.' ' unless (defined $GGSaddcmd);
	my $GGSgenomesize=0;
	local *GGSIN;
	
	my @GGoldfiles=glob "$RunDir/kmergenie/*";
	foreach (@GGoldfiles){
		unlink @GGoldfiles if (-e @GGoldfiles);
	}
	
	unless (defined $GGSfastq and -s $GGSfastq) {
		print STDERR "${GGSsubinfo}Error: invalid fastq input\n";
		return 0;
	}

	unless (exec_cmd_return("$GGSpath_kmergenie $GGSfastq $GGSaddcmd >>kmergenie.log 2>>kmergenie.err")) {
		print STDERR "${GGSsubinfo}Error: invalid fastq input\n";
		unlink glob("${GGSprefix}*");
		return 0;
	}
	if (-s "${GGSprefix}_report.html") {
		close GSSIN if (defined fileno(GSSIN));
		unless (open (GGSIN, "< ${GGSprefix}_report.html")) {
			print STDERR "${GGSsubinfo}Error: can not open kmergenie report: ${GGSprefix}_report.html\n";
			unlink glob("${GGSprefix}*");
			return 0;
		}
		while (my $GGSline=<GGSIN>) {
			if ($GGSline=~/Predicted\s+assembly\s+size:\s+(\d+)\s+bp/) {
				$GGSgenomesize=$1;
				last;
			}
		}
		close GSSIN;
	}
	if ($GGSgenomesize>0) {
		unlink glob("${GGSprefix}*");
		return (1, $GGSgenomesize);
	}
	else {
		print STDERR "${GGSsubinfo}Error: failed kmergenie genome size\n";
		unlink glob("${GGSprefix}*");
		return 0;
	}
}
