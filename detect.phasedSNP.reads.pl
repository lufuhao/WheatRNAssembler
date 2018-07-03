#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use Bio::DB::Sam;
use Scalar::Util;
use File::Copy 'move';
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150112

Requirements:
	Programs: hapcompass.jar hc2vcf.jar
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;

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
our ($help, $verbose, $ver);
our ($reference, $input_vcf);
our ($file_bam_aabbdd, $file_bam_aabb, $file_bam_aa, $file_bam_dd);
our ($file_cluster, $cluster_log, $fpkm_log);

=head6
GetOptions(
	"help|h!" => \$help,
	"bam|b:s" => \$reference,
	"vcf|f:s" => \$input_vcf,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;
=cut


### Defaults ########################################################
our $RootDir=$Bin;



our $debug=0;

our $vcf=$ARGV[0];
our $sam=$ARGV[1];
our $ref=$ARGV[2];
###Cluster QC
$cluster_log="Cluster.log" unless (defined $cluster_log);
#SUB(GroupVcf);
our $min_share_alignments=3;
###Express
our $path_express='express' unless (defined $path_express);
our $express_frag_len_mean=250 unless (defined $express_frag_len_mean);
our $express_frag_len_stddev=300 unless (defined $express_frag_len_stddev);
our $express_max_read_len=100 unless (defined $express_max_read_len);
###Freebayes
our $path_freebayes='freebayes' unless (defined $path_freebayes);
our $freebayes_min_coverage=3 unless (defined $freebayes_min_coverage);
our $freebayes_min_alternative_count=3 unless (defined $freebayes_min_alternative_count);
our $min_mapq=0 unless (defined $min_mapq);
#samtools
our $path_samtools='samtools';
#HapCompass
our $path_hapcompassjar="$RootDir/utils/hapcompass/hapcompass.jar";
our $path_hc2vcf="$RootDir/utils/hapcompass/hc2vcf.jar";
our $path_java='java';
#VCFtools
our $path_vcfmerge='vcf-merge';



### input and output ################################################
#Check important file inpput
die "Please specify reference file\n" unless (defined $reference and -s $reference);
die "Please specify cluster file\n" unless (defined $file_cluster and -s $file_cluster);
die "Please specify the file of bam list for AABBDD\n" unless (defined $file_bam_aabbdd and -s $file_bam_aabbdd);
die "Please specify the file of bam list for AABB\n" unless (defined $file_bam_aabb and -s $file_bam_aabb);
die "Please specify the file of bam list for AA\n" unless (defined $file_bam_aa and -s $file_bam_aa);
die "Please specify the file of bam list for DD\n" unless (defined $file_bam_dd and -s $file_bam_dd);

#Remove last-run cluster log 
unlink $cluster_log if (-e $cluster_log);
unlink $fpkm_log if (-e $fpkm_log);

#read AABBDD bam files
open (BAMAABBDD, "<$file_bam_aabbdd") || die "InputOutputError: can not open AABBDD $file_bam_aabbdd\n";
our @bam_AABBDD=();
while (my $line1=<BAMAABBDD>) {
	chomp $line1;
	if (defined $line1 and $line1 ne '' and -s $line1) {
		push (@bam_AABBDD, $line1);
	}
	else {
		die "InputOutputError: bam $line1 in AABBDD $file_bam_aabbdd not exist or empty\n";
	}
}
die "InputOutputError: empty AABBDD bam files" if (scalar(@bam_AABBDD)<1);
if ($verbose) {
	print "InputOutputReport: Reading AABBDD bams: ".scalar(@bam_AABBDD)."\n";
	foreach (@bam_AABBDD) {
		print "---> $_\n";
	}
}
close BAMAABBDD;

# read AABB bamfiles
open (BAMAABB, "<$file_bam_aabb") || die "InputOutputError: can not open AABB $file_bam_aabb\n";
our @bam_AABB=();
while (my $line2=<BAMAABB>) {
	chomp $line2;
	if (defined $line2 and $line2 ne '' and -s $line2) {
		push (@bam_AABB, $line2);
	}
	else {
		die "InputOutputError: bam $line2 in AABB $file_bam_aabb not exist or empty\n";
	}
}
die "InputOutputError: empty AABB bam files" if (scalar(@bam_AABB)<1);
close BAMAABB;
if ($verbose) {
	print "InputOutputReport: Reading AABB bams: ".scalar(@bam_AABB)."\n";
	foreach (@bam_AABB) {
		print "---> $_\n";
	}
}

# read AA bam files
open (BAMAA, "<$file_bam_aa") || die "InputOutputError: can not open AA $file_bam_aa\n";
our @bam_AA=();
while (my $line3=<BAMAA>) {
	chomp $line3;
	if (defined $line3 and $line3 ne '' and -s $line3) {
		push (@bam_AA, $line3);
	}
	else {
		die "InputOutputError: bam $line3 in AA $file_bam_aa not exist or empty\n";
	}
}
die "InputOutputError: empty AA bam files" if (scalar(@bam_AA)<1);
close BAMAA;
if ($verbose) {
	print "InputOutputReport: Reading AA bams: ".scalar(@bam_AA)."\n";
	foreach (@bam_AA) {
		print "---> $_\n";
	}
}

# read DD files
open (BAMDD, "<$file_bam_dd") || die "InputOutputError: can not open DD $file_bam_dd\n";
our @bam_DD=();
while (my $line4=<BAMDD>) {
	chomp $line4;
	if (defined $line4 and $line4 ne '' and -s $line4) {
		push (@bam_DD, $line4);
	}
	else {
		die "InputOutputError: bam $line4 in DD $file_bam_dd not exist or empty\n";
	}
}
die "InputOutputError: empty DD bam files" if (scalar(@bam_DD)<1);
close BAMDD;
if ($verbose) {
	print "InputOutputReport: Reading DD bams: ".scalar(@bam_DD)."\n";
	foreach (@bam_DD) {
		print "---> $_\n";
	}
}



### Main ############################################################
our $RunDir=getcwd;
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
our $fasta_index=&CdbFasta($reference);

# Load clustered fasta id list, 1 cluster perl line, saparated by space
open (CLUSTER, "<$file_cluster") || die "(Main)Error: can not open file: $file_cluster\n";
# Output cluster succeeds or not. Format: 
open (CLUSTERLOG, ">$cluster_log") || die "(Main)Error: can not write file: $cluster_log\n";
# output fpkm calculated
open (FPKMLOG, ">$fpkm_log") || die "(Main)Error: can not write file: $fpkm_log\n";
our $cluster_linenum=0;###for quick find which line/cluster has problem
our @cluster_seqids=();###
while (my $cluster_line=<CLUSTER>) {
	chomp $cluster_line;
	$cluster_linenum++;
	print STDOUT "\n\n\n\n\n##### Prcessing Cluster $cluster_linenum ###\n";
	print STDERR "\n\n\n\n\n##### Prcessing Cluster $cluster_linenum ###\n";
	@cluster_seqids=(); ###Empty this in case of abnormal duplicates
	@cluster_seqids=split(/\s+/, $cluster_line);
#Check if empty line
	if (scalar(@cluster_seqids)<1) {
		print STDERR "MainWarnings: line $cluster_linenum in $file_cluster ignored as empty\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	my $bam_extract_region=join " ", @cluster_seqids;
#extract fasta reference
	unless (mkdir ("$RunDir/Clust$cluster_linenum", 0766)) {
		print STDERR "(Main)Error: create folder $RunDir/Clust$cluster_linenum\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	if (&CdbYank($fasta_index, "Clust$cluster_linenum/ref.$cluster_linenum.fa", @cluster_seqids)) {
		print STDERR "(Main)Error: failed extract ClusterID $cluster_linenum: @cluster_seqids\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
#Empty folder
	unlink glob "$RunDir/AABBDD/*";
	unlink glob "$RunDir/AABB/*";
	unlink glob "$RunDir/AA/*";
	unlink glob "$RunDir/DD/*";
#Extract bam from @bam_AABBDD, @bam_AABB, @bam_AA, @bam_DD
	my $i=0;
	my $bam_failure=0;
	foreach (@bam_AABBDD) {
		$i++;
		my $cmd="$path_samtools view -h $_ $bam_extract_region > AABBDD/AABBDD.$cluster_linenum.$i.bam";
		if (&exec_cmd_return($cmd)) {
			print STDERR "SUB(Main): AABBDD bam extract $bam_extract_region (Cluster: $cluster_linenum) from $_ failed\n";
			$bam_failure=1;
		}
	}
	$i=0;
	foreach (@bam_AABB) {
		$i++;
		my $cmd="$path_samtools view -h $_ $bam_extract_region > AABB/AABB.$cluster_linenum.$i.bam";
		if (&exec_cmd_return($cmd)) {
			print STDERR "SUB(Main): AABBDD bam extract $bam_extract_region (Cluster: $cluster_linenum) from $_ failed\n";
			$bam_failure=1;
		}
	}
	$i=0;
	foreach (@bam_AA) {
		$i++;
		my $cmd="$path_samtools view -h $_ $bam_extract_region > AA/AA.$cluster_linenum.$i.bam";
		if (&exec_cmd_return($cmd)) {
			print STDERR "SUB(Main): AABBDD bam extract $bam_extract_region (Cluster: $cluster_linenum) from $_ failed\n";
			$bam_failure=1;
		}
	}
	$i=0;
	foreach (@bam_DD) {
		$i++;
		my $cmd="$path_samtools view -h $_ $bam_extract_region > DD/DD.$cluster_linenum.$i.bam";
		if (&exec_cmd_return($cmd)) {
			print STDERR "SUB(Main): AABBDD bam extract $bam_extract_region (Cluster: $cluster_linenum) from $_ failed\n";
			$bam_failure=1;
		}
	}
	if ($bam_failure) {
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	else {###Extract Bam and index
		if (&exec_cmd_return("$path_samtools merge Clust$cluster_linenum/AABBDD.$cluster_linenum.bam AABBDD/*.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}elsif (&exec_cmd_return("$path_samtools index Clust$cluster_linenum/AABBDD.$cluster_linenum.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}
		if (&exec_cmd_return("$path_samtools merge Clust$cluster_linenum/AABB.$cluster_linenum.bam AABB/*.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}elsif (&exec_cmd_return("$path_samtools index Clust$cluster_linenum/AABB.$cluster_linenum.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}
		if (&exec_cmd_return("$path_samtools merge Clust$cluster_linenum/AA.$cluster_linenum.bam AA/*.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}elsif (&exec_cmd_return("$path_samtools index Clust$cluster_linenum/AA.$cluster_linenum.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}
		if (&exec_cmd_return("$path_samtools merge Clust$cluster_linenum/DD.$cluster_linenum.bam DD/*.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}elsif (&exec_cmd_return("$path_samtools index Clust$cluster_linenum/DD.$cluster_linenum.bam")) {
			print CLUSTERLOG $cluster_line."\tFail\t1\n";
			next;
		}
	}
# Empty bamfiles
	unlink glob "$RunDir/AABBDD/*" if (-s "Clust$cluster_linenum/AABBDD.$cluster_linenum.bam");
	unlink glob "$RunDir/AABB/*" if (-s "Clust$cluster_linenum/AABB.$cluster_linenum.bam");
	unlink glob "$RunDir/AA/*" if (-s "Clust$cluster_linenum/AA.$cluster_linenum.bam");
	unlink glob "$RunDir/DD/*" if (-s "Clust$cluster_linenum/DD.$cluster_linenum.bam");
# Call SNP and prepahse
	if (&RunFreebayes("Clust$cluster_linenum/ref.$cluster_linenum.fa", "Clust$cluster_linenum/AABBDD.$cluster_linenum.bam", "Clust$cluster_linenum/AABB.$cluster_linenum.bam", "Clust$cluster_linenum/AA.$cluster_linenum.bam", "Clust$cluster_linenum/DD.$cluster_linenum.bam")==1) {
		#do sth
	};
}
close CLUSTER;
close CLUSTERLOG;
close FPKMLOG;

#my $readvcf=&ReadVcf($vcf);
#my $read_phased_vcf=&GroupVcf($readvcf, $sam, $ref);
#&AsignABDallele($readvcf);





#####################################################################
###                         sub functions                         ###
#####################################################################



###ExpressRPKM
###&ExpressFpkm($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, @EFsamfiles)
###Global:$path_express, $cluster_linenum, @cluster_seqids
###Dependancy: &exec_cmd, &mytime, &RetrvNoExt
sub ExpressFpkm {
	my ($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, @EFsamfiles)=@_;
	my @EFfpkms=();
	my $EFi=0;
#@EFsamfiles=(bam_AABBDD, bam_AABB, bam_AA, bam_DD) in this order
	foreach my $EFbamfile (@EFsamfiles) {
		my $EFbamfile_base=&RetrvNoExt($EFbamfile);
		unless (-d $EFbamfile_base) {
			mkdir ("$EFbamfile_base", 0777) || die "SUB(ExpressFpkm)Error: can not create folder $EFbamfile_base\n";
		}
		unlink glob "$EFbamfile_base/*.xprs";###delete last-run files
		if (-s $EFbamfile) {
			my $cmd="$path_express --frag-len-mean $EFfrag_len_mean --frag-len-stddev $EFfrag_len_stddev --max-read-len $EFmax_read_len --output-dir $EFbamfile_base $EFref $EFbamfile";
			&exec_cmd($cmd);
			open (FPKM, "$EFbamfile_base/results.xprs") || die "SUB(ExpressFpkm)Error: can not open file $EFbamfile_base/results.xprs\n";
			my $EFline=<FPKM>;###firstline is header
			while ($EFline=<FPKM>) {
				chomp $EFline;
				my @arr=split(/\t/, $EFline);
#$arr[1] is reference sequence ID
#$arr[10] is the estimated relative abundance in units of fragments per kilobase per million mapped
				${$EFfpkms[$EFi]}{$arr[1]}=$arr[10];
			}
			close FPKM;
			unlink glob "$EFbamfile_base/*.xprs";###delete last-run files
			$EFi++;
			####Delete directory###PAUSE
		}
		else {
			die "SUB(ExpressFpkm)Error: bam file $EFbamfile empty or non-exists\n";
		}
	}
#Check if each @cluster_seqids have a FPKM value, 0 if not;
#Output FPKM to $fpkm_log file
	print FPKMLOG "\n\n\n\n\n##### Cluster ID: $cluster_linenum #####\n";
	foreach my $EFseqid (@cluster_seqids) {
		my $EFfpkm_logprint=$EFseqid;
		for (my $i=0; $i<scalar(@EFfpkms); $i++) {
			if (! exists ${$EFfpkms[$i]}{$EFseqid}) {
				${$EFfpkms[$i]}{$EFseqid}=0;
			}
			$EFfpkm_logprint.="\t".${$EFfpkms[$i]}{$EFseqid};
			print "SUB(ExpressFpkm)Test: Bam $EFsamfiles[$i]: Ref $EFseqid: ${$EFfpkms[$i]}{$EFseqid}\n" if ($debug); ###Test###
		}
		print FPKMLOG $EFfpkm_logprint."\n";
	}
	return \@EFfpkms;
	
#Return format
#@EFfpkms=(
#			(ref1 => fpkm1, ref2 => fpkm2, ...), 
#			(ref1 => fpkm1, ref2 => fpkm2, ...), 
#			...)
#$EFfpkms[$ith]->{reference_id}=fpkm
}



###RunFreebayes, return merged vcf and raw AABBDD vcf
###&RunFreebayes($RFclusterID, ref_fasta, bam_AABBDD, bam_AABB, bam_AA, bam_DD)
###Global: $freebayes_min_coverage, $freebayes_min_alternative_count, $min_mapq,$express_frag_len_mean, $express_frag_len_stddev, $express_max_read_len, $cluster_linenum, @cluster_seqids, $RunDir, $path_vcfmerge
###Dependancy: &ReadSam,
sub RunFreebayes {
	my ($RFfile_reference, $RFbam_aabbdd, $RFbam_aabb, $RFbam_aa, $RFbam_dd)=@_;
	my $RFcmd='';
	my @RF_express=();
	my %RFploidy=0;
#calculate RPKM into @RF_express
	my $RFexpress_obj=&ExpressFpkm($RFfile_reference, $express_frag_len_mean, $express_frag_len_stddev, $express_max_read_len, $RFbam_aabbdd, $RFbam_aabb, $RFbam_aa, $RFbam_dd);
#AABBDD:	$RF_express[0]->{ref}=fpkm
#AABB:		$RF_express[1]->{ref}=fpkm
#AA:		$RF_express[2]->{ref}=fpkm
#DD:		$RF_express[3]->{ref}=fpkm
	@RF_express=@{$RFexpress_obj};
#generate AABBDD vcf for freebayes --variant-input
	$RFcmd="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --pooled-continuous --min-mapping-quality $min_mapq $RFbam_aabbdd | cut -f -8 > AABBDD.guide.vcf";
	&exec_cmd($RFcmd);
	unless (-s "AABBDD.guide.vcf") {
		print STDERR "SUB(RunFreebayes)Error: can not generate --variant-input file for cluster $cluster_linenum\n";
		return 1;
	}
###generate vcf for each reference
	unlink glob "AABBDD/*";
	unlink glob "AABB/*";
	unlink glob "AA/*";
	unlink glob "DD/*";
	foreach my $RFind_ref (@cluster_seqids) {
##COMMENT: check if AABBDD transcript
		if (! defined ${$RF_express[0]}{$RFind_ref} or ${$RF_express[0]}{$RFind_ref}==0) {
			print STDERR "SUB(RunFreebayes)Error: $RFind_ref FPKM problem: non-exist or 0\n";
			next;
		}
		${$RF_express[1]}{$RFind_ref}=0 if (! exists ${$RF_express[1]}{$RFind_ref});
		${$RF_express[2]}{$RFind_ref}=0 if (! exists ${$RF_express[2]}{$RFind_ref});
		${$RF_express[3]}{$RFind_ref}=0 if (! exists ${$RF_express[3]}{$RFind_ref});
##COMMENT: check ploidy
		my ($RFaabb_expressed, $RFaa_expressed, $RFdd_expressed)=(0, 0, 0);
		if (${$RF_express[1]}{$RFind_ref}==0) {
			$RFaabb_expressed=0;
		}
		elsif (${$RF_express[1]}{$RFind_ref}>0) {
			$RFaabb_expressed=1;
		}
		else {print STDERR "SUB(RunFreebayes)Error: $RFind_ref FPKM problem in AABB: unknown FPKM value: ${$RF_express[1]}{$RFind_ref}\n";}

		if (${$RF_express[2]}{$RFind_ref}==0) {
			$RFaa_expressed=0;
		}
		elsif (${$RF_express[2]}{$RFind_ref}>0) {
			$RFaa_expressed=1;
		}
		else {print STDERR "SUB(RunFreebayes)Error: $RFind_ref FPKM problem in AA: unknown FPKM value: ${$RF_express[2]}{$RFind_ref}\n";}

		if (${$RF_express[3]}{$RFind_ref}==0) {
			$RFdd_expressed=0;
		}
		elsif (${$RF_express[3]}{$RFind_ref}>0) {
			$RFdd_expressed=1;
		}
		else {print STDERR "SUB(RunFreebayes)Error: $RFind_ref FPKM problem in DD: unknown FPKM value: ${$RF_express[3]}{$RFind_ref}\n";}
##COMMENT: define ploidy
		if (exists $RFploidy{$RFind_ref}) {
			print STDERR "SUB(RunFreebayes)Error: $RFind_ref ploidy value already exists\n";
			next;
		}
		$RFploidy{$RFind_ref}=$RFaabb_expressed+$RFaa_expressed+$RFdd_expressed;
		$RFploidy{$RFind_ref}=2 if ($RFaa_expressed==0 and $RFdd_expressed==0);
		if ($RFploidy{$RFind_ref}<1 or $RFploidy{$RFind_ref}>3) {
			print STDERR "SUB(RunFreebayes)Error: $RFind_ref un-expected ploidy value: $RFploidy{$RFind_ref}\n";
		}
		elsif ($RFploidy{$RFind_ref}==1) {
			print STDERR "SUB(RunFreebayes)Warnings: $RFind_ref expected ploidy value: $RFploidy{$RFind_ref}\n";
			next;
		}
		print STDOUT "Reference: $RFind_ref\tPloidy: $RFploidy{$RFind_ref}\n";
##COMMENT: runFreebayes and read vcf into hash
		my $RFcmd_merge_vcf='';
		$RFcmd_merge_vcf="$path_vcfmerge";
		my $RFcmd_freebayes_aabbdd='';
		$RFcmd_freebayes_aabbdd="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy{$RFind_ref} --pooled-discrete --min-mapping-quality $min_mapq --only-use-input-alleles --variant-input AABBDD.guide.vcf --region $RFind_ref  $RFbam_aabbdd | gzip -c > AABBDD/AABBDD.$RFind_ref.vcf.gz";
		if (&exec_cmd_return($RFcmd_freebayes_aabbdd)) {
			print STDERR "SUB(RunFreebayes)Error: Running freebayes for AABBDD error\n";
			return 1;
		}
		if (! -s "AABBDD/AABBDD.$RFind_ref.vcf.gz") {
			print STDERR "SUB(RunFreebayes)Error: freebayes output for AABBDD error\n";
			return 1;
		}else {
			$RFcmd_merge_vcf.=" AABBDD/AABBDD.$RFind_ref.vcf.gz";
		}
		
		my %RFaabb_freebayes_vcf=();
		if (${$RF_express[1]}{$RFind_ref}>0) {
			my $RFcmd_freebayes_aabb='';
			$RFcmd_freebayes_aabb="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy 2 --pooled-discrete --min-mapping-quality $min_mapq --only-use-input-alleles --variant-input AABBDD.guide.vcf --region $RFind_ref  $RFbam_aabb | gzip -c > AABB/AABB.$RFind_ref.vcf.gz";
			if (&exec_cmd_return($RFcmd_freebayes_aabb)) {
				print STDERR "SUB(RunFreebayes)Error: Running freebayes for AABB error\n";
				return 1;
			}
			if (! -s "AABB/AABB.$RFind_ref.vcf.gz") {
				print STDERR "SUB(RunFreebayes)Error: freebayes output for AABB error\n";
				return 1;
			}else {
				my $RFaabbvcf=&ReadVcf("AABB/AABB.$RFind_ref.vcf.gz");
				%RFaabb_freebayes_vcf=%{$RFaabbvcf};
				$RFcmd_merge_vcf.=" AABB/AABB.$RFind_ref.vcf.gz";
			}
		}
		
		my %RFaa_freebayes_vcf=();
		if (${$RF_express[2]}{$RFind_ref}>0) {
			my $RFcmd_freebayes_aa='';
			$RFcmd_freebayes_aa="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy 2 --pooled-discrete --min-mapping-quality $min_mapq --only-use-input-alleles --variant-input AABBDD.guide.vcf --region $RFind_ref  $RFbam_aa | gzip -c > AA/AA.$RFind_ref.vcf.gz";
			if (&exec_cmd_return($RFcmd_freebayes_aa)) {
				print STDERR "SUB(RunFreebayes)Error: Running freebayes for AA error\n";
				return 1;
			}
			if (! -s "AA/AA.$RFind_ref.vcf.gz") {
				print STDERR "SUB(RunFreebayes)Error: freebayes output for AA error\n";
				return 1;
			}else {
				my $RFaavcf=&ReadVcf("AA/AA.$RFind_ref.vcf.gz");
				%RFaa_freebayes_vcf=%{$RFaavcf};
				$RFcmd_merge_vcf.=" AA/AA.$RFind_ref.vcf.gz";
			}
		}
		
		my %RFdd_freebayes_vcf=();
		if (${$RF_express[3]}{$RFind_ref}>0) {
			my $RFcmd_freebayes_dd='';
			$RFcmd_freebayes_dd="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy 2 --pooled-discrete --min-mapping-quality $min_mapq --only-use-input-alleles --variant-input AABBDD.guide.vcf --region $RFind_ref  $RFbam_dd | gzip -c > DD/DD.$RFind_ref.vcf.gz";
			if (&exec_cmd_return($RFcmd_freebayes_dd)) {
				print STDERR "SUB(RunFreebayes)Error: Running freebayes for DD error\n";
				return 1;
			}
			if (! -s "DD/DD.$RFind_ref.vcf.gz") {
				print STDERR "SUB(RunFreebayes)Error: freebayes output for DD error\n";
				return 1;
			}else {
				my $RFddvcf=&ReadVcf("DD/DD.$RFind_ref.vcf.gz");
				%RFdd_freebayes_vcf=%{$RFddvcf};
				$RFcmd_merge_vcf.=" DD/DD.$RFind_ref.vcf.gz";
			}
		}
##COMMENT: merge vcf for debuging
		$RFcmd_merge_vcf.=" > Clust$cluster_linenum/$RFind_ref.merge.vcf";
		if (&exec_cmd_return($RFcmd_merge_vcf)) {
			print STDERR "SUB(RunFreebayes)Error: $RFind_ref vcf-merge running error\n";
			return 1;
		} elsif (! -s "Clust$cluster_linenum/$RFind_ref.merge.vcf") {
			print STDERR "SUB(RunFreebayes)Error: $RFind_ref vcf-merge output error\n";
			return 1;
		}
##COMMENT: Correct VCF
		
		
		
##COMMENT: output prephased vcf for hapcompass


##COMMENT: 




	}###End of freebayes call for each reference in one cluster
###delete############################################################
=delete
#Concatenate vcfs in one cluster
	my @RFaabbdd_vcfs=glob "AABBDD/AABBDD.*.vcf.gz";
	if (scalar(@RFaabbdd_vcfs)>1) {
		if (&exec_cmd_return("vcf-concat AABBDD/AABBDD.*.vcf.gz | gzip -c > Clust$cluster_linenum/AABBDD.vcf.gz")) {
			return 1;
		}
		elsif (-s "Clust$cluster_linenum/AABBDD.vcf.gz") {
			if (&exec_cmd_return("tabix -p vcf Clust$cluster_linenum/AABBDD.vcf.gz")) {
				return 1;
			}
			elsif (! -s "AABBDD.vcf.gz.tbi") {
				return 1;
			}
		}
		else {
			return 1;
		}
	}
	elsif (scalar(@RFaabbdd_vcfs)==1) {
		if (move($RFaabbdd_vcfs[0], "Clust$cluster_linenum/AABBDD.vcf.gz")) {
			return 1;
		}
	}else {
		return 1;
	}

	my %RFreadvcf_aabb=();
	if (${$RF_express[1]}{$RFind_ref}>0) {
		my @RFaabb_vcfs=glob "AABB/AABB.*.vcf.gz";
		if (scalar(@RFaabb_vcfs)>1) {
			if (&exec_cmd_return("vcf-concat AABB/AABB.*.vcf.gz | gzip -c > Clust$cluster_linenum/AABB.vcf.gz")) {
				return 1;
			}
			elsif (-s "Clust$cluster_linenum/AABB.vcf.gz") {
				if (&exec_cmd_return("tabix -p vcf Clust$cluster_linenum/AABB.vcf.gz")) {
					return 1;
				}
				elsif (! -s "AABB.vcf.gz.tbi") {
					return 1;
				} else {
					my $RFaabb_vcf=&readVcf("Clust$cluster_linenum/AABB.vcf.gz");
					%RFreadvcf_aabb=%{$RFaabb_vcf};
				}
			}
			else {
				return 1;
			}
		}
		elsif (scalar(@RFaabb_vcfs)==1) {
			if (move($RFaabb_vcfs[0], "Clust$cluster_linenum/AABB.vcf.gz")) {
				return 1;
			}
		}else {
			return 1;
		}
	}
	
	my %RFreadvcf_aa=();
	if (${$RF_express[2]}{$RFind_ref}>0) {
		my @RFaa_vcfs=glob "AA/AA.*.vcf.gz";
		if (scalar(@RFaa_vcfs)>1) {
			if (&exec_cmd_return("vcf-concat AA/AA.*.vcf.gz | gzip -c > Clust$cluster_linenum/AA.vcf.gz")) {
				return 1;
			}
			elsif (-s "Clust$cluster_linenum/AA.vcf.gz") {
				if (&exec_cmd_return("tabix -p vcf Clust$cluster_linenum/AA.vcf.gz")) {
					return 1;
				}
				elsif (! -s "AA.vcf.gz.tbi") {
					return 1;
				} else {
					my $RFaa_vcf=&readVcf("Clust$cluster_linenum/AA.vcf.gz");
					%RFreadvcf_aa=%{$RFaa_vcf};
				}
			}
			else {
				return 1;
			}
		}
		elsif (scalar(@RFaa_vcfs)==1) {
			if (move($RFaa_vcfs[0], "Clust$cluster_linenum/AA.vcf.gz")) {
				return 1;
			}
		}else {
			return 1;
		}
	}
	
	my %RFreadvcf_dd=();
	if (${$RF_express[3]}{$RFind_ref}>0) {
		my @RFdd_vcfs=glob "DD/DD.*.vcf.gz";
		if (scalar(@RFdd_vcfs)>1) {
			if (&exec_cmd_return("vcf-concat DD/DD.*.vcf.gz | gzip -c > Clust$cluster_linenum/DD.vcf.gz")) {
				return 1;
			}
			elsif (-s "Clust$cluster_linenum/DD.vcf.gz") {
				if (&exec_cmd_return("tabix -p vcf Clust$cluster_linenum/DD.vcf.gz")) {
					return 1;
				}
				elsif (! -s "DD.vcf.gz.tbi") {
					return 1;
				} else {
					my $RFdd_vcf=&readVcf("Clust$cluster_linenum/DD.vcf.gz");
					%RFreadvcf_dd=%{$RFdd_vcf};
				}
			}
			else {
				return 1;
			}
		}
		elsif (scalar(@RFdd_vcfs)==1) {
			if (move($RFdd_vcfs[0], "Clust$cluster_linenum/DD.vcf.gz")) {
				return 1;
			}
		}else {
			return 1;
		}
	}
###ReadVCF produced
	if (-s "Clust$cluster_linenum/AABBDD.vcf.gz") {
		if (open (AABBDD_ORIVCF, "zcat Clust$cluster_linenum/AABBDD.vcf.gz|")) {
			print STDERR "";
			return 1;
		}
	}
	else {
		return 1;
	}
	while (my $RFline=<AAAABBDD_ORIVCF>) {
		
	}
	close AABBDD_ORIVCF;
=cut
###Delete############################################################
}



###Asigning SNP allele 
###&AssignVariationAllele()
###Global:
###Dependancy: 
sub AssignVariationAllele {
	my ($AVAvcf_obj, $AVAtest_abd, $AVAtest_ab, $AVAtest_a, $AVAtestd)=@_;
	
	my %AVAvcf=%{$AVAvcf_obj};
	foreach my $AVAchrom (keys %AVAvcf) {
		foreach my $AVApos (keys %{$AVAvcf{$AVAchrom}}) {
			print "SUB(AssignVariationAllele)Test: Chr: $AVAchrom\tPos: $AVApos\t Arr: @{${$AVAvcf{$AVAchrom}}{$AVApos}}\n";
			
		}
	}
	 ###PAUSE
}



###Read vcf into hash, return index
###&ReadVcf(vcf_file)
###Global: $debug, 
###Dependancy: 
###VCF
#0		1	2	3	4	5		6		7		8		9
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ParLeaf
#0	ParMerge_MblContig65
#1	1106
#2	.
#3	C
#4	T
#5	310.79
#6	.
#7		AB=0.44;ABP=3.79203;AC=1;AF=0.333333;AN=3;AO=11;CIGAR=1X;DP=25;DPB=25;DPRA=0;EPP=19.0002;EPPR=5.49198;GTI=0;LEN=1;MEANALT=1;MQM=41.4545;MQMR=41.8571;NS=1;NUMALT=1;ODDS=4.15888;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=431;QR=517;RO=14;RPL=8;RPP=7.94546;RPPR=5.49198;RPR=3;RUN=1;SAF=4;SAP=4.78696;SAR=7;SRF=8;SRP=3.63072;SRR=6;TYPE=snp;technology.ILLUMINA=1	
#8	GT:DP:RO:QR:AO:QA
#9	0/0/1:25:14:517:11:431
sub ReadVcf {
	my $RVvcf=shift;
	if ($RVvcf=~/\.gz$|\.gzip/i) {
		open (VCF, "zcat $RVvcf|") || die "SUB(ReadVcf)Error: can not open vcf.gz file\n";
	}
	elsif ($RVvcf=~/\.vcf$/i) {
		open (VCF, "cat $RVvcf|") || die "SUB(ReadVcf)Error: can not open vcf file\n";
	}
	else {
		die "SUB(ReadVcf)Error: can not guess vcf file format\n";
	}
	my $num_line=0;
	my %vcf=();
	while (my $line=<VCF>) {
		$num_line++;
		next if ($line =~/^#/);
		chomp $line;
		my @arr=();
		@arr=split(/\t/, $line);
#next line if column number <=10
		if (scalar(@arr)<10) {
			print STDERR "SUB(ReadVcf)Error: Wrong variation at $num_line in $RVvcf\n";
			next;
		}
#check if allele number >1
		my ($RVchr, $RVpos, $RVref, $RVvar)=('', '', '', '');
		($RVchr, $RVpos, $RVref, $RVvar)=($arr[0], $arr[1], $arr[3], $arr[4]);
		if (exists ${$vcf{$RVchr}}{$RVpos} ) {
			print STDERR "SUB(ReadVcf)Error: exists $RVchr - $RVpos at $num_line in $RVvcf\n";
			next;
		}
		@{${$vcf{$RVchr}}{$RVpos}}=($RVref, $RVvar);
		for (my $i=9; $i<scalar(@arr); $i++) {
			my @arr2=();
			@arr2=split(/:/, $arr[$i]);
#comment following two lines if not ignore those sites with only one var allele
#			unless (&IsZeroIn($arr2[2])) {
#				delete ${$vcf{$RVchr}}{$RVpos};
#				last;
#			}
#			unless (&IsZeroIn($arr2[4])){
#				delete ${$vcf{$RVchr}}{$RVpos};
#				last;
#			}
			push (@{${$vcf{$RVchr}}{$RVpos}}, $arr2[0]);
		}
	}
	close VCF;
	
	#	print "vcf input\n";###test###
	if ($debug) {
		foreach my $chrom (keys %vcf) {
			my @arr3=keys %{$vcf{$chrom}};
			@arr3=sort {$a<=>$b} @arr3;
			print "SUB(ReadVcf)Reference: $chrom\n";
			print "SUB(ReadVcf)Number of variations on $chrom: ", scalar(@arr3), "\n";
			foreach (@arr3) {
				print $chrom, "\t", "$_\t${${$vcf{$chrom}}{$_}}[0]\t${${$vcf{$chrom}}{$_}}[1]\t${${$vcf{$chrom}}{$_}}[2]\n";
			}
		}
	}
	return \%vcf;
###return
###%vcf=(chr1 => (pos1 => (ref, var, gen, gen2, ...),
###				  pos2 => (ref, var, gen, gen2, ...),
###				...
###				),
###		chr2 => (pos1 => (ref, var, gen, gen2, ...),
###				  pos2 => (ref, var, gen, gen2, ...),
###				...
###				),
###		...
###		);
###$vcf{chr}->{pos}=(ref, var, gen, gen2, ...);
###Be noted: var could be comma-delimited string, in case of two variants
}



###RunHapCompass
###&RunHapCompass()
###Global: $path_java, $path_hapcompassjar, $path_hc2vcf
###Require: hapcompass.jar hc2vcf.jar
###Dependancy: freebayes
###Note: 
sub RunHapCompass {
	my ($RHCref, $RHCbam, $RHCvcf, $RHCploidy, $RHCoutput_vcf_prefix)=@_;
##COMMIT: check input files for HapCompass
	unless (-s $RHCref and -s $RHCbam and -s $RHCvcf) {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass inputs error: \n";
		return 1;
	}
	unless ($RHCploidy>=2) {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass ploidy error: \n";
		return 1;
	}
	if (! defined $RHCoutput_vcf_prefix or $RHCoutput_vcf_prefix eq '') {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass output prefix error: \n";
		return 1;
	}
##COMMIT: Run HapCompass, output: $RHCoutput_vcf_prefix
##									_frags.txt
##									_MWER_solution.txt
##									_phasedSolution.txt
##									_reads.sam
##									_reduced_representation.sam
##									_reduced_representation.vcf
	my $RHCcmd_hapcompass='';
	$RHCcmd_hapcompass="$path_java -jar $path_hapcompassjar --reference $RHCref --bam $RHCbam --vcf $RHCvcf --ploidy $RHCploidy --output $RHCoutput_vcf_prefix";
	if (&exec_cmd_return($RHCcmd_hapcompass)) {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass running error\n";
		return 1;
	}
##COMMIT: check HapCompass output
	if (! -s $RHCoutput_vcf_prefix.'_MWER_solution.txt') {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass output *_MWER_solution.txt error\n";
		return 1;
	}
##COMMIT: run hc2vcf to convert hapcompass output to vcf
##			output $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf'
	my $RHCcmd_hc2vcf='';
	$RHCcmd_hc2vcf="$path_java -jar $path_hc2vcf $RHCvcf $RHCploidy";
	if (&exec_cmd_return($RHCcmd_hc2vcf)) {
		print STDERR "SUB(RunHapCompass)Reference: hc2vcf running error\n";
		return 1;
	}
	if (! -s $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf') {
		print STDERR "SUB(RunHapCompass)Reference: hc2vcf running error\n";
		return 1;
	}
	else {
		return $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf';
	}
}



### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependancy:
###Note: 1=Bio::DB::Sam objective; 2=reference name array, 3=???
sub ReadSam {
	my ($RSbam_file, $RSref_file, $RSret_code)=@_;
	my $RSsam_obj=Bio::DB::Sam->new(  -bam => "$RSbam_file", 
							-fasta => "$RSref_file",
							-expand_flags => 1,
							-split_splices => 1,
							-autoindex => 1,
						);
	if ($RSret_code==1) {
		return $RSsam_obj;
	}
	elsif ($RSret_code==2) {
		my @arr=$RSsam_obj->seq_ids;
		return @arr;
	}
	elsif ($RSret_code==3) {
		##do sth
	}
	else {
		die "SUB(ReadSam)Error: unknown return code\n";
	}
}



###group SNPs based on read length
###GroupVcf($ReadVcf, $sam, reference_fa)
###Global: $debug
###Dependancy:
sub GroupVcf {
	my ($GVinput, $GVsamfile, $GVref)=@_;
	my %GVvcf=%{$GVinput};#ReadVcf object
	my %snp_phased_by_reads=();
	my $GVsam=&ReadSam($GVsamfile, $GVref, 1);
	foreach my $chrom (keys %GVvcf) {
		my @positions=sort {$a<=>$b} (keys %{$GVvcf{$chrom}});
		print "SUB(GroupVcf)Test: Reference2: $chrom\n" if ($debug);
		print "SUB(GroupVcf)Test: Number of variations on $chrom: ", scalar(@positions), "\n" if ($debug);
		my %pos=();
		foreach (@positions) {
			print $chrom, "\t", "Pos: $_\tRef:${${$GVvcf{$chrom}}{$_}}[0]\tVar${${$GVvcf{$chrom}}{$_}}[1]\tGen:${${$GVvcf{$chrom}}{$_}}[2]", "\n" if ($debug);
			@{$pos{$_}}= $GVsam->get_features_by_location(-seq_id => "$chrom", -start => "$_", -end => "$_");
		}
		print "SUB(GroupVcf)Chrom: $chrom\nPositions: @positions\n" if ($debug);
		foreach my $position (@positions) {
			foreach my $position2 (@positions) {
				next if ($position2 <= $position);
				my %share_alignments=();
				my $num_share_alignments=0; 
##Calculate number of shared reads at two vcf SNP/InDel sites (including same id for mate-pair)
				foreach my $align (@{$pos{$position}}) {
					foreach my $align2 (@{$pos{$position2}}) {
						if ($align->name eq $align2->name) {
							${$share_alignments{$align->name}}{$position}=$align;
							${$share_alignments{$align2->name}}{$position2}=$align2;
							$num_share_alignments++;
						}
					}
				}
				next if ($num_share_alignments <$min_share_alignments);#Ignore these two sites if number of reads less than expected
##Retrieve SNP/InDel SNP type for each shared reads
				my @share_readids=keys %share_alignments;
				print "Chr $chrom Pos $position vs Pos $position2\n";#." Shared IDs:\n@share_readids\n";
				my %GVtest=();
				my %alleleA=(); my %alleleB=();
				foreach my $readid (@share_readids) {
					my @readpos=keys %{$share_alignments{$readid}};
					@readpos=sort {$a <=> $b} @readpos;
					my $sam1=${$share_alignments{$readid}}{$position};
					my $align2geno=&ReadVariantType($position, ${${$GVvcf{$chrom}}{$position}}[0], ${${$GVvcf{$chrom}}{$position}}[1], ${${$GVvcf{$chrom}}{$position}}[2], $sam1);
					my $sam2=${$share_alignments{$readid}}{$position2};
					my $align2geno2=&ReadVariantType($position2, ${${$GVvcf{$chrom}}{$position2}}[0], ${${$GVvcf{$chrom}}{$position2}}[1], ${${$GVvcf{$chrom}}{$position2}}[2], $sam2);
					print "SUB(GroupVcf)Test: $readid\nposi: $position\t$position2\ngeno: $align2geno\t$align2geno2\n" if ($debug);
					if ($align2geno ne '?' and $align2geno2 ne '?') {
						$GVtest{"$align2geno~$align2geno2"}++;
						$alleleA{$align2geno}++;
						$alleleB{$align2geno2}++;
					}
				}
##Detect if any allele at these two sites are locked = appear at the same time in these reads
				foreach my $GVphased_type (keys %GVtest) {
#					print $GVphased_type."\n";###test###
					my $GVtest2=0;
					if ($GVtest{$GVphased_type} >=$min_share_alignments) {
						my ($alleleA, $alleleB)=split(/~/, $GVphased_type);
						BLOCK1: {foreach my $GVa (keys %alleleA) {
							foreach my $GVb (keys %alleleB) {
								if ((($GVa eq $alleleA and $GVb ne $alleleB ) or ($GVa ne $alleleA and $GVb eq $alleleB)) and (exists $GVtest{"$GVa~$GVb"})) {
									$GVtest2++;
									last BLOCK1;
								}
							}
						}}
						if ($GVtest2==0) {
						print "Read-phased SNP:\tChr: ".$chrom."\tPos: ".$position.'~'.$position2."\t Alleles: ".$GVphased_type."\tNum: ".$GVtest{$GVphased_type}."\n";
						@{${${snp_phased_by_reads{$chrom}}{$position}}{$position2}}=($alleleA, $alleleB);
						@{${${snp_phased_by_reads{$chrom}}{$position2}}{$position}}=($alleleB, $alleleA);
					}else {next;}
					}
				}
			}
		}
	}
	return \%snp_phased_by_reads;
###Format
###%snp_phased_by_reads=(	chr1 => (	pos1	=>	(	pos20	=>	(Allele_pos1, Allele_pos20),
###														pos21	=>	(Allele_pos1, Allele_pos21),
###														...;
###														)
###										pos2	=>	(	pos31	=>	(Allele_pos2, Allele_pos31),
###														pos32	=>	(Allele_pos2, Allele_pos32),
###														...,
###														)
###									...
###								)
###							chr2	=> (...)
###						)
###$snp_phased_by_reads->{chr}->{pos1}->{pos2}=(Allele_pos1, Allele_pos2);
}



###Read variant type based on the vcf genotypes
###ReadVariantType(Pos, RefAllele, VarAllele, GenoType, Bio::DB::Sam alignment)
###Global: $debug
###Dependancy: Bio::DB::Sam
sub ReadVariantType {
	my ($RVTpos, $RVTref, $RVTvar, $RVTgeno, $RVTsam_align)=@_;
#	print "SUB(ReadVariantType)Input: $RVTpos, $RVTref, $RVTvar, $RVTgeno\n";###test###
#	print "DeliverHash: ", $RVTsam_align->name."\t".$RVTsam_align->seq_id."\t".$RVTsam_align->query->dna."\n";###test###
	my $RVTquery_seq=$RVTsam_align->query->dna;
	my $RVTcigar=$RVTsam_align->cigar_str;
	my $RVTmdstring=$RVTsam_align->get_tag_values('MD');
	print "SUB(ReadVariantType)Test: CIGAR: $RVTcigar\t\tMDstring: $RVTmdstring\n" if ($debug);
	my $returnString='';
##COMMENT: Separate genotype, Assign number to each genotype: ref=0, var1=1, var2=2 ...
	my %RVTallele2geno=();
	$RVTallele2geno{$RVTref}=0;
#	print "SUB(ReadVariantType)Test: ref: $RVTref, $RVTallele2geno{$RVTref}\n";###test###
	my $RVTi=1;
	my @RVTvars=split(/,/, $RVTvar);
	foreach (@RVTvars) {
		$RVTallele2geno{$_}=$RVTi;
		$RVTi++;
	}
#	foreach (keys %RVTallele2geno) {print "SUB(ReadVariantType)Test: ref2: $_, $RVTallele2geno{$_}\n";} ###test###
##COMMENT: Determine SNP/InDel vcf position for a read
	my $RVTrefPos = $RVTsam_align->start;
	print "SUB(ReadVariantType)Test: Alignments starts: $RVTrefPos\n" if ($debug);
	my $RVTseqPos = 0;
	my $RVTcigarOperations = &SplitCigar($RVTsam_align->cigar_str);
	my ($RVTlast_is_insert, $RVTlast_is_deletion, $RVTlast_is_chop, $RVTlast_is_match)=(0, 0, 0, 0);
	foreach my $operation (@$RVTcigarOperations) {
		my $cig_length = $operation->[0];
		my $cig_op = $operation->[1];
		print "SUB(ReadVariantType)Test: Cigar: ($cig_length, $cig_op)\n" if ($debug);
		my ($RVTadd_ref, $RVTadd_seq)=(0, 0);
		($RVTlast_is_insert, $RVTlast_is_deletion, $RVTlast_is_chop, $RVTlast_is_match)=(0, 0, 0, 0);
		my $test_type=0;
		if ($cig_op =~ /^D$/) {
			$returnString .= "Deletion, $RVTrefPos, $cig_length\n";
			$RVTadd_ref=$cig_length;
			$test_type=3;
			$RVTlast_is_deletion=1;
		}
		elsif ($cig_op =~ /^I$/) {
			my $insertedBases = substr($RVTquery_seq, $RVTseqPos, $cig_length);
			$returnString .= "Insertion, $RVTrefPos, $insertedBases\n";
			$RVTadd_seq=$cig_length;
			$test_type=2;
			$RVTlast_is_insert=1;
		}
		elsif ($cig_op =~ /^M$/) {
			$RVTadd_ref= $cig_length;
			$RVTadd_seq= $cig_length;
			$test_type=1;
			$RVTlast_is_match=1;
		}
		elsif ($cig_op =~ /^S$/) {
			$RVTadd_seq=$cig_length;
			$test_type=4;
			$RVTlast_is_chop=1;
			## Don't increment refPos
		}
		elsif ($cig_op =~ /^N$/) {
			$RVTadd_ref= $cig_length;
			$test_type=5;
		}
		elsif ($cig_op =~ /^H$/) {
			$test_type=6;
		}
		elsif ($cig_op =~ /^P$/) {
			$test_type=7;
		}
		else {
			print STDERR "SUB(ReadVariantType)Error: unrecognized Cigar State: $cig_op (Cigar: ". $RVTsam_align->cigar_str. " at position $RVTpos of Reference sequence ".$RVTsam_align->seq_id."\n";
			return '?';
		}
		if (($RVTrefPos <= $RVTpos) and ($RVTpos <($RVTrefPos+$RVTadd_ref))) {
			$RVTseqPos+=abs($RVTpos-$RVTrefPos);
			last;
		}
		else {
			next;
		}
		$RVTrefPos+=$RVTadd_ref;
		$RVTseqPos+=$RVTadd_seq;
	}
	print $RVTsam_align->start." $RVTpos $RVTseqPos\n" if ($debug);
#	print substr($RVTquery_seq, $RVTseqPos, 1)."\n";###test###
##COMMENT: Determine this read belone to which allele type;

	my @returnarr=();
#	my @length_allele2geno=map {length($_)} keys %RVTallele2geno;
#	@length_allele2geno=sort {$b<=>$a} @length_allele2geno;



	foreach my $RVTind_allele (keys %RVTallele2geno) {
#		print "SUB(ReadVariantType)Test: allele2geno: $RVTind_allele, $RVTallele2geno{$RVTind_allele}\n";###test###
		if (! defined $RVTind_allele or $RVTind_allele eq '') {
			print STDERR "SUB(ReadVariantType)Error: empty ref/var values in ".$RVTsam_align->seq_id."\t".$RVTpos."\t".$RVTref."\t". $RVTvar."\t".$RVTgeno."\n";
			next;
		}
=old_algorithm
		if ($RVTind_allele eq substr($RVTquery_seq, $RVTseqPos, length($RVTind_allele))) {
			push (@returnarr, $RVTallele2geno{$RVTind_allele});
			print "SUB(ReadVariantType)Test: allele2geno: $RVTind_allele, $RVTallele2geno{$RVTind_allele}\n" if ($debug);
		}
=cut
		my $EVTsubstr_length=length($RVTind_allele);
		if ((length($RVTsam_align->query->dna)-$RVTseqPos+1)< length($RVTind_allele)) {#if reads not long enough to cover allele
			$EVTsubstr_length=length($RVTsam_align->query->dna)-$RVTseqPos+1;
		}
		elsif (length($RVTind_allele) == length ($RVTref)) { 
			if ($RVTind_allele eq substr($RVTquery_seq, $RVTseqPos, $EVTsubstr_length) and ($RVTlast_is_match==1)) {
				push (@returnarr, $RVTallele2geno{$RVTind_allele});
				print "SUB(ReadVariantType)Test: allele2geno: $RVTind_allele, $RVTallele2geno{$RVTind_allele}\n" if ($debug);
			}
		}
		elsif (($RVTlast_is_deletion==1) or ($RVTlast_is_insert==1) or ($RVTlast_is_chop==1)) {
			#There is a deletion in variant
			if ($RVTind_allele eq substr($RVTquery_seq, $RVTseqPos, $EVTsubstr_length)) {
				push (@returnarr, $RVTallele2geno{$RVTind_allele});
				print "SUB(ReadVariantType)Test: allele2geno2: $RVTind_allele, $RVTallele2geno{$RVTind_allele}\n" if ($debug);
			}
		}
	}
#	print "@returnarr\n";###test
	print "SUB(ReadVariantType)Test: Return: ", (scalar(@returnarr)==1) ? $returnarr[0]:'?', "\n" if ($debug);
	return (scalar(@returnarr)==1) ? $returnarr[0]:'?';	  
}



###split cigar into double array ((2, M), (5, D), ....)
###&SplitCigar(SamCigar)
###Global: None
###Dependancy: none
sub SplitCigar($) {
	my $SCcigar_string = shift;
	my @SCreturn_cigar_arr;
	my (@SCcigars) = ($SCcigar_string =~ /(\d*\w)/g);
	foreach (@SCcigars) {
		my @SCoperation = ($_ =~ /(\d+)(\w)/);
		push @SCreturn_cigar_arr, \@SCoperation;
	}
	return \@SCreturn_cigar_arr;
}



###Assign ABD allele for wheat
###&AsignABDallele($readvcf)
###Global: None
###Dependancy: none
sub AsignABDallele {
	print "###SUB(AsignABDallele)###\n";
	my $ASAvcf=shift;
	my %ASAreadvcf=%{$ASAvcf};
	foreach my $ASAchrom (keys %ASAreadvcf) {
		my @ASApositions=sort {$a<=>$b} (keys %{$ASAreadvcf{$ASAchrom}});
#		print "SUB(GroupVcf)Test: Reference3: $ASAchrom\n";
#		print "SUB(GroupVcf)Test: Number of variations on $ASAchrom: ", scalar(@ASApositions), "\n";
		foreach my $ASApos (@ASApositions) {
#			print "SUB(GroupVcf)Test: pos: $ASApos\n";
#			print "SUB(GroupVcf)Test: @{${$ASAreadvcf{$ASAchrom}}{$ASApos}}\n";
			my ($ASAref, $ASAvar, @ASAgenos)=@{${$ASAreadvcf{$ASAchrom}}{$ASApos}};
			####PAUSE
		}
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
###Dependancy: none
sub IsZeroIn {
	my $IZIstr=shift;
	my @IZIarr=split(/,/, $IZIstr);
	my $ISIzero=1;
	foreach (@IZIarr) {
		$ISIzero=0 if ($_ ==0);
	}
	return $ISIzero;
}



#&mytime()
sub mytime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;
}
###Process command
#&exec_cmd(cmd)
sub exec_cmd {
	my ($cmd) = @_;
	print &mytime()."CMD: $cmd\n";
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
#	if ($return_code == -1) {
#		print “failed to execute: $!\n”;
#	}
#	elsif ($return_code & 127) {
#		printf “child died with signal %d, %s coredump\n”, ($return_code & 127),  ($return_code & 128) ? ‘with’ : ‘without’;
#	}
#	else {
#		printf “child exited with value %d\n”, $return_code >> 8;
#	}
	if ($return_code) {
#		print "Error, cmd: $cmd died with ReturnCode $return_code\n";
		die "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}
sub exec_cmd_return {
	my ($cmd) = @_;
	print &mytime()."CMD: $cmd\n";
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
#	if ($return_code == -1) {
#		print “failed to execute: $!\n”;
#	}
#	elsif ($return_code & 127) {
#		printf “child died with signal %d, %s coredump\n”, ($return_code & 127),  ($return_code & 128) ? ‘with’ : ‘without’;
#	}
#	else {
#		printf “child exited with value %d\n”, $return_code >> 8;
#	}
	if ($return_code) {
		print STDERR "Error, cmd: $cmd died with ReturnCode $return_code\n";
#		die "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}



###cdbfasta index fasta, return fasta index name for cdbyank
###&CdbFasta(fasta_file)
###Global:
###Dependancy: &exec_cmd
sub CdbFasta {
	my $CSfasta_file=shift;
	if (-s "$CSfasta_file.cdbz" and -s "$CSfasta_file.cdbz.cidx") {
		return "$CSfasta_file.cdbz.cidx";
	}
	elsif (-s $CSfasta_file) {
		my $cmd="cdbfasta $CSfasta_file -z $CSfasta_file.cdbz";
		&exec_cmd($cmd);
		if (-s "$CSfasta_file.cdbz" and -s "$CSfasta_file.cdbz.cidx") {
			return "$CSfasta_file.cdbz.cidx";
		}
		else {
			die "SUB(CdbFasta)Error: creating fasta index failed\n";
		}
	}
	else {
		die "SUB(CdbFasta)Error: fasta not exist or empty\n";
	}
}



###Retrieve fasta sequences using input
###&CdbYank(cdbfasta_index, $CYoutput, @CYseq_ids);
###Global:
###Dependancy: &exec_cmd
sub CdbYank {
	my ($CYindex, $CYoutput, @CYseq_ids)=@_;
	my $CYseqids_join=join('\n', @CYseq_ids);
	my $cmd="echo -e \"$CYseqids_join\" | cdbyank $CYindex -o $CYoutput -w";
	if (&exec_cmd_return($cmd)) {
		return 1;
	}
	if (-s $CYoutput) {
		return 0;
	}
	else {
		print STDERR "SUB(CdbYank)Error: could not extract @CYseq_ids from $CYindex, and output $CYoutput\n";
		return 1;
		}
}



###Get filename without extension
###&RetrvNoExt(file)
###Golobal: none
###Dependancy: none
sub RetrvNoExt {
	my $RNEfile=shift @_;
	chomp $RNEfile;
	my ($RNEreturn, $RNEfilename)=('', '');
	($RNEfilename=$RNEfile)=~s/.*\///s;
	($RNEreturn=$RNEfilename)=~s/^(\S+)\.\w+$/$1/;
	die "SUB(RetrvNoExt)Error: empty file name\n" if (!defined $RNEreturn or ($RNEreturn eq ''));
	return $RNEreturn;
}




