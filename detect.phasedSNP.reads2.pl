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
Version: LUFUHAO20150116

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
our $geno_delimiter=',';
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
	my $RFaabbdd_bam_obj=&ReadSam($RFbam_aabbdd, $RFfile_reference, 1);
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
		}
		else {
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
			}
			else {
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
			}
			else {
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
		}
		elsif (! -s "Clust$cluster_linenum/$RFind_ref.merge.vcf") {
			print STDERR "SUB(RunFreebayes)Error: $RFind_ref vcf-merge output error\n";
			return 1;
		}
##COMMENT: assign alleles to A B D
		my $RFaabbdd_vcf_obj=&ReadVcf("AABBDD/AABBDD.$RFind_ref.vcf.gz");
		my $RFaabb_vcf_obj=&ReadVcf("AABB/AABB.$RFind_ref.vcf.gz");
		my $RFaa_vcf_obj=&ReadVcf("AA/AA.$RFind_ref.vcf.gz");
		my $RFdd_vcf_obj=&ReadVcf("DD/DD.$RFind_ref.vcf.gz");
		my ($RFtest_assign_true, $test_b_assign, $RFassigned_allele2gnm, $RFassign_gnm2allele)=&AssignVariationAllele($RFaabbdd_vcf_obj, $RFaabb_vcf_obj, $RFaa_vcf_obj, $RFdd_vcf_obj);
		if ($RFtest_assign_true==1) {
			return 1
		}
		&GroupVcf($RFaabbdd_bam_obj, $RFaabbdd_vcf_obj, $RFaabb_vcf_obj, $RFaa_vcf_obj, $RFdd_vcf_obj);
##COMMENT: output prephased vcf for hapcompass


##COMMENT: SNP fix




	}###End of freebayes call for each reference in one cluster
}



###Asigning SNP allele 
###&AssignVariationAllele()
###Global:
###Dependancy: &ExtractAllele
sub AssignVariationAllele {
	my ($AVAaabbdd_vcf_index, $AVAaabb_vcf_index, $AVAaa_vcf_index, $AVAdd_vcf_index)=@_;
	my %AVAaabbdd_vcf_obj=%{$AVAaabbdd_vcf_index};
	my %AVAaabb_vcf_obj=%{$AVAaabb_vcf_index};
	my %AVAaa_vcf_obj=%{$AVAaa_vcf_index};
	my %AVAdd_vcf_obj=%{$AVAdd_vcf_index};
	my $AVAtest_AA_expressed=1;
	$AVAtest_AA_expressed=0 if (scalar(keys %{$AVAaa_vcf_index})==0);
	my $AVAtest_DD_expressed=1;
	$AVAtest_DD_expressed=0 if (scalar(keys %{$AVAdd_vcf_index})==0);
	my $AVAtest_bb_expressed=0;
#Format: %AVAgenome2allele=(chr => (pos => ('A' => allele, 'B' => Allele, 'D' => allele)))
	my %AVAgenome2allele==();
#Format: %AVAallele2genome=(chr => (pos => (allele => 'A', allele => 'B', allele => 'D')))
	my %AVAallele2genome=();
	foreach my $AVAchrom (keys %AVAaabbdd_vcf_obj) {
		foreach my $AVApos (sort {$a<=>$b} keys %{$AVAaabbdd_vcf_obj{$AVAchrom}}) {
			print "SUB(AssignVariationAllele)Test: Chr: $AVAchrom\tPos: $AVApos\t Arr: @{${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}\n";
			my ($AVAaabbdd_genohash_index, $AVAaabbdd_genoarr_index)=&ExtractAllele(${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAaabbdd_vcf_obj{$AVAchrom}}{$AVApos}}[2]);
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
					my ($AVAdd_genohash_index, $AVAdd_genoarr_index)=&ExtractAllele(${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAdd_vcf_obj{$AVAchrom}}{$AVApos}}[2]);
					if (scalar(@{$AVAdd_genoarr_index})==1 and ${$AVAdd_genoarr_index}[0]=~/^\d+$/) {
						$AVAdd_allele=${$AVAdd_genoarr_index}[0];
						if (exists ${$AVAaabbdd_genohash_index}{$AVAdd_allele} and ${$AVAdd_genohash_index}{$AVAdd_allele} eq ${$AVAaabbdd_genohash_index}{$AVAdd_allele}) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}=$AVAdd_allele;
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{$AVAdd_allele}.='D';
						}
						else {
							print STDERR "SUB(AssignVariationAllele)Warning: non AABBDD allele $AVAdd_allele in DD\n";
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'});
##COMMENT: defined AA allele
				if (exists ${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}) {
					my ($AVAaa_genohash_index, $AVAaa_genoarr_index)=&ExtractAllele(${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAaa_vcf_obj{$AVAchrom}}{$AVApos}}[2]);
					if (scalar(@$AVAaa_genoarr_index)==1 and $AVAaa_allele=~/^\d+$/) {
						$AVAaa_allele=${$AVAaa_genoarr_index}[0];
						if (exists ${$AVAaabbdd_genohash_index}{$AVAaa_allele} and ${$AVAaa_genohash_index}{$AVAaa_allele} eq ${$AVAaabbdd_genohash_index}{$AVAaa_allele}) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}=$AVAaa_allele;
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{$AVAaa_allele}.='A';
						}
						else {
								print STDERR "SUB(AssignVariationAllele)Warning: non AABBDD allele $AVAaa_allele in AA\n";
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'});
##COMMENT: defined B allele---Difficult
				if ($AVAtest_AA_expressed==1 and $AVAaa_allele=~/^\d$/ and exists ${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}) {
					my ($AVAaabb_genohash_index, $AVAaabb_genoarr_index)=&ExtractAllele(${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[0], ${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[1], ${${$AVAaabb_vcf_obj{$AVAchrom}}{$AVApos}}[2]);
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
				if (1) {
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
	return (\%AVAgenome2allele);
}



sub ExtractAllele {
	my ($EAref, $EAvar, $EAgeno)=@_;
	my %EAgeno_hash=();
	$EAgeno_hash{0}=$EAref;
	my @EAvar_arr=split(/$geno_delimiter/,$EAvar);
	for (my $EAi=0; $EAi<scalar(@EAvar_arr);$EAi++) {
		$EAgeno_hash{$EAi+1}=$EAvar_arr[$EAi];
	}
	my %EAgeno2_hash=();
	my @EAgeno_arr=();
	if ($EAgeno eq '.') {
		push (@EAgeno_arr, '?');
	}
	elsif ($EAgeno =~m/\//) {
		@EAgeno_arr=split(/\//, $EAgeno);
		foreach (@EAgeno_arr) {
			if ($_=~m/^\d$/) {
				$EAgeno2_hash{_}++;
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
	return (\%EAgeno_hash, \@EAgeno2_arr);
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
		foreach my $RVchrom (keys %vcf) {
			my @arr3=keys %{$vcf{$RVchrom}};
			@arr3=sort {$a<=>$b} @arr3;
			print "SUB(ReadVcf)Reference: $RVchrom\n";
			print "SUB(ReadVcf)Number of variations on $RVchrom: ", scalar(@arr3), "\n";
			foreach (@arr3) {
				print $RVchrom, "\t", "$_\t${${$vcf{$RVchrom}}{$_}}[0]\t${${$vcf{$RVchrom}}{$_}}[1]\t${${$vcf{$RVchrom}}{$_}}[2]\n";
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
###GroupVcf($ReadVcf_obj, $sam, reference_fa)
###Global: $debug
###Dependancy: &ReadVcf, &ReadVariantType, &ReadSam, $min_share_alignments
sub GroupVcf {
	my ($GVinput, $GVsamfile, $GVref)=@_;
	my %GVvcf=%{$GVinput};#ReadVcf object
	my %snp_phased_by_reads=();
	my $GVsam=&ReadSam($GVsamfile, $GVref, 1);
	foreach my $chrom (keys %GVvcf) {
		my @positions=sort {$a<=>$b} (keys %{$GVvcf{$chrom}});
		print "SUB(GroupVcf)Test: Reference2: $chrom\n" if ($debug);
		print "SUB(GroupVcf)Test: Number of variations on $chrom: ", scalar(@positions), "\n" if ($debug);
###COMMENT: retrieve alignment at each position
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
##COMMENT:Calculate number of shared reads at two vcf SNP/InDel sites (including same id for mate-pair)
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
##COMMENT: Retrieve SNP/InDel SNP type for each shared reads
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
##COMMENT: Detect if any allele at these two sites are locked = appear at the same time in these reads
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
			###do magic
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
###Dependancy: Bio::DB::Sam, $geno_delimiter
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
	my @RVTvars=split(/$geno_delimiter/, $RVTvar);
	foreach (@RVTvars) {
		$RVTallele2geno{$_}=$RVTi;
		$RVTi++;
	}
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
	my $RVTrefPos = $RVTsam_align->start;
	print "SUB(ReadVariantType)Test: Alignments starts: $RVTrefPos\n" if ($debug);
	my $RVTseqPos = 0;
	my $RVTthis_read_var='';
	my $EVTcapture_start_pos=0;
	my $RVTcapture_end_pos=0;
	my $RVTtest_var_start=0;
	my $RVTtest_var_end=0;
	my $RVTref_var_end=$RVTpos+length($RVTref);
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
		elsif ($cig_op =~ /^[SH]{1}$/) {
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
			if ($RVTcigar_type=~/^[56]{1}$/) {
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
	print $RVTsam_align->start." $RVTpos $RVTseqPos\n" if ($debug);
	my $RVTcapture_string=substr($RVTquery_seq, $EVTcapture_start_pos, $RVTcapture_end_pos-$EVTcapture_start_pos);
	print "Captured String: $RVTcapture_string\n" if ($debug);
	my @returnarr=();
	if ($RVTtest_var_end) {
		if (exists $RVTallele2geno{$RVTcapture_string}) {
			return $RVTallele2geno{$RVTcapture_string};
		}
		else {
			return 'NotExist';
		}
	}
	else {
		foreach my $RVTind_allele (keys %RVTallele2geno) {
			if ($RVTind_allele=~/^$RVTcapture_string/) {
				push (@returnarr, $RVTallele2geno{$RVTind_allele});
			}
		}
	}
#	print "@returnarr\n";###test
	return (scalar(@returnarr)==1) ? $returnarr[0]:'M';	###Modify  
}



###split cigar into double array ((2, M), (5, D), ....)
###&SplitCigar(SamCigar)
###Global: None
###Dependancy: none
sub SplitCigar($) {
	my $SCcigar_string = shift;
	my @SCreturn_cigar_arr=();
	while ($SCcigar_string =~ /(\d+)([MIDNSHP=X])/g) {
		my @SCoperation=($1, $2);
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
###Dependancy: $geno_delimiter
sub IsZeroIn {
	my $IZIstr=shift;
	my @IZIarr=split(/$geno_delimiter/, $IZIstr);
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




