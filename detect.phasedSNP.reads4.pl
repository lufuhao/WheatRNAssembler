#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use Bio::DB::Sam;
use Scalar::Util;
use File::Path;
use File::Copy 'move';
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150116

Requirements:
	Programs: hapcompass.jar hc2vcf.jar, gzip, gunzip, cat, zcat
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
our $bam_genome_tag='zw';
our $bam_chromo_tag='zc';
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
our $path_tabix='tabix';
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
our $path_vcfconcat='vcf-concat';


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
##COMMENT: Check if empty line
	if (scalar(@cluster_seqids)<1) {
		print STDERR "MainWarnings: line $cluster_linenum in $file_cluster ignored as empty\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
##COMMENT: extract fasta reference
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
##COMMENT: Empty folder
	unlink glob "$RunDir/AABBDD/*";
	unlink glob "$RunDir/AABB/*";
	unlink glob "$RunDir/AA/*";
	unlink glob "$RunDir/DD/*";
##COMMENT: Extract bam from @bam_AABBDD, @bam_AABB, @bam_AA, @bam_DD
	my @cluster_bam_files=();
	if (&ExtactBam(\@bam_AABBDD, \@cluster_seqids, "Clust$cluster_linenum/AABBDD.$cluster_linenum.bam")) {
		print STDERR "(Main)Error: AABBDD Bam extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	
	if (&ExtactBam(\@bam_AABB, \@cluster_seqids, "Clust$cluster_linenum/AABB.$cluster_linenum.bam")) {
		print STDERR "(Main)Error: AABB Bam extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	if (&ExtactBam(\@bam_AA, \@cluster_seqids, "Clust$cluster_linenum/AA.$cluster_linenum.bam")) {
		print STDERR "(Main)Error: AA Bam extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	if (&ExtactBam(\@bam_DD, \@cluster_seqids, "Clust$cluster_linenum/DD.$cluster_linenum.bam")) {
		print STDERR "(Main)Error: DD Bam extraction failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	@cluster_bam_files=("Clust$cluster_linenum/AABBDD.$cluster_linenum.bam", "Clust$cluster_linenum/AABB.$cluster_linenum.bam", "Clust$cluster_linenum/AA.$cluster_linenum.bam", "Clust$cluster_linenum/DD.$cluster_linenum.bam");
##COMMENT: run Express to evaluate FPKMs
	my ($test_express, $ref_expressed_index)=&ExpressFpkm("Clust$cluster_linenum/ref.$cluster_linenum.fa", $express_frag_len_mean, $express_frag_len_stddev, $express_max_read_len, \@cluster_bam_files);
	if ($test_express) {
		print STDERR "(Main)Error: express FPKM failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	my ($test_aabbdd_exporessed, $test_aabb_expressed, $test_aa_expressed, $test_dd_expressed)=@{$ref_expressed_index};
##COMMENT: Call SNP and prepahse
	if (&RunFreebayes("Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "Clust$cluster_linenum/AABBDD.$cluster_linenum.bam", 3, "Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz", 0)) {
		print STDERR "(Main)Error: freebayes running AABBDD failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	my $cmd="gunzip -c Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz > Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf";
	if (&exec_cmd_return("gunzip -c Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz > Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf")) {
		print STDERR "(Main)Error: gunzip uncompress Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
	elsif (! -s "Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf") {
		print STDERR "(Main)Error: gunzip output Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
#run freebayes for AABB
	if (&RunFreebayes("Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "Clust$cluster_linenum/AABB.$cluster_linenum.bam", 2, "Clust$cluster_linenum/AABB.$cluster_linenum.vcf.gz", "Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf")) {
		print STDERR "(Main)Error: freebayes running AABB failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
#run_freebayes for AA
		if (&RunFreebayes("Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "Clust$cluster_linenum/AA.$cluster_linenum.bam", 2, "Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz", "Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf")) {
		print STDERR "(Main)Error: freebayes running AA failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
#run_freebayes for DD
	if (&RunFreebayes("Clust$cluster_linenum/ref.$cluster_linenum.fa", \@cluster_seqids, "Clust$cluster_linenum/DD.$cluster_linenum.bam", 2, "Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz", "Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf")) {
		print STDERR "(Main)Error: freebayes running DD failed\n";
		print CLUSTERLOG $cluster_line."\tFail\t1\n";
		next;
	}
##COMMENT: Group alleles into subgenome
	my $aabbdd_bam_obj=&ReadSam("Clust$cluster_linenum/AABBDD.$cluster_linenum.bam", "Clust$cluster_linenum/ref.$cluster_linenum.fa", 1);
	my $aabbdd_vcf_obj=&ReadVcf("Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz");
	my $aabb_vcf_obj=&ReadVcf("Clust$cluster_linenum/AABBDD.$cluster_linenum.vcf.gz");
	my $aa_vcf_obj=&ReadVcf("Clust$cluster_linenum/AA.$cluster_linenum.vcf.gz");
	my $dd_vcf_obj=&ReadVcf("Clust$cluster_linenum/DD.$cluster_linenum.vcf.gz");
	my ($test_groupvcf, $fixed_allele_hashindex, $rnaseq_allele_hashindex, $allele2readids_hashindex, $readid2genome_hashindex, $sharedReadid_hashindex)=&GroupVcf($aabbdd_bam_obj, $aabbdd_vcf_obj, $aabb_vcf_obj, $aa_vcf_obj, $dd_vcf_obj);
	&FillInVariations($aabbdd_bam_obj, $aabbdd_vcf_obj, $aabb_vcf_obj, $aa_vcf_obj, $dd_vcf_obj, $fixed_allele_hashindex);
	if ($debug) {
	
	}


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



###extract batch of sequence alignment from Bams, and index
###&ExtactBam(bam_arr, sequence_arr, output)
###Global:$cluster_linenum, $path_samtools
###Dependency:
###Note: 
sub ExtactBam {
	my ($EBfiles_bam_index, $EBseq_obj, $EBoutput)=@_;
	my $EBbam_failure=0;
	return 1 if (scalar(@{$EBseq_obj})<1);###return error if no sequences
##COMMENT: extract bam from each
	my $EBsam_seq=join (' ', @{$EBseq_obj});
	my $EBi=0;
	my @EBtemp_bams=();
	print STDERR "SUB(ExtractBam)Info: sequence to be extracted: $EBsam_seq\n";###for test###
	foreach my $EBind_bam (@{$EBfiles_bam_index}) {
		$EBi++;
		my $EBcmd="$path_samtools view -bh $EBind_bam $EBsam_seq > $EBoutput.$EBi.bam";
		if (&exec_cmd_return($EBcmd)) {
			print STDERR "SUB(ExtractBam)Error: bam extract $EBsam_seq (Cluster: $cluster_linenum) from $EBind_bam failed\n";
			$EBbam_failure=1;
			return 1;
		}
		else {
			push (@EBtemp_bams, "$EBoutput.$EBi.bam");
		}
	}
##COMMENT: merge bams 
	my $EBmerge_bams=join(' ', @EBtemp_bams);
	my $EBcmd2="$path_samtools merge -r $EBoutput $EBmerge_bams";
	if (&exec_cmd_return($EBcmd2)) {
		print STDERR "SUB(ExtractBam)Error: bam merge (Cluster: $cluster_linenum) from $EBmerge_bams failed\n";
		return 1;
	}
	else {
		if (-s $EBoutput) {
			map {unlink($_)} @{$EBfiles_bam_index};###delete temporary files
		}
		else{
			print STDERR "SUB(ExtractBam)Error: bam merge (Cluster: $cluster_linenum) from $EBmerge_bams not found\n";
			return 1;
		}
	}
	my $EBcmd3="$path_samtools index $EBoutput";
	if (&exec_cmd_return($EBcmd3)) {
		print STDERR "SUB(ExtractBam)Error: bam index (Cluster: $cluster_linenum) from $EBmerge_bams failed\n";
		return 1;
	}
	return 0;
}



###ExpressRPKM
###&ExpressFpkm($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, @EFsamfiles)
###Global:$path_express, $cluster_linenum, @cluster_seqids
###Dependancy: &exec_cmd, &mytime, &RetrvNoExt
sub ExpressFpkm {
	my ($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, $EFsamfiles_index)=@_;
	my @EFbamfiles=@{$EFsamfiles_index};
	my @EFfpkms=();
	my $EFi=0;
	my $EFexpress_fail=0;
#Format: @return_arr=(AABBDD(1/0), AABB(1/0), AA(1/0), DD(1/0))
	my @return_arr=();
#@EFsamfiles=(bam_AABBDD, bam_AABB, bam_AA, bam_DD) in this order
	foreach my $EFbamfile (@EFbamfiles) {
		my $EFbamfile_base=&RetrvNoExt($EFbamfile);
		unless (-d "Express/$EFbamfile_base") {
			unless (mkdir ("Express/$EFbamfile_base", 0777)) {
				print STDERR "SUB(ExpressFpkm)Error: can not create folder $EFbamfile_base\n";
				return 1;
			}
		}
		unlink glob "$EFbamfile_base/*.xprs";###delete last-run files
		if (-s $EFbamfile) {
			my $EFcmd="$path_express --frag-len-mean $EFfrag_len_mean --frag-len-stddev $EFfrag_len_stddev --max-read-len $EFmax_read_len --output-dir Express/$EFbamfile_base $EFref $EFbamfile";
			if (&exec_cmd_return($EFcmd)) {
				return 1;
			}
			elsif (! -s "Express/$EFbamfile_base/results.xprs") {
				return 1;
			}
			unless (open (FPKM, "Express/$EFbamfile_base/results.xprs")) {
				print STDERR "SUB(ExpressFpkm)Error: can not open file Express/$EFbamfile_base/results.xprs\n";
				return 1;
			}
			my $EFline=<FPKM>;###firstline is header
			my $EFtest_expressed=0;
			while ($EFline=<FPKM>) {
				chomp $EFline;
				my @arr=split(/\t/, $EFline);
#$arr[1] is reference sequence ID
#$arr[10] is the estimated relative abundance in units of fragments per kilobase per million mapped
				if ($arr[10]>0) {###express if any reference FPKM >0 in this cluster
					$EFtest_expressed=1;
				}
				${$EFfpkms[$EFi]}{$arr[1]}=$arr[10];
			}
			close FPKM;
			push (@return_arr, $EFtest_expressed);
			remove_tree("Express/$EFbamfile_base");#delete last-run files
			#unlink glob "$EFbamfile_base/*.xprs";###delete last-run files
			####Delete directory###PAUSE
			$EFi++;
		}
		else {
			print STDERR "SUB(ExpressFpkm)Error: bam file $EFbamfile empty or non-exists\n";
			return 1;
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
			print "SUB(ExpressFpkm)Test: Bam $EFbamfiles[$i]: Ref $EFseqid: ${$EFfpkms[$i]}{$EFseqid}\n" if ($debug); ###Test###
		}
		print FPKMLOG $EFfpkm_logprint."\n";
	}
	#return \@EFfpkms;
#Return format
#@EFfpkms=(
#			(ref1 => fpkm1, ref2 => fpkm2, ...), 
#			(ref1 => fpkm1, ref2 => fpkm2, ...), 
#			...)
#$EFfpkms[$ith]->{reference_id}=fpkm
	unless (scalar(@EFbamfiles) == scalar(@return_arr)) {
		print 
		return 1;
	}
	return (0, \@return_arr);
}



###RunFreebayes, return merged vcf and raw AABBDD vcf
###&RunFreebayes(ref_fasta, seqids_arrindex, bam, ploidy, output.vcf.gz, guided_vcf)
###Global: $freebayes_min_coverage, $freebayes_min_alternative_count, $min_mapq,$express_frag_len_mean, $express_frag_len_stddev, 
###Dependancy: 
sub RunFreebayes {
	my ($RFfile_reference, $RFseqids_arrindex, $RFfile_bam, $RFploidy, $RFoutput, $RFguide_vcf)=@_;
	my $RFcmd='';
	my @RFseqids=@{$RFseqids_arrindex};
	if (scalar(@RFseqids)<1) {
		print "SUB(RunFreebayes)Error: no seqs for freebayes";
		return 1;
	}
	if ($RFguide_vcf ==0) {
		$RFcmd="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy --pooled-discrete --min-mapping-quality $min_mapq $RFfile_bam | gzip -c > $RFoutput";
		if (&exec_cmd_return($RFcmd)) {
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
			$RFcmd="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy --pooled-discrete --min-mapping-quality $min_mapq --only-use-input-alleles --variant-input $RFguide_vcf --region $RFseqids[0]  $RFfile_bam | gzip -c > $RFoutput";
			if (&exec_cmd_return($RFcmd)) {
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
				$RFcmd="$path_freebayes --fasta-reference $RFfile_reference --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --ploidy $RFploidy --pooled-discrete --min-mapping-quality $min_mapq --only-use-input-alleles --variant-input $RFguide_vcf --region $RFind_seq  $RFfile_bam | gzip -c > Freebayes/$RFind_seq.vcf.gz";
				if (&exec_cmd_return($RFcmd)) {
					print STDERR "SUB(RunFreebayes)Error3:  freebayes running error\n";
					return 1;
				}
				if (! -s "Freebayes/$RFind_seq.vcf.gz") {
					print STDERR "SUB(RunFreebayes)Error3: freebayes output error\n";
					return 1;
				}
				push (@RFtemp_vcf, "Freebayes/$RFind_seq.vcf.gz");
			}
			$RFcmd="$path_vcfconcat ".join(' ', @RFtemp_vcf)." | gzip -c > $RFoutput";
			if (&exec_cmd_return($RFcmd)) {
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
	my $RFcmd2="$path_tabix -p vcf $RFoutput";
	if (&exec_cmd_return($RFcmd2)) {
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
	my %AVAgenome2allele=();
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
	return (\%AVAgenome2allele);
}



###Extract allele based VCF ref, var, geno
###&ExtractAllele(ref, var, geno)
###Example: &ExtractAllele('AA', 'TT,GG', 0/1/1)
###Global: 
###Dependency:
###Note: 
sub ExtractAllele {
	my ($EAref, $EAvar, $EAgeno)=@_;
#Format: %EAgeno_hash=(0 => 'AA', 1 => 'TT', 2 => 'GG')
	my %EAgeno_hash=();
	$EAgeno_hash{0}=$EAref;
	my @EAvar_arr=split(/$geno_delimiter/,$EAvar);
	for (my $EAi=0; $EAi<scalar(@EAvar_arr);$EAi++) {
		$EAgeno_hash{$EAi+1}=$EAvar_arr[$EAi];
	}
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
###%vcf=(chr1 => (pos1 => (ref, var, gen/gen2, ...),
###				  pos2 => (ref, var, gen/gen2, ...),
###				...
###				),
###		chr2 => (pos1 => (ref, var, gen/gen2, ...),
###				  pos2 => (ref, var, gen/gen2, ...),
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
###&GroupVcf($ReadSam_obj, AABBDD_Vcf_obj, AABB_vcf_obj, AA_vcf_obj, DD_vcf_obj)
###Global: $debug
###Dependancy:&ReadVcf, &ReadVariantType, $min_share_alignments, $min_mapq, $bam_genome_tag
sub GroupVcf {
	my ($GVaabbdd_samobj, $GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj)=@_;
	my $GVfixAlleles_index=&AssignVariationAllele($GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj);
	my $GVtest_failed=0;###Running Control ###PAUSE
	my %GVvcf=%{$GVaabbdd_vcfobj};#ReadVcf object
#Format: %GVfixAlleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %GVfixAlleles_rnaseq=%{$GVfixAlleles_index};
#Format: %GVshare_readIDs = (chr => (readid1 => 1, readid2 => 1)), and then exclude those IDs assigned to alleles
	my %GVshare_readIDs=();
#Format: %GVreadid_at_pos=(chr => (pos => @readids)); All reads at each position
	my %GVreadid_at_pos=();
#Format: %GVreadid_by_allele=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
	my %GVreadid_by_allele=();
#Format: %GVreadid2genome=(readID=> ('A' =>1, 'B' =>1, 'C'=1))
	my %GVreadid2genome=();
#Format: %GVreadid2chrmosome=(readID => ('1AL' => 1; 3B => 1, ...)
	my %GVreadid2chrmosome=();

	foreach my $GVchrom (keys %GVvcf) {
		my @positions=sort {$a<=>$b} (keys %{$GVvcf{$GVchrom}});
		print "SUB(GroupVcf)Test: Reference2: $GVchrom\n" if ($debug);
		print "SUB(GroupVcf)Test: Number of variations on $GVchrom: ", scalar(@positions), "\n" if ($debug);
###COMMENT: retrieve all alignments mapped to this reference sequence
		my @GVchr_alignments=$GVaabbdd_samobj->get_features_by_location(-seq_id => "$GVchrom");
		foreach (@GVchr_alignments) {
			if ($_->qual>= $min_mapq) {#filter mapping quality
				${$GVshare_readIDs{$GVchrom}}{$_->name}++;
			}
		}
###COMMENT: retrieve alignment at each position, uniform read ID assigned to each subgenome
		my %GVpos=();
		foreach my $GVind_pos (@positions) {
			print "SUB(GroupVcf)Test:\t".$GVchrom, "\t", "Pos: $GVind_pos\tRef:${${$GVvcf{$GVchrom}}{$GVind_pos}}[0]\tVar:${${$GVvcf{$GVchrom}}{$GVind_pos}}[1]\tGen:${${$GVvcf{$GVchrom}}{$GVind_pos}}[2]", "\n" if ($debug);
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
							${${GVreadid2chrmosome{$GVind_align->name}}}{$1}++;
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
						print "BestCount: $GVbest_count\n";###for test###
					}
					elsif (${$GVallele_assign{$GVsubgenome}}{$GValleleat_pos} == $GVmax){
						$GVbest_count++;
					}
				}
				$GVbest_allele='?' if ($GVbest_count != 1);
				print $GVsubgenome."\t".$GVbest_allele."\n";
				if ($GVsubgenome eq 'A' and $GVaa_allele ne '?') {
					${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'A'}=$GVbest_allele;
				}
				elsif ($GVsubgenome eq 'B' and $GVbb_allele ne '?') {
					${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'B'}=$GVbest_allele;
				}
				elsif ($GVsubgenome eq 'B' and $GVdd_allele ne '?') {
					${${$GVfixAlleles_rnaseq{$GVchrom}}{$GVind_pos}}{'D'}=$GVbest_allele;
				}
			}
		}
		print "SUB(GroupVcf)Test: number of share reads ids on $GVchrom: ".scalar(keys %{$GVshare_readIDs{$GVchrom}})."\n" if ($debug);
	}
	return ($GVtest_failed, \%GVfixAlleles_rnaseq, $GVfixAlleles_index, \%GVreadid_by_allele, \%GVreadid2chrmosome, \%GVshare_readIDs);
}



###Read variant type based on the vcf genotypes
###ReadVariantType(Pos, RefAllele, VarAllele, GenoType, Bio::DB::Sam alignment)
###Global: $debug
###Dependancy: Bio::DB::Sam, $geno_delimiter
sub ReadVariantType {
	my ($RVTpos, $RVTref, $RVTvar, $RVTgeno, $RVTsam_align)=@_;
#	print "SUB(ReadVariantType)Input: $RVTpos, $RVTref, $RVTvar, $RVTgeno\n";###test###
#	print "DeliverHash: ", $RVTsam_align->name."\t".$RVTsam_align->seq_id."\t".$RVTsam_align->query->dna."\n";###test###
	my $RVTref_seqid=$RVTsam_align->seq_id;
	my $RVTreadid=$RVTsam_align->name;
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
	print $RVTsam_align->start." $RVTpos $RVTseqPos\n" if ($debug);
#	print "Ref: $RVTref_seqid\tPos: $RVTpos\tRef:$RVTref\tVar: $RVTvar\tReads: $RVTreadid\tLength: ".length($RVTquery_seq)."\t $EVTcapture_start_pos - $RVTcapture_end_pos\n";###For test###

	if ($RVTcapture_end_pos>length($RVTquery_seq)) {
		print STDERR "Readlength: ".length($RVTquery_seq)."\tSubstrEnd: $RVTcapture_end_pos\n";
		print STDERR "SUB(ReadVariantType)Error: use read ($RVTreadid) length at $RVTpos instead\n";###For test###
		$RVTcapture_end_pos=length($RVTquery_seq);
	}
	unless ($EVTcapture_start_pos < length($RVTquery_seq) and $EVTcapture_start_pos < $RVTcapture_end_pos) {
		print STDERR "SUB(ReadVariantType)Error: Unknown2\n";###For test###
		return "unknown";
	}
	my $RVTcapture_string=substr($RVTquery_seq, $EVTcapture_start_pos, $RVTcapture_end_pos-$EVTcapture_start_pos);
	unless (defined $RVTcapture_string and $RVTcapture_string ne '') {
		print STDERR "SUB(ReadVariantType)Error: substr failed\n";###For test###
		return "unknown2";
	}
#	print "Captured String: $RVTcapture_string\n";# if ($debug);
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






### FillInVariations
###&FillInVariations()
###Global:
###Dependency:
###Note
sub FillInVariations {
	##PAUSE
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
	print "#####\n".&mytime()."CMD: $cmd\n";
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
		die "SUB(exec_cmd)Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print STDERR "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}
sub exec_cmd_return {
	my ($cmd) = @_;
	print "#####\n".&mytime()."CMD: $cmd\n";
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
		print STDERR "SUB(exec_cmd_return)Error, cmd: $cmd died with ReturnCode $return_code\n";
#		die "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print STDERR "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
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
###&CdbYank(cdbfasta_index, $CYoutput, $CYseq_ids_arr_index);
###Global:
###Dependancy: &exec_cmd_return
sub CdbYank {
	my ($CYindex, $CYoutput, $CYseq_ids_index)=@_;
	my $CYseqids_join=join('\n', @{$CYseq_ids_index});
	my $CYcmd="echo -e \"$CYseqids_join\" | cdbyank $CYindex -o $CYoutput -w";
	if (&exec_cmd_return($CYcmd)) {
		return 1;
	}
	if (-s $CYoutput) {
		if (&exec_cmd_return("$path_samtools faidx $CYoutput")) {
			print STDERR "SUB(CdbYank)Error: could not index $CYseqids_join from $CYindex, and output $CYoutput\n";
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		print STDERR "SUB(CdbYank)Error: could not extract $CYseqids_join from $CYindex, and output $CYoutput\n";
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




