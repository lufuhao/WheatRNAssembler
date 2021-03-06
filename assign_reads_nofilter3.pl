#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Sam;
use constant USAGE =><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150223

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--fasta	<Fasta>
		[Msg] Reference sequences in fasta
	--bam	<Bam>
		[Msg] Sorted Bam files
	--vcf	<VCF>
		[Msg] VCF file (Freebayes)
	--output|-o	<Include>
		[Opt] Output for included read ID, 
		Default: Bam_basename.include
	--exclude|-e	<Exclude>
		[Opt] Output for excluded Read ID,
		Default: Bam_basename.exclude
	--notsure|-n
		[Opt] Output for not sure (included and excluded) read ID
		Default: Bam_basename.notsure
	--minmapq|-q	<Int>
		[Opt] Minimum mapping quality
	--debug
		[Opt] Output detailed info for debugging
	--verbose
		Detailed output for trouble-shooting;
	--version|v!
		Print current SCRIPT version;

Example:
	perl $0 --fasta ref.fa --bam my.bam --vcf my.vcf -q 2

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
our ($help, $verbose, $debug, $ver);
our ($file_fasta, $file_bam, $file_vcf, $file_include, $file_exclude, $file_notsure);
our ($min_share_alignments, $min_mapq);
GetOptions(
	"help|h!" => \$help,
	"fasta:s" => \$file_fasta,
	"bam:s" => \$file_bam,
	"vcf:s" => \$file_vcf,
	"output|o:s" => \$file_include,
	"exclude|e:s" => \$file_exclude,
	"notsure|n:s" => \$file_notsure,
#	"minnumalign:i" => \$min_share_alignments,
	"minmapq|q:i" => \$min_mapq,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
$verbose=0 unless (defined $verbose);
$debug=0 unless (defined $debug);
#$min_share_alignments=3 unless (defined $min_share_alignments);
$min_mapq=0 unless (defined $min_mapq);
our $geno_delimiter=',';



### input and output ################################################
die USAGE unless (defined $file_fasta and -s $file_fasta and defined $file_bam and -s $file_bam and defined $file_vcf and $file_vcf);
my $output_prefix=&RetrvNoExt($file_bam);
$file_include=$output_prefix.".include" unless (defined $file_include);
$file_exclude=$output_prefix.".exclude" unless (defined $file_exclude);
$file_notsure=$output_prefix.".notsure" unless (defined $file_notsure);
unlink ($file_include) if (-e $file_include);
unlink ($file_exclude) if (-e $file_exclude);
unlink ($file_notsure) if (-e $file_notsure);

print "##### Summary #####\n";
print "Fasta: $file_fasta\nBam: $file_bam\nVcf: $file_vcf\nMinMAPQ: $min_mapq\n";
print "#Output: \nInclude: $file_include\nExclude: $file_exclude\nNotSure: $file_notsure\n";
print "#Extra: \nDebug: $debug\nVerbose: $verbose\n";
print "##########\n";




### Main ############################################################

#&ReadVcf($vcffile);

open (INCLUDE, ">$file_include") || die "MainError: write to $file_include\n";
open (EXCLUDE, ">$file_exclude") || die "MainError: write to $file_exclude\n";
&GetRefAlleleReads($file_vcf, $file_bam, $file_fasta);
close INCLUDE;
close EXCLUDE;



#####################################################################
###                         sub functions                         ###
#####################################################################
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
		foreach my $chrom (sort (keys %vcf)) {
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



###group SNPs based on read length
###GetRefAlleleReads($ReadVcf, $sam, reference_fa)
###Global: $debug, %include, %exclude, INCLUDE, EXCLUDE
###Dependancy:
sub GetRefAlleleReads {
	my ($GRARinput, $GRARsamfile, $GRARref)=@_;
	print "SUB(GetRefAlleleReads)Info: starting...\n" if ($debug);
	my %GRARvcf=%{&ReadVcf($GRARinput)};#ReadVcf object
	my $GRARsam=&ReadSam($GRARsamfile, $GRARref, 1);
	foreach my $GRARchrom (keys %GRARvcf) {
		my @GRARpositions=sort {$a<=>$b} (keys %{$GRARvcf{$GRARchrom}});
		print "SUB(GetRefAlleleReads)Test: Reference2: $GRARchrom\n" if ($debug);
		print "SUB(GetRefAlleleReads)Test: Number of variations on $GRARchrom: ", scalar(@GRARpositions), "\n" if ($debug);
		my %GRARpos=();
		foreach (@GRARpositions) {
			print $GRARchrom."\t"."Pos: $_\t"."Ref:${${$GRARvcf{$GRARchrom}}{$_}}[0]\t"."Var${${$GRARvcf{$GRARchrom}}{$_}}[1]\t"."Gen:${${$GRARvcf{$GRARchrom}}{$_}}[2]"."\n" if ($debug);
			my $GRARstart=$_;###add
			my $GRARend=$_;###add
			if (${${$GRARvcf{$GRARchrom}}{$_}}[0]=~/^[a-zA-Z]+$/) {###add
				$GRARend=$GRARend + length(${${$GRARvcf{$GRARchrom}}{$_}}[0])-1;###addd
			}###add
			@{$GRARpos{$_}}= $GRARsam->get_features_by_location(-seq_id => "$GRARchrom", -start => "$GRARstart", -end => "$GRARend");###mod
			print "SUB(GetRefAlleleReads)Test: ReadIDs: ".@{$GRARpos{$_}}."\n" if ($debug);
		}
		print "SUB(GetRefAlleleReads)Chrom: $GRARchrom\nPositions: @GRARpositions\n" if ($debug);
		foreach my $GRARposition (@GRARpositions) {
			print "Chr: $GRARchrom\tPosition: $GRARposition"."\t".${${$GRARvcf{$GRARchrom}}{$GRARposition}}[0]."\t".${${$GRARvcf{$GRARchrom}}{$GRARposition}}[1]."\n" if ($verbose or $debug);
			foreach my $GRARalign (@{$GRARpos{$GRARposition}}) {
#				print "Chr: $GRARchrom\tPosition: $GRARposition\t".$GRARalign->name."\n";
#				next if ($GRARalign->qual < $min_mapq);
#				print "FLAG: ".$GRARalign->flag."\n";
#				print "FLAG_str: ".$GRARalign->flag_str."\n";
				my @GRARqualscore=$GRARalign->query->qscore;
#				print "Read Quality: @GRARqualscore"."\n";
				my ($align2geno, $GRARqual)=&ReadVariantType($GRARposition, ${${$GRARvcf{$GRARchrom}}{$GRARposition}}[0], ${${$GRARvcf{$GRARchrom}}{$GRARposition}}[1], ${${$GRARvcf{$GRARchrom}}{$GRARposition}}[2], $GRARalign, 2);
				print "SUB(GetRefAlleleReads)Test: Chr: $GRARchrom Pos: $GRARposition Ref: ${${$GRARvcf{$GRARchrom}}{$GRARposition}}[0] Var: ${${$GRARvcf{$GRARchrom}}{$GRARposition}}[1] Genotype: ".$GRARalign->name.' '.$align2geno."\n\n\n" if ($verbose or $debug);###for test###
				if ($align2geno eq '0') {
#					$include{$GRARalign->name}++;
					print INCLUDE $GRARalign->name."\n";
				}
				elsif ($align2geno =~m/^\d+$/ and $align2geno>0) {
#					$exclude{$GRARalign->name}++;
					print EXCLUDE $GRARalign->name."\n";
				}
			}
		}
	}
	print "SUB(GetRefAlleleReads)Info: ENDing...\n" if ($debug);
}



### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependancy: Bio::DB::Sam
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



sub ReadVariantType {
	my ($RVTpos, $RVTref, $RVTvar, $RVTgeno, $RVTsam_align, $RVTret_code)=@_;
	$RVTret_code=1 unless (defined $RVTret_code);
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
#	print "SUB(ReadVariantType)Test: Captured String: $RVTcapture_string\n";# if ($debug); ### For test ###
	my @returnarr=();
	my $RVTreturn_type;
	EVTBLOCK1: {if ($RVTtest_var_end) {
		if ($RVT5endloss==0) {
			if (exists $RVTallele2geno{$RVTcapture_string}) {
				$RVTreturn_type=$RVTallele2geno{$RVTcapture_string};
				last EVTBLOCK1;
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
#	print "@returnarr\n";### For test ###
	
	if (scalar(@returnarr)==1) {
		$RVTreturn_type= $returnarr[0];
	}
	else {
		return 'Multiple';
	}
#	print "ReturnType: $RVTreturn_type\n";### For test ###
	}###EVTBLOCK1
	
	if ($RVTret_code ==1) {
		return $RVTreturn_type;
	}
	elsif ($RVTret_code ==2) {
#		print "ReturnCode: 2\n"; ### For test ###
		my @RVTqscore=$RVTsam_align->query->qscore;
#		print "QScore: @RVTqscore\n";### For test ###
		my @RVTqsstr=map {chr($_+33)} @RVTqscore;
#		print "ReadID: $RVTreadid\nReadSeq: $RVTquery_seq\nReadQual: ".join('',@RVTqsstr)."\n";### For test ###
		
		my $RVTreturn_score='D';
		my $EVTqs_string='';
#		print "Start-End: $EVTcapture_start_pos-$RVTcapture_end_pos\n";### For test ###
		if (($RVTcapture_end_pos-$EVTcapture_start_pos) >0) {
			my $RVTcount=0;
			my $RVTsum=0;
#			print "Capture quality\n";### For test ###
#			print "Start-end: $EVTcapture_start_pos - $RVTcapture_end_pos\n";### For test ###
			for (my $RVTi=$EVTcapture_start_pos; $RVTi<$RVTcapture_end_pos; $RVTi++) {
				if (0) {### For test ###
					$EVTqs_string.=chr($RVTqscore[$RVTi]+33);
#					print "Capture quality: $EVTqs_string\n";
				}
				$RVTsum+=$RVTqscore[$RVTi];
				$RVTcount++;
#				print "Converting $RVTi\n";### For test ###
			}
			my $RVTreturn_value=int($RVTsum/$RVTcount);
			$RVTreturn_score=chr($RVTreturn_value+33);
			if (1) {### For test ###
#				print "QualityString: $EVTqs_string\nAverage_value: $RVTreturn_value\nAverageCode: $RVTreturn_score\n";
			}
		}
		return ($RVTreturn_type, $RVTreturn_score);
	}
	else {
		##to be defined
	}
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



###Get min/max value from a arr
###&GetValues (min/max, @num_arr)
###Global: 
###Dependency:
###Nnote:
sub GetValues {
	my ($GVindex, @GVnumbers)=@_;
	@GVnumbers=sort {$a<=>$b} @GVnumbers;
	if ($GVindex=~/^min$/i) {
		return defined $GVnumbers[0] ? $GVnumbers[0]:'err';
	}
	elsif ($GVindex=~/^max$/i) {
		return defined $GVnumbers[-1] ? $GVnumbers[-1]:'err';
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



###Retrieve filebasename without extension
###& RetrvNoExt(file)
###Global:

sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my $RNE_new='';
	my $RNE_base='';
	($RNE_base=$RNE_ori)=~ s/.*\///s;
	($RNE_new=$RNE_base)=~s/^(\S+)\.\w+$/$1/;
	return $RNE_new;
}

