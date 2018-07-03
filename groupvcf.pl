#!/usr/bin/env perl
use warnings;
use strict;
use Bio::DB::Sam;
use constant USAGE =><<EOH;


Usage:

perl $0 ref bam vcf

EOH

die USAGE if (scalar(@ARGV)!=3);



###Default
our $debug=0;
our $min_share_alignments=3;
our $min_mapq=0 unless (defined $min_mapq);
our $bam_genome_tag='zw';
our $bam_chromo_tag='zc';
our $geno_delimiter=',';

###main
our $reference_file=$ARGV[0];
our $sam_file=$ARGV[1];
our $vcf_file=$ARGV[2];
my $samobj=&ReadSam($sam_file, $reference_file, 1);
my $vcfobj=&ReadVcf($vcf_file);
&GroupVcf($vcfobj,  $samobj);



###group SNPs based on read length
###&GroupVcf($ReadSam_obj, AABBDD_Vcf_obj, AABB_vcf_obj, AA_vcf_obj, DD_vcf_obj)
###Global: $debug
###Dependancy:&ReadVcf, &ReadVariantType, &ReadSam, $min_share_alignments, $min_mapq, $bam_genome_tag
sub GroupVcf {
	my ($GVaabbdd_samobj, $GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj)=@_;
	my $GVfixAlleles_index=&AssignVariationAllele($GVaabbdd_vcfobj, $GVaabb_vcfobj, $GVaa_vcfobj, $dd_vcf_obj);

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
							${${GVreadid2chrmosome{$GVind_align->name}}=$1;
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
	return \%GVfixAlleles_rnaseq, \%GVreadid_by_allele, \%GVreadid2chrmosome;
}



###Get min/max value from a arr
###&GetValues (min/max, @num_arr)
###Global: 
###Dependency:
###Nnote:
sub GetValues {
	my ($GVindex, @GVnumbers)=@_;
	return @GVnumbers if (scalar(@GVnumbers)==1);
	my @GVreturn_arr=();
	if ($GVindex=~/^min$/i) {
		@GVnumbers=sort {$a<=>$b} @GVnumbers;
	}
	elsif (($GVindex=~/^max$/i)) {
		@GVnumbers=sort {$b<=>$a} @GVnumbers;
	}
	else {
		print STDERR "SUB(GetValues)ERROR: unknown index\n";
	}
	push (@GVreturn_arr, $GVnumbers[0]);
	for (my $GVi=1;$GVi<scalar(@GVnumbers); $GVi++) {
		if ($GVnumbers[$GVi] == $GVnumbers[0]) {
			push (@GVreturn_arr, $GVnumbers[$GVi]);
		}
		else {
			last;
		}
	}
	return @GVreturn_arr;
}



###Compare two array, return unique or share
###&CompareTwoList($list1_arr_index, $list2_arr_index, 1/2/3/4)
###Global:
###Dependency:
###Note: 1=list1_unique, 2=share, 3=list2_unique, 4=ALL
sub CompareTwoList {
	my ($CTLlist1_arrindex, $CTLlist2_arrindex, $CTLmode)=@_;
	my @CTLlist1=@{$CTLlist1_arrindex};
	my @CTLlist2=@{$CTLlist2_arrindex};
	my %CTLunique1=();
	my %CTLunique2=();
	my %CTLshared=();
	foreach (@CTLlist1){
		$CTLunique1{$_}=1;
	}
	foreach (@CTLlist2) {
		if (exists $CTLunique1{$_}) {
			$CTLshared{$_}=1;
			delete $CTLunique1{$_};
		}
		else {
			$CTLunique2{$_}=1;
		}
	}
	my @CTLlist1unique=keys %CTLunique1;
	my @CTLlist2unique=keys %CTLunique2;
	my @CTLlist_share=keys %CTLshared;
	if ($CTLmode==1) {
		return \@CTLlist1unique;
	}
	elsif ($CTLmode==2) {
		return \@CTLlist_share;
	}
	elsif ($CTLmode==3) {
		return \@CTLlist2unique;
	}
	else {
		return (\@CTLlist1unique, \@CTLlist_share,\@CTLlist2unique);
	}
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
