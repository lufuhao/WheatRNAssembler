#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::VcfKit qw/LoadVcf/;
use constant USAGE =><<EOH;

usage: $0 merge.caller.vcf good.vcf bad.vcf allele.sum

EOH
my $debug=0;
my $vcfinput=$ARGV[0];
my $vcfoutput=$ARGV[1];
my $vcfproblem=$ARGV[2];
my $vcfsum=$ARGV[3];
die "Error: invalid vcfinput\n" unless (defined $vcfinput and -s $vcfinput);

my ($test_cmd, $vcfobj)=&LoadVcf($vcfinput, 0, 0);
unless ($test_cmd) {
	die "Error: loading VCF error\n";
}

if ($vcfinput=~/\.vcf$/i) {
	print "Info: VCF input in VCF format\n" if ($debug);
	open (VCFIN, "< $vcfinput") || die "Error: can not open vcfinput: $vcfinput\n";
}
elsif ($vcfinput=~/\.vcf\.gz$/i) {
	print "Info: VCF input in VCF.GZ format\n" if ($debug);
	open (VCFIN, "zcat $vcfinput |") || die "Error: can not open vcfinput: $vcfinput\n";
}
open (VCFGOOD, " | bgzip > $vcfoutput") || die "Error: can not write good.vcf: $vcfoutput\n";
open (VCFBAD, " > $vcfproblem") || die "Error: can not write bad.vcf: $vcfproblem\n";
open (VCFSUM, " > $vcfsum") || die "Error: can not write vcfoutput: $vcfsum\n";
my $linenum=0;
my $assignedAA=0;
my $assignedDD=0;
while (my $line=<VCFIN>) {
	$linenum++;
	if ($line=~/^#/) {
		print VCFGOOD $line;
		next;
	}
	else {
		my $aabbdd_origeno='undef';
		my $aabb_origeno='undef';
		my $aa_origeno='undef';
		my $dd_origeno='undef';
		my $aabb_revgeno='undef';
		my $aa_revgeno='undef';
		my $dd_revgeno='undef';
		my $aabb_numallele='undef';
		my $aa_numallele='undef';
		my $dd_numallele='undef';
		chomp $line;
		my @arr1=split(/\t/, $line);
		if (scalar(@arr1)<13) {
			print STDERR "Warnings: colnum<13 at line $linenum of VCF $vcfinput\n";
			print VCFBAD $line."\n";
			next;
		}
		if (exists ${$vcfobj}{$arr1[0]} and exists ${$vcfobj}{$arr1[0]}{$arr1[1]} and exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABBDD'} and exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABBDD'}{'genotypes'}) {
			unless (${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABBDD'}{'genotypes'} =~/\//) {
				print STDERR "Warnings: unknown geno at line $linenum of VCF $vcfinput\n";
				print VCFBAD $line."\n";
				next;
			}
			$aabbdd_origeno=${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABBDD'}{'GT'};
			print VCFGOOD "$arr1[0]\t$arr1[1]\t$arr1[2]\t$arr1[3]\t$arr1[4]\t$arr1[5]\t$arr1[6]\t$arr1[7]\t$arr1[8]\t$aabbdd_origeno\n";
			my @arr2=split(/\t/, $aabbdd_origeno);
			my %genohash1=();
			map {$genohash1{$_}++} @arr2;
			
			if (exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABB'} and exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABB'}{'genotypes'}) {
				$aabb_origeno=${$vcfobj}{$arr1[0]}{$arr1[1]}{'AABB'}{'GT'};
				my ($test_code1, $aabb_numallele, $geno1, $uniqallele1)=&SplitGeno($aabb_origeno);
				if ($test_code1==2) {
					$aabb_revgeno='.';
				}
				elsif ($test_code1==0) {
					$aabb_revgeno=$geno1;
				}
			}
			else {
				print STDERR "Warnings: no AABB geno at line $linenum of VCF $vcfinput\n";
			}
			if (exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'AA'} and exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'AA'}{'genotypes'}) {
				$aa_origeno=${$vcfobj}{$arr1[0]}{$arr1[1]}{'AA'}{'GT'};
				my ($test_code1, $aa_numallele, $geno1, $uniqallele1)=&SplitGeno($aa_origeno);
				if ($test_code1==2) {
					$aa_revgeno='?';
				}
				elsif ($test_code1==0) {
					if ($aa_numallele==1 and defined $uniqallele1 and exists $genohash1{$uniqallele1}) {
						$assignedAA++;
					}
					$aa_revgeno=$geno1;
				}
			}
			else {
				print STDERR "Warnings: no AA geno at line $linenum of VCF $vcfinput\n";
			}
			if (exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'DD'} and exists ${$vcfobj}{$arr1[0]}{$arr1[1]}{'DD'}{'genotypes'}) {
				$dd_origeno=${$vcfobj}{$arr1[0]}{$arr1[1]}{'DD'}{'GT'};
				my ($test_code1, $dd_numallele, $geno1, $uniqallele1)=&SplitGeno($dd_origeno);
				if ($test_code1==2) {
					$dd_revgeno='?';
				}
				elsif ($test_code1==0) {
					if ($dd_numallele==1 and defined $uniqallele1 and exists $genohash1{$uniqallele1}) {
						$assignedDD++;
					}
					$dd_revgeno=$geno1;
				}
			}
			else {
				print STDERR "Warnings: no DD geno at line $linenum of VCF $vcfinput\n";
			}
		}
		print VCFSUM "$arr1[0]\t$arr1[1]\t$arr1[3]\t$arr1[4]\t$aabbdd_origeno\t$aabb_origeno\t$aa_origeno\t$dd_origeno\t$aabb_numallele\t$aa_numallele\t$dd_numallele\t$aabb_revgeno\t$aa_revgeno\t$dd_revgeno\n";
	}
}




sub SplitGeno {
	my $SGgeno=shift;
	my $SGallelenum=0;
	my $SGfix_geno='';
	unless ($SGgeno =~ /\//) {
		return (2, $SGallelenum, $SGfix_geno);
	}
	my @SGarr=split(/\//, $SGgeno);

	my %SGgenohash=();
	
	map {$SGgenohash{$_}++} @SGarr;
	
	my @SGarr2=keys %SGgenohash;
	$SGallelenum=scalar(@SGarr2);
	if ($SGallelenum==1) {
		$SGfix_geno=$SGarr2[0].'/'.$SGarr2[0];
		return (0, $SGallelenum, $SGfix_geno, $SGarr2[0]);
	}
	elsif ($SGallelenum>1) {
		$SGfix_geno=join('/', @SGarr2);
		return (0, $SGallelenum, $SGfix_geno);
	}
	else {
		return (1, $SGallelenum, $SGfix_geno);
	}
}
