#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150506

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Assign AABBDD alleles based on it's ancestor allele

Options:
	--help|-h
		Print this help/usage;
	--input|-i	<File>
		[Msg] Fasta ID file
	--bin|-b	<INT>
		[Opt] Load INT sequence/batch, default: 100
	--vcf1|-1	[VCF.GZ_File]
	--vcf2|-2	[VCF.GZ_File]
	--vcf3|-3	[VCF.GZ_File]
	--vcf4|-4	[VCF.GZ_File]
	--vcfqual|-q	<INT>
		[Opt] Minimum VCF quality, default: 20
	--genodelimiter	<STR>
		[Opt] default ','
	--line2start	<INT>
		[Opt] Default: 0
	--output|-o	<File>
		[Msg] Output file
	--path2tabix	</path/to/tabix>
	--path2bgzip	</path/to/bgzip>
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
my ($help, $verbose, $debug, $version);
my ($input, $binsize, $output, $vcfaabbdd, $vcfaabb, $vcfaa, $vcfdd);
my ($vcfqual, $line2start, $geno_delimiter, $path_tabix, $path_bgzip);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"bin|b:i" => \$binsize,
	"vcf1|1:s" => \$vcfaabbdd,
	"vcf2|2:s" => \$vcfaabb,
	"vcf3|3:s" => \$vcfaa,
	"vcf4|4:s" => \$vcfdd,
	"vcfqual|q:i" => \$vcfqual,
	"genodelimiter:s" => \$geno_delimiter,
	"line2start:i" => \$line2start,
	"output|o:s" => \$output,
	"path2tabix:s" => \$path_tabix,
	"path2bgzip:s" => \$path_bgzip,
#	!:s:i
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$vcfqual=20 unless (defined $vcfqual);
die "Error: invalid VCF MAPQ vcfqual|q: $vcfqual\n" unless ($vcfqual =~ m/^\d+$/);
$geno_delimiter=',' unless (defined $geno_delimiter);
$binsize=100 unless (defined $binsize);
die "Error: invalid binsize --bin/-b $binsize\n" unless ($binsize =~ m/^\d+$/);
$line2start=0 unless (defined $line2start);
die "Error: invalid line number to start --line2start $line2start\n" unless ($line2start =~ m/^\d+$/);
$path_tabix='tabix' unless (defined $path_tabix);
$path_bgzip='bgzip' unless (defined $path_bgzip);



### input and output ################################################
die "Error: invalid Fasta ID file: $input\n" unless (-s $input);
die "Error: invalid VCF1 file: $vcfaabbdd\n" unless (-s $vcfaabbdd);
die "Error: invalid VCF2 file: $vcfaabb\n" unless (-s $vcfaabb);
die "Error: invalid VCF3 file: $vcfaa\n" unless (-s $vcfaa);
die "Error: invalid VCF4 file: $vcfdd\n" unless (-s $vcfdd);
###VCF1 check
if ($vcfaabbdd=~/\.vcf$/) {
	if (&CompressVcf($vcfaabbdd, "$vcfaabbdd.gz")) {
		die "Error: Compress VCF error: $vcfaabbdd\n";
	}
	$vcfaabbdd="$vcfaabbdd.gz";
}
if ($vcfaabbdd=~/\.vcf\.gz$/) {
	unless ( -s "$vcfaabbdd.tbi") {
		if (&IndexVcfGz($vcfaabbdd)) {
			die "Error: Compress VCF error: $vcfaabbdd\n";
		}
	}
}
else {
	die "Error: unknown VCF format: $vcfaabbdd\n";
}
###VCF2 check
if ($vcfaabb=~/\.vcf$/) {
	if (&CompressVcf($vcfaabb, "$vcfaabb.gz")) {
		die "Error: Compress VCF error: $vcfaabb\n";
	}
	$vcfaabb="$vcfaabb.gz";
}
if ($vcfaabb=~/\.vcf\.gz$/) {
	unless ( -s "$vcfaabb.tbi") {
		if (&IndexVcfGz($vcfaabb)) {
			die "Error: Compress VCF error: $vcfaabb\n";
		}
	}
}
else {
	die "Error: unknown VCF format: $vcfaabb\n";
}
###VCF3 check
if ($vcfaa=~/\.vcf$/) {
	if (&CompressVcf($vcfaa, "$vcfaa.gz")) {
		die "Error: Compress VCF error: $vcfaa\n";
	}
	$vcfaa="$vcfaa.gz";
}
if ($vcfaa=~/\.vcf\.gz$/) {
	unless ( -s "$vcfaa.tbi") {
		if (&IndexVcfGz($vcfaa)) {
			die "Error: Compress VCF error: $vcfaa\n";
		}
	}
}
else {
	die "Error: unknown VCF format: $vcfaa\n";
}
###VCF4 check
if ($vcfdd=~/\.vcf$/) {
	if (&CompressVcf($vcfdd, "$vcfdd.gz")) {
		die "Error: Compress VCF error: $vcfdd\n";
	}
	$vcfdd="$vcfdd.gz";
}
if ($vcfdd=~/\.vcf\.gz$/) {
	unless ( -s "$vcfdd.tbi") {
		if (&IndexVcfGz($vcfdd)) {
			die "Error: Compress VCF error: $vcfdd\n";
		}
	}
}
else {
	die "Error: unknown VCF format: $vcfdd\n";
}



### Main ############################################################
my $total_numline=0;
open (ID, "<", $input) || die "Error: can not open ID: $input\n";
while (my $line1=<ID>) {
	$total_numline++;
}
close ID;
open (OUTPUT, '>', $output) || die "Error: can not write output: $output\n";
while ($line2start<=$total_numline) {
	my @faids=();
	my $line2end=$line2start+$binsize;
	@faids=&getfastaid($line2start, $line2end);
	my $all_fasid='';
	if (scalar(@faids) >0) {
		$all_fasid=join(' ', @faids);
	}
	else {
		print STDERR "No any sequence IDs at Line: $line2start-$line2end\n";
		print STDOUT "Problematic Line: ".$line2start."-".$line2start+$binsize."\n";
	}
###Extract VCF
	if (&ExtractVcf($vcfaabbdd, $all_fasid, "$vcfaabbdd.$line2start-$line2end.vcf.gz")) {
		die "Error: VCF extraction: $vcfaabbdd $line2start-$line2end\n";
	}
	if (&ExtractVcf($vcfaabb, $all_fasid, "$vcfaabb.$line2start-$line2end.vcf.gz")) {
		die "Error: VCF extraction: $vcfaabb $line2start-$line2end\n";
	}
	if (&ExtractVcf($vcfaa, $all_fasid, "$vcfaa.$line2start-$line2end.vcf.gz")) {
		die "Error: VCF extraction: $vcfaa $line2start-$line2end\n";
	}
	if (&ExtractVcf($vcfdd, $all_fasid, "$vcfdd.$line2start-$line2end.vcf.gz")) {
		die "Error: VCF extraction: $vcfdd $line2start-$line2end\n";
	}
###READ VCF into hash
	my ($test_rv1, $aabbdd_vcfobj)=&ReadVcf("$vcfaabbdd.$line2start-$line2end.vcf.gz");
	if ($test_rv1) {
		die "Error: Reading VCF error 1\n";
	}
	my ($test_rv2, $aabb_vcfobj)=&ReadVcf("$vcfaabb.$line2start-$line2end.vcf.gz");
	if ($test_rv2) {
		die "Error: Reading VCF error 2\n";
	}
	my ($test_rv3, $aa_vcfobj)=&ReadVcf("$vcfaa.$line2start-$line2end.vcf.gz");
	if ($test_rv3) {
		die "Error: Reading VCF error 3\n";
	}
	my ($test_rv4, $dd_vcfobj)=&ReadVcf("$vcfdd.$line2start-$line2end.vcf.gz");
	if ($test_rv4) {
		die "Error: Reading VCF error 4\n";
	}
### Merge AABBDD AABB AA DD
	foreach my $chrom (@faids) {
		unless (exists ${$aabbdd_vcfobj}{$chrom}) {
			print STDERR "Info: no variations for $chrom";
			next;
		}
		foreach my $posit (sort {$a <=> $b } keys %{${$aabbdd_vcfobj}{$chrom}}) {
			my $geno1=${${${$aabbdd_vcfobj}{$chrom}}{$posit}}[2];
			die "Error: unknown genotype $geno1 at chr:Pos $chrom:$posit\n" unless ($geno1=~/^\d+.*\d+$/);
			my @genoallele1=();
			@genoallele1=split(/\//, $geno1);
			die "Error: num_allele !=3 at genotype $geno1 at chr:Pos $chrom:$posit\n" unless (scalar(@genoallele1)==3);
			my %genonum1=();
			map {$genonum1{$_}++} @genoallele1;
			###
			my ($genoA, $genoB, $genoD)=('?', '?', '?');
			my ($geno2, $geno3, $geno4)=('?', '?', '?');
			###AABB
			if (exists ${${${$aabb_vcfobj}{$chrom}}{$posit}}[2]) {
				$geno2=${${${$aabb_vcfobj}{$chrom}}{$posit}}[2];
				my ($test2_sg, $num2allele, $unique2allele)=&SplitGeno($geno2);
				if ($test2_sg ==0 and defined $num2allele and $num2allele ==1) {
					#dosth
				}
			}
			###AA
			if (exists ${${${$aa_vcfobj}{$chrom}}{$posit}}[2]) {
				$geno3=${${${$aa_vcfobj}{$chrom}}{$posit}}[2];
				my ($test3_sg, $num3allele, $unique3allele)=&SplitGeno($geno3);
				if ($test3_sg ==0 and defined $num3allele and $num3allele ==1) {
					$genoA=$unique3allele if (exists $genonum1{$unique3allele});
				}
			}
			if (exists ${${${$dd_vcfobj}{$chrom}}{$posit}}[2]) {
				$geno4=${${${$dd_vcfobj}{$chrom}}{$posit}}[2];
				my ($test4_sg, $num4allele, $unique4allele)=&SplitGeno($geno4);
				if ($test4_sg ==0 and defined $num4allele and $num4allele ==1) {
					$genoD=$unique4allele if (exists $genonum1{$unique4allele});
				}
			}
			print OUTPUT "$chrom\t$posit\t$genoA\t$genoB\t$genoD\t$geno1\t$geno2\t$geno3\t$geno4\n";
		}
	}
	print STDOUT "Success Line: ".$line2start."-".$line2start+$binsize."\n";
	$line2start=$line2end;
	unlink "$vcfaabbdd.$line2start-$line2end.vcf.gz";
	unlink "$vcfaabb.$line2start-$line2end.vcf.gz";
	unlink "$vcfaa.$line2start-$line2end.vcf.gz";
	unlink "$vcfdd.$line2start-$line2end.vcf.gz";
}
close OUTPUT;



#####################################################################
###                         sub functions                         ###
#####################################################################

###Read vcf into hash, return index
###&ReadVcf(vcf_file)
###Global: $debug, $vcfqual, $geno_delimiter
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
	my $RVtest_cmd=0;
	my %RVret_vcf=();
##COMMENT: close filehandle before reopen
	if (defined fileno(VCF)) {
		unless (close VCF) {
			print STDERR "SUB(ReadVcf)Error: can not close VCF filehandle for $RVvcf\n";
			$RVtest_cmd=1;
			return ($RVtest_cmd, \%RVret_vcf);
		}
	}
##COMMENT: guess format
	if ($RVvcf=~/\.gz$|\.gzip/i) {
		unless (open (VCF, "zcat $RVvcf|")) {
			print STDERR "SUB(ReadVcf)Error: can not open vcf.gz file\n";
			$RVtest_cmd=1;
			return ($RVtest_cmd, \%RVret_vcf);
		}
	}
	elsif ($RVvcf=~/\.vcf$/i) {
		unless (open (VCF, "cat $RVvcf|")) {
			print STDERR "SUB(ReadVcf)Error: can not open vcf file\n";
			$RVtest_cmd=1;
			return ($RVtest_cmd, \%RVret_vcf);
		}
	}
	else {
		print STDERR "SUB(ReadVcf)Error: can not guess vcf file format\n";
		$RVtest_cmd=1;
		return ($RVtest_cmd, \%RVret_vcf);
	}
##COMMENT: reading VCF
	my $RVnum_line=0;
	while (my $RVline=<VCF>) {
		$RVnum_line++;
		next if ($RVline =~/^#/);
		chomp $RVline;
		my @RVarr=();
		@RVarr=split(/\t/, $RVline);
##COMMENT: comment next line if column number <10
		if (scalar(@RVarr)<10) {
			print STDERR "SUB(ReadVcf)Error: Wrong variation at $RVnum_line in $RVvcf\n";
			next;
		}
##COMMENT: 
		my ($RVchr, $RVpos, $RVref, $RVvar)=('', '', '', '');
		($RVchr, $RVpos, $RVref, $RVvar)=($RVarr[0], $RVarr[1], $RVarr[3], $RVarr[4]);
		if (exists ${$RVret_vcf{$RVchr}}{$RVpos} ) {
			print STDERR "SUB(ReadVcf)Error: repeated $RVchr - $RVpos at $RVnum_line in $RVvcf\n";
			$RVtest_cmd=1;
			return ($RVtest_cmd, \%RVret_vcf);
		}
		@{${$RVret_vcf{$RVchr}}{$RVpos}}=($RVref, $RVvar);
##COMMENT: Parse vcf
		for (my $i=9; $i<scalar(@RVarr); $i++) {
			#GT:DP:RO:QR:AO:QA=0/0/1:90:47:1816:42:1060
			if ($RVarr[$i] eq '.') {
				push (@{${$RVret_vcf{$RVchr}}{$RVpos}}, '.');
				next;
			}
			my @RVarr2=();
			@RVarr2=split(/:/, $RVarr[$i]);
			#geno: $RVarr2[0]
			#readdepth: $RVarr2[1]
			#ref count: $RVarr2[2]
			#var count: $RVarr2[4]
##COMMENT: following two lines if not ignore those sites with only one var allele
#			unless (&IsZeroIn($RVarr2[2])) {
#				delete ${$RVret_vcf{$RVchr}}{$RVpos};
#				last;
#			}
#			unless (&IsZeroIn($RVarr2[4])){
#				delete ${$RVret_vcf{$RVchr}}{$RVpos};
#				last;
#			}
			if (defined $RVarr[5] and $RVarr[5]<$vcfqual) {##COMMENT: check mapping quality; would take the allele with max count
				my @arr3=split(/$geno_delimiter/, $RVvar);
				my @arr4=split(/$geno_delimiter/, $RVarr2[4]);
				unless (scalar(@arr3) == scalar(@arr4)) {
					print STDERR "SUB(ReadVcf)Error: VCF allele count $RVchr - $RVpos at $RVnum_line in $RVvcf\n";
					$RVtest_cmd=1;
					return ($RVtest_cmd, \%RVret_vcf);
				}
				my %RVallele_count=();
				for (my $RVj=0; $RVj<scalar(@arr4); $RVj++) {
					$RVallele_count{$RVj+1}=$arr4[$RVj];
				}
				my $RVbest_allele=0;
				my $RVmax_count=$RVarr2[2];
				foreach (keys %RVallele_count) {
					if ($RVallele_count{$_} > $RVmax_count) {
						$RVbest_allele=$_;
					}
				}
				my @RVarr5=split(/\//, $RVarr2[0]);
				my @RVarr6=();
				foreach (@RVarr5) {
					push (@RVarr6, $RVbest_allele);
				}
				$RVarr2[0]=join('/', @RVarr6);
			}
			push (@{${$RVret_vcf{$RVchr}}{$RVpos}}, $RVarr2[0]);
		}
	}
	close VCF;
	
	#	print "vcf input\n";###test###
	if ($debug) {
		foreach my $RVchrom (keys %RVret_vcf) {
			my @RVarr3=keys %{$RVret_vcf{$RVchrom}};
			@RVarr3=sort {$a<=>$b} @RVarr3;
			print "SUB(ReadVcf)Reference: $RVchrom\n";
			print "SUB(ReadVcf)Number of variations on $RVchrom: ", scalar(@RVarr3), "\n";
			foreach (@RVarr3) {
				print $RVchrom, "\t", "$_\t${${$RVret_vcf{$RVchrom}}{$_}}[0]\t${${$RVret_vcf{$RVchrom}}{$_}}[1]\t${${$RVret_vcf{$RVchrom}}{$_}}[2]\n";
			}
		}
	}
	return ($RVtest_cmd, \%RVret_vcf);
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
###$vcf{chr}->{pos}=(ref, var, gen/gen2, ...);
###Be noted: var could be comma-delimited string, in case of two variants
}



### Read sequence ids from line_start to line_end
### &getfastaid(linenum_start, linenum_end)
### Global: $input
### Dependent: 
### Note:
sub getfastaid {
	my ($GFIstart, $GFIend)=@_;
	my @GFIfastaid=();
	open (IDLIST, "<", $input) || die "SUB(getfastaid)Error: can not open IDlist: $input\n";
	my $GFInumline=0;
	while (my $GFIline=<IDLIST>) {
		$GFInumline++;
		chomp $GFIline;
		if ( $GFInumline>$GFIstart and $GFInumline<=$GFIend ) {
			if ($GFIline=~m/^#/) {
				next;
			}
			else {
				(my $GFIind_id=$GFIline)=~s/\s+.*$//;
				push (@GFIfastaid, $GFIind_id);
			}
		}
		else {
			next;
		}
	}
	close IDLIST || die "SUB(getfastaid)Error: can not close IDLIST\n";
	return @GFIfastaid;
}



###	compress vcf
###	CompressVcf(input.vcf, output.vcf.gz)
### Global:
### Dependent: $path_bgzip
### Note:
sub CompressVcf {
	my ($CVvcffile, $CVvcfoutgz)=@_;
	die "SUB(CompressVcf)Error: file not found: $CVvcffile\n" unless ( -s $CVvcffile );
	if ( -e $CVvcfoutgz ) {
		print STDERR "SUB(CompressVcf)Error: output exists: $CVvcfoutgz; deleting...\n";
		unlink $CVvcfoutgz;
	}
	if (&exec_cmd_return("cat $CVvcffile | $path_bgzip > $CVvcfoutgz")) {
		print STDERR "SUB(CompressVcf)Error: bgzip running: $CVvcfoutgz\n";
		return 1;
	}
	elsif ( -s $CVvcfoutgz ) {
		if (&IndexVcfGz($CVvcfoutgz)) {
			die "SUB(CompressVcf)Error: tabix index error: $CVvcfoutgz\n";
		}
		else {
			return 0;
		}
	}
	else {
		print STDERR "SUB(CompressVcf)Error: bgzip output not found: $CVvcfoutgz\n";
		return 1;
	}
}



### Split VCF geno 0/1/1
### SplitGeno (0/1/1)
### Global:
### Dependency:
### Note:
### Return: test_cmd, num_allele, unique_allele
sub SplitGeno {
	my $SGgeno=shift;
	return 1 unless ($SGgeno=~/^\d+.*\d+$/);
	my @SGarr=split(/\//, $SGgeno);
	return 1 unless (scalar(@SGarr) == 2);
	my %SGhash=();
	map {$SGhash{$_}++} @SGarr;
	my $SGnum_alleles=scalar(keys %SGhash);
	if ($SGnum_alleles==1) {
		return (0, 1, $SGarr[0]);
	}
	elsif ($SGnum_alleles==2) {
		return (0, 2, $SGarr[0]);
	}
	else {
		die "SUB(SplitGeno)Error: unknown genotype: $SGgeno\n";
	}
}



###	Index vcf
###	CompressVcf(input.vcf.gz)
### Global:
### Dependent: $path_tabix
### Note:
sub IndexVcfGz {
	my $IVGvcfgz_file=shift;
	die "SUB(IndexVcfGz)Error: file not found: $IVGvcfgz_file\n" unless ( -s $IVGvcfgz_file);
	unlink "$IVGvcfgz_file.tbi" if ( -e "$IVGvcfgz_file.tbi");
	if (&exec_cmd_return("$path_tabix -p vcf $IVGvcfgz_file")) {
		print STDERR "SUB(IndexVcfGz)Error: tabix index running: $IVGvcfgz_file\n";
		return 1;
	}
	elsif ( -s "$IVGvcfgz_file.tbi") {
		return 0;
	}
	else {
		print STDERR "SUB(IndexVcfGz)Error: tabix output not found: $IVGvcfgz_file.tbi\n";
		return 1;
	}
}




### Read sequence ids from line_start to line_end
### &ExtractVcf(source.vcf.gz, ids[space delimited], output.vcf.gz)
### Global:
### Dependent: $path_tabix, $path_bgzip
### Note:
sub ExtractVcf {
	my ($EVvcffile, $EVids, $EVvcfout)=@_;
	die "SUB(ExtractVcf)Error: file not found: $EVvcffile\n" unless ( -s $EVvcffile );
	die "SUB(ExtractVcf)Error: empty FastaIDs to extract\n" if ($EVids=~m/^\s*$/);
	if ( -e $EVvcfout ) {
		print STDERR "SUB(ExtractVcf)Error: output exists: $EVvcfout; deleting...\n";
		unlink $EVvcfout;
	}
	if (&exec_cmd_return("$path_tabix $EVvcffile $EVids | $path_bgzip > $EVvcfout")) {
		print STDERR "SUB(ExtractVcf)Error: tabix extract running: $EVvcfout\n";
		return 1;
	}
	elsif ( -s $EVvcfout ) {
		return 0;
	}
	else {
		print STDERR "SUB(ExtractVcf)Error: output not found: $EVvcfout\n";
		return 1;
	}
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
		print STDERR "Finished command: $cmd\tat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
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
		print STDERR "Finished command: $cmd\tat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}
