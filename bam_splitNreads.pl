#!/usr/bin/env perl
use strict;
use warnings;
#use Getopt::Long;
#use Cwd;
#use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS:

	samtools view -h xx.bam | perl $0 | samtools calmd - xx.reference.fa | samtools sort -o - xx.sort.bam
	
Version: LUFUHAO20150331

Requirements:
	Programs: Perl

Descriptions:
	Split CIGAR containing N into separate parts, fields affected:
		FLAG: complementary aplignment will be with 0x800 flag
		POS:  base on CIGAR
		CIGAR:recalculate, introduce H
		TLEN: recalculate
		SEQ:  substring based on CIGAR
		QUAL: substring based on CIGAR
	
	Note: 
		Need to recalculate MD using samtools calmd
		output in sam format

Example:
	samtools view -h xx.bam | perl $0 | samtools calmd - xx.reference.fa | samtools sort -o - xx.sort.bam

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
my $help=0;
my $version=0;
my $debug=0;
#our ($help, $verbose, $debug, $version);
#our ($input, $output);

#GetOptions(
#	"help|h!" => \$help,
#	"input|i:s" => \$input,
#	"output|o:s" => \$output,
#	!:s:i
#	"debug!" => \$debug,
#	"verbose!" => \$verbose,
#	"version|v!" => \$version) or die USAGE;
#($help or $version) and die USAGE;



### Defaults ########################################################
my $complementaryflag=2048;
###If want to include the N len in splitted CIGAR, set $includingNlen=1;
my $includingNlen=0; 



### input and output ################################################
while (my $line=<>) {
	chomp $line;
	if ($line=~/^\@/) {
		print $line."\n";
		next;
	}
	if ($line=~/RG:A:/)	{
		$line=~s/RG:A:/RG:Z:/g;
	}
	my @arr1=split(/\t/, $line);
###Verify CIGAR length
#$arr1[0]	readid
#$arr1[1]	FLAG		***
#$arr1[2]	reference
#$arr1[3]	pos			***
#$arr1[4]	mapq
#$arr1[5]	CIGAR		***
#$arr1[6]	mate name
#$arr1[7]	mate pos
#$arr1[8]	insert		***
#$arr1[9]	seq			***
#$arr1[10]	qual		***
#$arr1[11]	Tag

### Arrays to be changed
	my @flags=();
	my @positions=();
	my @cigars=();
	my @tlens=();
	my @sequence=();
	my @quals=();

### split cigar
	my $seq_length=length($arr1[9]);
	my ($cigar_len, $readpartlen_arrindex, $refpartlen_arrindex, $npartlen_arrindex)=&SplitCigar($arr1[5]);
	if ($cigar_len != $seq_length) {###verify cigar len == sequence len
		print STDERR $line."\n";
		next;
	}

###Parse different partial alignments
	if ($arr1[5]=~/N/) {
		my $total_parts=$arr1[5]=~s/N/N/g;
		$total_parts+=1;
		my @cigar_parts=split(/\d+N/, $arr1[5]);
		my @reflen_arr=@{$refpartlen_arrindex}; undef $refpartlen_arrindex;
#		print STDERR "ref_len\n"; map {print STDERR $_."\n"} @reflen_arr;### For test ###
		my @readpartlen=@{$readpartlen_arrindex}; undef $readpartlen_arrindex;
		my @npartlen=@{$npartlen_arrindex}; undef $npartlen_arrindex;
		push (@npartlen, 0);
		unless ($total_parts == scalar(@cigar_parts) or $total_parts == scalar(@reflen_arr) or $total_parts == scalar(@readpartlen) or $total_parts == scalar(@npartlen)) {
			die "Split N error at $line\n";
		}
		my $new_pos=$arr1[3];
		my $part_reflen=0;
		my $cigar_prefix=0;
		my $part_nlen=0;
		my $new_tlen= ($arr1[7]>=$arr1[3]) ? $arr1[8] : $arr1[7]-$arr1[3];
		my $newseq_start=0;
		for (my $i=0; $i<$total_parts; $i++) {
			my $cigar_suffix=0;
			if ($i==0) {
				push (@flags, $arr1[1]);
				push (@positions, $new_pos);
			}
			else {
				push (@flags, $arr1[1]+$complementaryflag);
				$new_pos+=$part_reflen;
				$new_pos+=$part_nlen;
				push (@positions, $new_pos);
			}
			my $newcigar=shift @cigar_parts;
			my $part_readlen=shift @readpartlen;
			$part_reflen=shift @reflen_arr;
			if ($i==$total_parts-1) {
				$cigar_suffix= 0;
			}
			else {
				if (scalar(@readpartlen)>0) {
					foreach (@readpartlen) {
						die "Invalid number: $_ at line: $line\n" unless (/^\d+$/);
						$cigar_suffix+=$_;
					}
				}
				if ($includingNlen) {
					if (scalar(@npartlen)>0) {
						foreach (@npartlen) {
							die "Invalid number2: $_ at line: $line\n" unless (/^\d+$/);
							$cigar_suffix+=$_;
						}
					}
				}
			}
			$cigar_suffix= ($cigar_suffix>0) ? $cigar_suffix : 0;
			$newcigar=($cigar_prefix>0) ? $cigar_prefix.'H'.$newcigar : $newcigar;
			$newcigar=($cigar_suffix>0) ? $newcigar.$cigar_suffix.'H' : $newcigar;
			push (@cigars, $newcigar);
			if ($arr1[7]<$arr1[3]) {
				$new_tlen-=$part_reflen+$part_nlen;
				push (@tlens, $new_tlen);
			}
			$part_nlen=($i==$total_parts-1) ? 0 : shift @npartlen;
			
			$cigar_prefix+=$part_readlen;
			$cigar_prefix+=$part_nlen if ($includingNlen);
			if ($arr1[7]>=$arr1[3]) {
				push (@tlens, $new_tlen);
				$new_tlen-=$part_reflen+$part_nlen;
			}
			my $newseq=substr($arr1[9], $newseq_start, $part_readlen);
			my $newqual=substr($arr1[10], $newseq_start, $part_readlen);
			if (defined $newseq and $newseq ne '' and defined $newqual and $newqual ne '') {
				push (@sequence, $newseq);
				push (@quals, $newqual);
				$newseq_start+=$part_readlen;
			}
			else {
				die "Empty sequence or qual at line: $line\n";
			}

		}
		for (my $j=0; $j<$total_parts; $j++) {
			my @arr2=@arr1;
###Verify CIGAR length
#$arr2[0]	readid
#$arr2[1]	FLAG		***
#$arr2[2]	reference
#$arr2[3]	pos			***
#$arr2[4]	mapq
#$arr2[5]	CIGAR		***
#$arr2[6]	mate name
#$arr2[7]	mate pos
#$arr2[8]	insert		***
#$arr2[9]	seq			***
#$arr2[10]	qual		***
#$arr2[11]	Tag
			$arr2[1]=$flags[$j];
			$arr2[3]=$positions[$j];
			$arr2[5]=$cigars[$j];
			$arr2[8]=$tlens[$j];
			$arr2[9]=$sequence[$j];
			$arr2[10]=$quals[$j];
			print join("\t",@arr2), "\n";
		}
	}
	else {
		print $line."\n";
	}
	if ($debug) {### For test ###
		print "\n\n\n###SUMMARY #####\n";
		print "FLAG: @flags\n";
		print "Posi: @positions\n";
		print "CIGA: @cigars\n";
		print "Tlen: @tlens\n";
		print "Sequ: @sequence\n";
		print "Qual: @quals\n";
		print "###SUMMARY #####\n\n\n";
	}


}
exit 0;




### Main ############################################################




#####################################################################
###                         sub functions                         ###
#####################################################################
###split cigar into double array ((2, M), (5, D), ....)
###&SplitCigar(SamCigar)
###Global: None
###Dependancy: none
sub SplitCigar($) {
	my $SCcigar_string = shift;
	my @SCreturn_cigar_arr=();
	my $SCtaglength=0;
##Format: store the length of reads mapped to each part 30M*N40M*N20M => (30, 40, 20)
	my @SCsplitCigarLength=();
	my $SCindquery_length=0;
##Format: store the length of corresponding ref separated by N cigar: 2D1I30M*N40M*N20M => (32,40,20)
	my @SCrefCigarLength=();
##Format: store the N lenth
	my @SCnlen=();
	my $SCindref_length=0;
	while ($SCcigar_string =~ /(\d+)([MIDNSHP=X])/g) {
		my @SCoperation=($1, $2);
		my $SCcigar_num=$1;
		my $SCcigar_chr=$2;
		push @SCreturn_cigar_arr, \@SCoperation;
		if ($SCcigar_chr =~ m/^(M)|(I)|(S)|(=)|(X)$/) {
#		Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			$SCtaglength+=$SCcigar_num;
			$SCindquery_length+=$SCcigar_num;
		}
		if ($SCcigar_chr =~ m/^(M)|(D)|(=)|(X)$/) {
			$SCindref_length+=$SCcigar_num;
		}
		if ($SCcigar_chr =~ m/^N$/) {
			push (@SCsplitCigarLength, $SCindquery_length);
			push (@SCrefCigarLength, $SCindref_length);
			push (@SCnlen,  $SCcigar_num);
			$SCindref_length=0;
			$SCindquery_length=0;
		}
	}
	push (@SCsplitCigarLength, $SCindquery_length);###last N part
	push (@SCrefCigarLength, $SCindref_length);
	return ($SCtaglength, \@SCsplitCigarLength, \@SCrefCigarLength, \@SCnlen);
}
=example
3NG5HQ1:156:C2C0AACXX:7:1204:7945:91165	99	5594_1al	147	40	26S42M12D2M316N21M	=	608	195	GAGAAAACAAGAGTTCTCCAATCGACACAATACATCGAGGGAGCAGTAAATAAACTCCCTTTCTTATGATATGTATATGATCTTCCTTTCC	HHGJJIIIJJJGGFHIGIIIIJGGIFIIJIIG@FGEHGIIIEGGHHEEHFFFFFFFCDDEDD@CDACDDEECD>BDECDDDCED@ACDDDC	RG:Z:ParLeaf1	MD:Z:9A32^ACACTAGGAATA18C4	NH:i:1	HI:i:1NM:i:14	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:Z:M
3NG5HQ1:156:C2C0AACXX:7:1204:7945:91165	99	5594_1al	147	40	26S42M12D2M337H	=	608	195	GAGAAAACAAGAGTTCTCCAATCGACACAATACATCGAGGGAGCAGTAAATAAACTCCCTTTCTTATGAT	HHGJJIIIJJJGGFHIGIIIIJGGIFIIJIIG@FGEHGIIIEGGHHEEHFFFFFFFCDDEDD@CDACDDE	RG:Z:ParLeaf1	MD:Z:9A32^ACACTAGGAATA18C4	NH:i:1	HI:i:1	NM:i:14	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M
3NG5HQ1:156:C2C0AACXX:7:1204:7945:91165	2147	5594_1al	519	40	372H21M	=	608	-177	ATGTATATGATCTTCCTTTCC	ECD>BDECDDDCED@ACDDDC	RG:Z:ParLeaf1	MD:Z:9A32^ACACTAGGAATA18C4	NH:i:1	HI:i:1	NM:i:14	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M



3NG5HQ1:156:C2C0AACXX:7:1101:16090:17761	163	192028_1al	1014	40	66M80N24M	=	1054	131	TCCCTGGTATCGAAGAGGTCAACATCTTTAAGGATGATGTGGTTATTCAGTTTCTCAATCCTAAAGTGCAAGCTTCGATTGGTGCTAATA	HHHJJJJHIJJIJJJJJIIIJJJIJJJJJJJJIIIJJJIFFHIJHIIIJIJJJIIJIHIIJIJHCEEHED@CDF@CCBEDBDCDCC<>AC	RG:Z:ParLeaf1	MD:Z:90	NH:i:1	HI:i:1	NM:i:0	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+
3NG5HQ1:156:C2C0AACXX:7:1101:16090:17761	83	192028_1al	1054	40	26M80N57M225N8M	=	1014	-131	GGTTATTCAGTTTCTCAATCCTAAAGTGCAAGCTTCGATTGGTGCTAATACATGGGTGGTCAGTGGAACTCCACAGACAAAGAAACTGCAA	:DDEDDDDDDDDDEEEEEEFFFFFFEHHHHJJJJJJJJJHIJJJJJJIJIJJJJJJJJJIJJJJJJJIJJJJHJJJJJJJJJJJIJJJHHH	RG:Z:ParLeaf1	MD:Z:91	NH:i:1	HI:i:1	NM:i:0	SM:i:40XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:Z:M
3NG5HQ1:156:C2C0AACXX:7:1101:16090:17761	83	192028_1al	1054	40	26M370H=	1014	-66	GGTTATTCAGTTTCTCAATCCTAAAG	:DDEDDDDDDDDDEEEEEEFFFFFFE	RG:Z:ParLeaf1	MD:Z:91	NH:i:1	HI:i:1	NM:i:0	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M
3NG5HQ1:156:C2C0AACXX:7:1101:16090:17761	2131	192028_1al	1160	40	106H57M233H	=	1014	-203	TGCAAGCTTCGATTGGTGCTAATACATGGGTGGTCAGTGGAACTCCACAGACAAAGAAACTGCAA HHHHJJJJJJJJJHIJJJJJJIJIJJJJJJJJJIJJJJJJJIJJJJHJJJJJJJJJJJIJJJHHH	RG:Z:ParLeaf1	MD:Z:91	NH:i:1	HI:i:1	NM:i:0	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M
3NG5HQ1:156:C2C0AACXX:7:1101:16090:17761	2131	192028_1al	1442	40	388H8M	=	1014	-436	AACTGCAA	JIJJJHHH	RG:Z:ParLeaf1	MD:Z:91	NH:i:1	HI:i:1	NM:i:0SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M


3NG5HQ1:156:C2C0AACXX:7:1103:10759:98553	83	192028_1al	1054	40	26M80N57M225N8M	=	659	-224	GGTTATTCAGTTTCTGAATCCTAAAGTGCAAGCTTCGATTGGTGCTAATACATGGGTGGTCAGTGGAACTCCACAGACAAAGAAACTGCAA	CDDEDDDDDDDDDEEEEEEFFFFFFEHHHHJJJJJIJJJIIJJJIJJIIIJJIIJJJJJJJJJJIJIJJJJJIJJJJJJJJJJJJJJJHHH	RG:Z:ParLeaf1	MD:Z:15C75	NH:i:1	HI:i:1	NM:i:1SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:Z:M
3NG5HQ1:156:C2C0AACXX:7:1103:10759:98553	83	192028_1al	1054	40	26M370H=	659	-421	GGTTATTCAGTTTCTGAATCCTAAAG	CDDEDDDDDDDDDEEEEEEFFFFFFE	RG:Z:ParLeaf1	MD:Z:15C75	NH:i:1	HI:i:1	NM:i:1	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M
3NG5HQ1:156:C2C0AACXX:7:1103:10759:98553	2131	192028_1al	1160	40	106H57M233H	=	659	-558	TGCAAGCTTCGATTGGTGCTAATACATGGGTGGTCAGTGGAACTCCACAGACAAAGAAACTGCAA	HHHHJJJJJIJJJIIJJJIJJIIIJJIIJJJJJJJJJJIJIJJJJJIJJJJJJJJJJJJJJJHHH	RG:Z:ParLeaf1	MD:Z:15C75	NH:i:1	HI:i:1	NM:i:1	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M
3NG5HQ1:156:C2C0AACXX:7:1103:10759:98553	2131	192028_1al	1442	40	388H8M=659	-791	AACTGCAA	JJJJJHHH	RG:Z:ParLeaf1	MD:Z:15C75	NH:i:1	HI:i:1NM:i:1	SM:i:40	XQ:i:40	X2:i:0	XO:Z:CU	XS:A:+	PG:A:M
=cut
