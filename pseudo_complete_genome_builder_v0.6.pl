#!/usr/bin/perl
use strict;
use List::Util qw/max min/;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
my $sca = "";
my $ref = "";
my $max = 0;
my $min = 0;
my $size = 3000;
my $filter = 3000;
my $chrnum = 1;
my $assembler = "spades";
my $trimlen = 150;
GetOptions('s=s' => \$sca, 'f=s' => \$ref, 'max=f' => \$max, 'min=f' => \$min, 'a=s' => \$assembler, 'l=f' => \$size, 'p=f' => \$filter, 'n=f' => \$chrnum, 't=f' => \$trimlen);
system("rm variant_call* pseudo_chromosome_assembly_chr*.fasta 1line.fa tempseq.fa* seq seq2 sequence sequence2 pseudo_chromosome_assembly.fasta temp_sort* tempfile mapped.paf filtered.paf");
my $tmp2 = $ref . "_";
system("rm -rf temp  report $tmp2* && mkdir temp && mkdir report");

#現在はreference配列は1つでしか機能しない。
#referenceのchr名はchr1かchromosome1にしか対応しない。
#現在、invはpositionのみcallされる。
#chr1とchr2まで対応している。chr2の方が少しcoverageが低いなら、coverageチェックで除外されてしまう。気をつける必要がある。

&title;

&change;

&map;

&filterpaf;

&get_position;

&sort1;
&sort2;

&exort_seq_chr1;
&exort_seq_chr2 if($chrnum == 2);

&concatenate_chr1;
&concatenate_chr2 if($chrnum == 2);

&seqret;

&variant_call;

&variant_filter1;
&variant_filter2 if($chrnum == 2);
system("rm $tmp2*");
exit;






#################################################################################################################################################
### SUBROUTINES ###
sub title {
	####################################################################################################################################
	#
	#						pseudo_complete_genome_builder version xxx
	#
	#						Kazuma Uesaka
	#						University of Nagoya
	#						1 September 2019
	#		
	#						A Perl scripts to recover pseudo-chromome assembly
	#
	#	
	#		pseudo_complete_genome_builder: Sensitive reference based pseudo-chromome assembly tool.
	#		Kazuma Uesaka, Tatsuo Omata, and Kunio Ihara
	#
	#						Dependency:
	#						 minimap2
	#						 paftools
	#						 samtools
	#						 seqkit
	#						 EMBOSS seqret
	#
	#						Input:
	#							reference file and de novo assembled contigs
	#
	#						Outnput:
	#							
	#						
	#						Usage:
	#						perl pseudo_complete_genome_builder.pl -f reference.fasta -s scaffolds.fasta -max 50 -min 40 -a spades -l 3000 -n 2 -p 3000 -t 150
	#
	#
	#	ver0.1 2019-11-12
	#	ver0.2 2019-11-14 unicycler supported, length filter added, chr2 supported
	#	ver0.3 2019-11-15 filtering paf from small alignment 
	#	ver0.4 2019-11-15 vag fix. variant call supported
	#	ver0.5 2019-11-19 duplicationa and MNP call compressed
	#	ver0.6 2019-11-20 5'end and 3'end trimming supported
	####################################################################################################################################

	print "\n\n##################################################################################################################################\n";
	print "Program: pseudo_complete_genome_builder\n";
	print "version 0.1\n\n";
	print "\tUsage\:\n";
	print "\tperl pseudo_complete_genome_builder.pl -f reference.fasta -s scaffolds.fasta -max 50 -min 40 -n 1 -p 3000 -l 3000 -t 150 -a spades\n\n";
	print "\tperl pseudo_complete_genome_builder.pl -f reference.fasta -s assembly.fasta -max 1.1 -min 0.9 -n 1 -p 3000 -l 3000 -t 150 -a unicycler\n\n";
	print "\tInput/output options:\n\n";
	print "\t -f	reference fasta (Required)\n";
	print "\t -s	contigs/scaffolds (Required)\n";
	print "\t -max	coverage max (Required)\n";
	print "\t -min	coverage minimum (Required)\n";
	print "\t -l	minimum contig length (bp) (default: 3000)\n";
	print "\t -g	minimum alignment paf (bp) (default: 3000)\n";
	print "\t -n	reference chromosome number (max 2) (default: 1)\n";
	print "\t -t	trimmng size from both ends of assembly (default: 150)\n";
	print "\t -a <spades|unicycler>	assembler (default: spades)\n";
	print "######################################################################################################################################\n\n";
	my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\t";
	print "Starts pseudo_complete_genome_builder\n\n\n";
	system("sleep 2s");

}


#-----------------------------------------------------------------------------------------------------------------------------------
sub change {
	#unicyclerはヘッダーの39とかの番号しか読み込まれないので、今のうちのにスペースを置換しておく。
	open IN, "<$sca" or die "cant open sscaffolds file!\n";
	open (OUT, '>temp/temp');
	while (my $line1 = <IN>) {
		chomp($line1);
		$line1=~ s/ /\_/g if($line1 =~ "\>");
		print OUT "$line1\n";
	}
	system("sleep 1s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub map {
	my $sum = $trimlen + $trimlen;
	system("seqkit seq -m $sum temp/temp |seqkit subseq -r $trimlen\:-$trimlen - > temp/tempseq.fa");
	system("sleep 1s");
	system("minimap2 -t 4 -cx asm5 $ref temp/tempseq.fa > temp/mapped.paf 2> temp/log");
	system("sleep 1s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub filterpaf {
	open INPUT, "<temp/mapped.paf" or die "cant open paf file1!\n";
	open (OUT, '>temp/filtered.paf');
	while (my $line = <INPUT>) {
		my @array = split(/\t/, $line);
		my $cover = abs($array[3] - $array[2]);
		next if($cover < $filter);#pafの短いアラインメントはスキップ
		print OUT "$line";
	}
	system("sleep 1s");
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub get_position {
open INPUT1, "<temp/filtered.paf" or die "cant open paf file1!\n";
open INPUT1B, "<temp/filtered.paf" or die "cant open paf file2!\n";
open (OUT, '>temp/tempfile');

print OUT "#scaffold_name\tstrand\tscaffold_start\tscaffold_end\tref_start\tref_end\tref_chr\n";
my $file = "";
my $line1B = <INPUT1B>;
while (my $line1 = <INPUT1>) {
	my $line1B = <INPUT1B>;
	chomp($line1);chomp($line1B);#改行を除く
	my @array1 = split(/\t/, $line1);
	my @array1B = split(/\t/, $line1B);
	next if($array1[0] eq $array1B[0]);
	open INPUT2, "<temp/filtered.paf" or die "cant open paf file3!\n";
	my $maxL = 1;
	my $strand = "+";
	while (my $line2 = <INPUT2>) {
		my @array2 = split(/\t/, $line2);
		next unless($array1[0] eq $array2[0]);
		$array1[7] = $array2[7] if($array2[7] < $array1[7]);#ref left most
		$array1[8] = $array2[8] if($array2[8] > $array1[8]);#ref right most
		
		$array1[2] = $array2[2] if($array2[2] < $array1[2]);#ref left most
		$array1[3] = $array2[3] if($array2[3] > $array1[3]);#ref right most
		
		#strand
		my $length = $array2[3] - $array2[2];
		$strand = $array2[4] if($maxL < $length);
		$maxL = $length if($maxL < $length);
	}
	#coverage check
	if($assembler eq "spades"){
		#NODE_17_length_6341_cov_72.379305
		my @check = split(/\_/, $array1[0]);
		next if ($check[5] > $max);
		next if ($check[5] < $min);
	}elsif($assembler eq "unicycler"){
		
		my @check = split(/\=/, $line1);
		$check[2] =~ s/\x$//g;
		next if ($check[2] > $max);
		next if ($check[2] < $min);
	}
	#length filter
	next if(abs($array1[2] - $array1[3]) < $size);
	print OUT "$array1[0]\t$strand\t$array1[2]\t$array1[3]\t$array1[7]\t$array1[8]\t$array1[5]\n";
}
system("sleep 1s");
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub sort1 {
	my $tsv_file = "temp/tempfile";
	open(IN,"$tsv_file");
	my @InFile = <IN>;
	close(IN);
	my @sortdata = sort { (split(/\t/,$a))[4] <=> (split(/\t/,$b))[4]} @InFile;

	open(OUT, '>temp/temp_sort0');
	print OUT @sortdata;
	close(OUT);
}

sub sort2 {
	open TEMP, "<temp/temp_sort0" or die "cant open temp_sort0!\n";
	open (OUT1, '>temp/temp_sort_chr1');
	open (OUT2, '>temp/temp_sort_chr2');
	my $line = <TEMP>;
	print OUT1 "$line";
	print OUT2 "$line";
	while ($line = <TEMP>) {
		chomp($line);
		my @chromosome = split(/\t/, $line);
		print OUT1 "$line\n" if($chromosome[6] eq "chr1" || $chromosome[6] eq "chromosome1");
		print OUT2 "$line\n" if($chromosome[6] eq "chr2" || $chromosome[6] eq "chromosome2");
	}
	close(OUT1);close(OUT2);
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub exort_seq_chr1 {
	open TEMP, "<temp/temp_sort_chr1" or die "cant open tempfile!\n";
	open TEMP2, "<temp/temp_sort_chr1" or die "cant open tempfile!\n";
	my $count = 0;
	my $line = <TEMP>;
	my $line2 = <TEMP2>;my $line2 = <TEMP2>;#1行先読み
	while ($line = <TEMP>) {
		$line2 = <TEMP2>;
		chomp($line);chomp($line2);
		my @array = split(/\t/, $line);
		my @array2 = split(/\t/, $line2);
		$count++;
		if($array[4] >0 && $count==1){ #refの配列から始まっているなら
			my $refstart = 1;
			my $refend = $array[4];
			system("samtools faidx $ref $array[6]\:$refstart\-$refend > temp/seq");
		}
		
		#ここが大事なところ、referenceとマッッチした"scaffoldsの"指定領域を取り出す
		my $start = $array[2];#scaffold_start
		my $end = $array[3];#scaffold_end
		system("samtools faidx temp/tempseq.fa $array[0]\:$start\-$end > temp/tempseq");
		system("sleep 0.5s");
		
		#rev compでreferenceにmapされているものは取り出したsca配列をrev compにしてから連結する。
		if($array[1] eq "-"){
			system("seqkit seq -v -pr temp/tempseq >> temp/seq");
		}elsif($array[1] eq "+"){
			system("cat temp/tempseq >> temp/seq");
		}
		system("sleep 0.5s");
		system("rm temp/tempseq");
		system("sleep 0.5s");
		
		#scaffoldの切れ目をreference配列で埋める
		my $distance = $array2[4] - $array[5];
		if($distance > 1){
			my $refstart = $array[5] + 1;
			my $refend = $array2[4];
			system("echo \>$ref\_$refstart\-$refend >> temp/seq");
			system("samtools faidx $ref $array[6]\:$refstart\-$refend >> temp/seq");
		}
		system("sleep 0.5s");
	}
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub exort_seq_chr2 {
	open TEMP, "<temp/temp_sort_chr2" or die "cant open tempfile!\n";
	open TEMP2, "<temp/temp_sort_chr2" or die "cant open tempfile!\n";
	my $count = 0;
	my $line = <TEMP>;
	my $line2 = <TEMP2>;my $line2 = <TEMP2>;#1行先読み
	while ($line = <TEMP>) {
		$line2 = <TEMP2>;
		chomp($line);chomp($line2);
		my @array = split(/\t/, $line);
		my @array2 = split(/\t/, $line2);
		$count++;
		if($array[4] >0 && $count==1){ #refの配列から始まっているなら
			my $refstart = 1;
			my $refend = $array[4];
			system("samtools faidx $ref $array[6]\:$refstart\-$refend > temp/seq2");
		}
		
		#ここが大事なところ、referenceとマッッチした"scaffoldsの"指定領域を取り出す
		my $start = $array[2];#scaffold_start
		my $end = $array[3];#scaffold_end
		system("samtools faidx temp/tempseq.fa $array[0]\:$start\-$end > temp/tempseq2");
		system("sleep 0.5s");
		
		#rev compでreferenceにmapされているものは取り出したsca配列をrev compにしてから連結する。
		if($array[1] eq "-"){
			system("seqkit seq -v -pr temp/tempseq2 >> temp/seq2");
		}elsif($array[1] eq "+"){
			system("cat temp/tempseq2 >> temp/seq2");
		}
		system("sleep 0.5s");
		system("rm temp/tempseq2");
		system("sleep 0.5s");
		
		#scaffoldの切れ目をreference配列で埋める
		my $distance = $array2[4] - $array[5];
		if($distance > 1){
			my $refstart = $array[5] + 1;
			my $refend = $array2[4];
			system("echo \>$ref\_$refstart\-$refend >> temp/seq2");
			system("samtools faidx $ref $array[6]\:$refstart\-$refend >> temp/seq2");
		}
		system("sleep 0.5s");
	}
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub concatenate_chr1 {
#fastaファイルを1行1ファイルにする。
	open INPUT3, "<temp/seq" or die "cant open seq file1!\n";
	open (OUT3, '>temp/sequence');
	my $line3 = <INPUT3>;
	print OUT3 "\>pseudo_chromosome_assembly_chr1\n" if($line3 =~ "\>");#先頭行を出力
	while ($line3 = <INPUT3>) {
		chomp($line3);
		next if($line3 =~ "\>");
		print OUT3 "$line3";
	}
	print OUT3 "\n";
	system("sleep 0.5s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub concatenate_chr2 {
	open INPUT3, "<temp/seq2" or die "cant open seq2 file. This happens when chr2 mapped contig not exit\n";
	open (OUT3, '>temp/sequence2');
	my $line3 = <INPUT3>;
	print OUT3 "\>pseudo_chromosome_assembly_chr2\n" if($line3 =~ "\>");
	while ($line3 = <INPUT3>) {
		chomp($line3);
		next if($line3 =~ "\>");
		print OUT3 "$line3";
	}
	print OUT3 "\n";
	system("sleep 0.5s");
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub seqret {
	system("seqret temp/sequence report/pseudo_chromosome_assembly_chr1.fasta");
	system("seqret temp/sequence2 report/pseudo_chromosome_assembly_chr2.fasta") if($chrnum == 2);
	system("sleep 1s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub variant_call {
	system("minimap2 -t 4 -cx asm5 --cs $ref report/pseudo_chromosome_assembly_chr1.fasta | sort -k6,6 -k8,8n -  | paftools.js call - > temp/variant_call_chr1.txt 2> report/variant_call_summary_chr1");
	system("minimap2 -t 4 -cx asm5 --cs $ref report/pseudo_chromosome_assembly_chr2.fasta | sort -k6,6 -k8,8n -  | paftools.js call - > temp/variant_call_chr2.txt 2> report/variant_call_summary_chr2") if($chrnum == 2);
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub variant_filter1 {
	open INPUT1, "<temp/variant_call_chr1.txt" or die "cant open file!\n";
	open (OUT1, '>report/variant_call_summary_chr1_filtered.txt');
	print OUT1 "\#type\tref_chr\tref_start\tref_end\tassembly_chr\tassembly_start\tassembly_end\tstrand\tref_size\talt_size\tref\talt\n";
	
	my $cycle = 0;
	while (my $line1 = <INPUT1>) {
		chomp($line1);
		my @array1 = split(/\t/, $line1);
		next if($array1[4] eq "");
		
		open INPUT2, "<temp/variant_call_chr1.txt" or die "cant open file!\n";
		my $con = 0;
		my $double = 1;
		while (my $line2 = <INPUT2>){
			chomp($line2);
			my @array2 = split(/\t/, $line2);
			my $diffL = $array2[2] - $array1[2];
			my $diffR = $array2[3] - $array1[3];
			next if($array2[2] < $array1[2]);
			
			if($diffL == 0 && $diffR == 0){
				$con++;
				next;
			}elsif($diffL == $double && $diffR == $double){
				$double++;
				$array1[6] .= $array2[6];
				$array1[7] .= $array2[7];
				next;
			}
		}
		my $reflen = length($array1[6]);
		my $altlen = length($array1[7]);
		$cycle++;
		#unique line
		if($double > 1){
			print OUT1 "MNPs\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t $reflen\t$altlen\t$array1[6]\t$array1[7]\n";
			for (my $count = 1; $count < $double; $count++){
				my $line1 = <INPUT1>;
			}
			next;
		}elsif($con == 1){
			print OUT1 "unique\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t $reflen\t$altlen\t$array1[6]\t$array1[7]\n";
			next;
		}elsif($con > 1){
		#dulicates
			print OUT1 "dups\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t $reflen\t$altlen\t$array1[6]\t$array1[7]\n";
			for (my $count = 1; $count < $con; $count++){
				my $line1 = <INPUT1>;
			}
		}
	}
}


#-----------------------------------------------------------------------------------------------------------------------------------
sub variant_filter2 {
	open INPUT1, "<temp/variant_call_chr2.txt" or die "cant open file!\n";
	open (OUT1, '>report/variant_call_summary_chr2_filtered.txt');
	print OUT1 "\#type\tref_chr\tref_start\tref_end\tassembly_chr\tassembly_start\tassembly_end\tstrand\tref_size\talt_size\tref\talt\n";
	
	my $cycle = 0;
	while (my $line1 = <INPUT1>) {
		chomp($line1);
		my @array1 = split(/\t/, $line1);
		next if($array1[4] eq "");
		
		open INPUT2, "<temp/variant_call_chr2.txt" or die "cant open file!\n";
		my $con = 0;
		my $double = 1;
		while (my $line2 = <INPUT2>){
			chomp($line2);
			my @array2 = split(/\t/, $line2);
			my $diffL = $array2[2] - $array1[2];
			my $diffR = $array2[3] - $array1[3];
			next if($array2[2] < $array1[2]);
			
			if($diffL == 0 && $diffR == 0){
				$con++;
				next;
			}elsif($diffL == $double && $diffR == $double){
				$double++;
				$array1[6] .= $array2[6];
				$array1[7] .= $array2[7];
				next;
			}
		}
		my $reflen = length($array1[6]);
		my $altlen = length($array1[7]);
		$cycle++;
		#unique line
		if($double > 1){
			print OUT1 "MNPs\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t $reflen\t$altlen\t$array1[6]\t$array1[7]\n";
			for (my $count = 1; $count < $double; $count++){
				my $line1 = <INPUT1>;
			}
			next;
		}elsif($con == 1){
			print OUT1 "unique\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t $reflen\t$altlen\t$array1[6]\t$array1[7]\n";
			next;
		}elsif($con > 1){
		#dulicates
			print OUT1 "dups\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t $reflen\t$altlen\t$array1[6]\t$array1[7]\n";
			for (my $count = 1; $count < $con; $count++){
				my $line1 = <INPUT1>;
			}
		}
	}
}
