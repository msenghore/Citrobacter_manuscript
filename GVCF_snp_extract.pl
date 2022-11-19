#!/usr/bin/perl
use strict;
#use warnings;

# open a text file with the snp calls
my $input = $ARGV[0];
my $prefix = $ARGV[1];
my $pos_file = $prefix.".snp_positions.txt";
my $maj_file = $prefix.".snp_table.txt";
my $aratio_file = $prefix.".alt_allele_frequency.txt";

open (VCF, "<$input") or die$!;
#open (POS,">> $pos_file");
open (MAJ,"> $maj_file");
open (ARAT,"> $aratio_file");

#say POS "Position,Quality";

my @vcf=<VCF>;
close VCF or die$!;
my @names;
#Parse the list file and take each line as a variable
foreach my $line (@vcf) {
	chomp $line;
#ignore information lines at top
	#next if /^##/;
#Use header line to store sample names and print them
	
	if 	($line =~ m/CHROM/){
		my @headers=split("\t",$line);
		my $length = @headers;
		#print "My vcf has $length header fields\n";
		#Create a loop through the line to place each sample name in a array (9 info fields)
		my @samples;
		foreach my $i (10..$length) {
			my $j = $i - 1;
			my $k = $i - 9;
			my $name = $headers[$j];
			#print "Sample $k is $name\n";
			#Push the sample ids into the array @names
			push(@names,$name);	
		} 
	#print "pos\t", join("\t", @names),"\n";
	say MAJ "pos\t", join("\t", @names);
	say ARAT "pos,", join(",", @names);
	}
	#Parse remaining lines of code assign information fields to variables
	if 	($line =~ m/#/){
	}else{
		my @fields=split("\t",$line);
		my $length2 = @fields;
		my $qual = $fields[5];
		my $ref = $fields[3];
		my $alt = $fields[4];
		my $pos = $fields[1];
		my @DPs;
		my @allele_count;
		my @major;
		my @minor;
		my @AR;
		#print "My data are $pos\t$ref\t$alt\n";
		#For each samples, take each field and assign DP, AD and push allele call to array in a hash
		#print POS "$pos,$qual\n";
		foreach my $k (10..$length2){
			my $l = $k - 1;
			my $n = $k - 10;
			my @stats = split(":",$fields[$l]);
			my $call;
			my $allele_ratio;
			my $min_base;
			my $maj_base;
			my $AR1;
			my $AR2;
			my $maj_base_ratio;
			
			#Assign a value to the DP and save it in an array called @DPs
			my $DP = $stats[2];
			push( @DPs,$DP);
			
			#Extract data on Allelic depth and choose most common allele
			my @AD = split(",",$stats[1]);
			my $allele_no = @AD;
			push (@allele_count,$allele_no);
			
			#Identify the majority allele and assign the ref or alt base, assign the allele ratio
			
			if ($DP > 0 ){
			# Step one filter return N if DP is less than 100
			#$AR1 = $AD[0] / $DP;
			$AR2 = ($AD[1]) / $DP;
			#$AR1 = sprintf "%.3f", $AR1;
			$AR2 = sprintf "%.3f", $AR2;
			push (@AR,$AR2);
			print "$pos,$names[$n],$AR2\n";
			}
			
			if ($DP >= 20 ){
				if ($AR2 < 0.5){
					$maj_base = $ref;
					$min_base = $alt;
					#push (@AR,$AR1);
				} elsif ($AR2 > 0.75){
					$maj_base = $alt;
					$min_base = $ref;
					#push (@AR,$AR2);
				}else{$maj_base ="N"}
				push (@major,$maj_base);
				
			} else{ 
				push (@major,"-");
				#push (@AR,"0");
				
			}
			
		}
		#print "$pos\t", join("\t", @DPs),"\n";
		#print "$pos\t", join("\t", @major),"\n";
		#print "$pos\t", join("\t", @AR),"\n";
		say MAJ "$pos\t", join("\t", @major);
		say ARAT "$pos,", join(",", @AR);
	}
}
close POS;
close VCF;
close MAJ;
close ARAT;	
	 