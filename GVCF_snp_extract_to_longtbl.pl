#!/usr/bin/perl
use strict;
#use warnings;

# open a text file with the snp calls
my $input = $ARGV[0];


open (VCF, "<$input") or die$!;


#say POS "Position,Quality";

print "Pos,Ref,Alt,ID,Chain,Mouse_No,DP,Allele_Ratio,SNP_type\n";

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
		#print "My data are $pos\t$ref\t$alt\n";
		#For each samples, take each field and assign DP, AD and push allele call to array in a hash
		#print POS "$pos,$qual\n";
		foreach my $k (10..$length2){
			my $l = $k - 1;
			my $n = $k - 10;
			my @stats = split(":",$fields[$l]);
			my $AR;
			my $type;
			my $ID = $names[$n];
			my @chain = split("P",$ID);
			
			#Assign a value to the DP and save it in an array called @DPs
			my $DP = $stats[2];
			my @AD = split(",",$stats[1]);

			
			#Based on Allelic ratio, designate SNPs as Reference (R) Het (H) or Homozygous (Z)
			
			if ($DP < 20 ){
			$AR = 0.000;
			$type = "-";
			}
			
			elsif ($DP >= 20 ){
			$AR = ($AD[1]) / $DP;
			$AR = sprintf "%.3f", $AR;
			
				if ($AR < 0.1){
					$type = "R";
				} elsif ($AR > 0.9){
					$type = "Z";
				}else{
				$type = "H";
				}				
			} 
			print "$pos,$ref,$alt,$ID,$chain[0],$chain[1],$DP,$AR,$type\n";
			
		}
		#print "$pos\t", join("\t", @DPs),"\n";
		#print "$pos\t", join("\t", @major),"\n";
		#print "$pos\t", join("\t", @AR),"\n";
		#say MAJ "$pos\t", join("\t", @major);
		#say ARAT "$pos,", join(",", @AR);
	}
}

close VCF;

	 