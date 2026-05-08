#!/usr/bin/perl

# exract the allele frequencies and counts per chromosome positions for species list

use warnings;
use strict;


my $vcfPrefix = $ARGV[0];
my $posDir = "/rds/project/rds-8b3VcZwY7rY/projects/cichlid/Bettina/data/allele_frequency_for_eggspot_candidate_snps_per_species/posFiles/";
my $speciesDir ="/rds/project/rds-8b3VcZwY7rY/projects/cichlid/Bettina/data/allele_frequency_for_eggspot_candidate_snps_per_species/speciesLists_callset_v3/";
my @chromosomes = (1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23);
my @speciesFile = glob('/rds/project/rds-8b3VcZwY7rY/projects/cichlid/Bettina/data/allele_frequency_for_eggspot_candidate_snps_per_species/speciesLists_callset_v3/*');


if (!@ARGV) {
	print"perl get_AF_per_chromosomes.pl <vcfPrefix>\n";
	print "\tvcfPrefix      path and prefix without the .chrN.bcf.gz\n";
	exit;
}



### process each species
foreach my $popFile (@speciesFile) {
    
    my $species = "NA";
    if($popFile =~ m/^$speciesDir(\S+)_individuals.txt$/) {
        $species = $1;
    }
    
    if($species ne "NA" & $species ne "callset_v3") {
        print "$species\n";
 
        my $resFile = "callset_v3_results/$species"."_AF_result.txt";
        
        open(OUT, ">", $resFile) or die $!;
        print OUT "CHROM\tPOS\tREF\tALT\tref|het|alt\tAF\n";
        
        ### iterate over the chromosomes
        foreach my $n (@chromosomes) {
    
            my $posFile = $posDir."chr".$n."_posFile.txt";
            my $vcfFile = "$vcfPrefix"."$n".".bcf.gz";
#            my $vcfFile = "$vcfPrefix"."$n".".vcf.gz";
#            print "$posFile\n$vcfFile\n\n";
    
	        my $bcf_command = "bcftools view -R $posFile -S $popFile $vcfFile | bcftools view -a | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed -e 's/\\//|/g'";

	
        	open(my $fh, '-|', $bcf_command) or die $!;
        	while (my $line = <$fh>) {
        		chomp $line;
        		countAndPrint($line);
        	}
        	close $fh;
	

        }
    }
    
    close OUT;
    
}



###### subroutines

sub countAndPrint {
	my ($line) = @_;
	my @tabs = split(/\t/,$line);
	my @GT;
	my $ref_cnt = 0;
	my $het_cnt = 0;
	my $alt_cnt = 0;
	my $chrom = shift(@tabs);
	my $pos = shift(@tabs);
	my $ref = shift(@tabs);
	my $alt = shift(@tabs);
	my $af = 0;
	my $not_biallelelic_flag = 0;
	# check if site is biallelic
	if($alt =~ m/[,]/) {
	    $not_biallelelic_flag = 1;
	    print OUT "$chrom\t$pos\t$ref\t$alt\tnot_biallelic\n";
	}
	if($not_biallelelic_flag == 0) {
    	# loop over sample columns and count
    	foreach my $GT (@tabs) {
    		if($GT =~ m/^0\|0$/) {
    			$ref_cnt++;
    			$GT = 0;
    			push(@GT, $GT);
    		}
    		if($GT =~ m/^0\|\.$/) {
    			$ref_cnt++;
    			$GT = 0;
    			push(@GT, $GT);
    		}
    		if($GT =~ m/^0\|1$/) {
    			$het_cnt++;
    			$GT = 1;
    			push(@GT, $GT);
    		}
    		if($GT =~ m/^1\|0$/) {
    			$het_cnt++;
    			$GT = 1;
    			push(@GT, $GT);
    		}
    		if($GT =~ /^1\|1$/) {
    			$alt_cnt++;
    			$GT = 2;
    			push(@GT, $GT);
    		}
    		if($GT =~ /^\.\|1$/) {
    			$alt_cnt++;
    			$GT = 2;
    			push(@GT, $GT);
    		}
    		if($GT =~ /^\.|\.$/) {
    			$GT = "NA";
    			push(@GT, $GT);
    		}
    	}
    	my $gt_string = join("\t", @GT);
    	
    	my $total = ($ref_cnt + $het_cnt + $alt_cnt) * 2;
    	if($total > 0) {
    	    $af = ($het_cnt + ($alt_cnt * 2)) / $total;
    	    $af = sprintf("%.4f", $af);
    	}
    #	print "$pos\t$ref\t$alt\t$ref_cnt|$het_cnt|$alt_cnt\t$af\t$gt_string\n";
    	print OUT "$chrom\t$pos\t$ref\t$alt\t$ref_cnt|$het_cnt|$alt_cnt\t$af\n";
	}
}
