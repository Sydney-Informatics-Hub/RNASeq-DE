#! /usr/bin/perl

use strict;
use warnings;

#########################################################
#
# Platform: NCI Gadi HPC
#
# Author: Tracy Chew
# tracy.chew@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the scientific and technical assistance
# <or e.g. bioinformatics assistance of <PERSON>> of Sydney Informatics
# Hub and resources and services from the National Computational
# Infrastructure (NCI), which is supported by the Australian Government
# with access facilitated by the University of Sydney.
#
#########################################################

# Create matrix of counts for DESeq2
# rows = genes, columns = sample
# Count data is from HTSeq
# Only ../$dataset\_htseq-count/<sampleid>.counts are included
# sampleids are from column 2 of the config file
# Take list of dataset and sample IDs so that matrix can be customised
# For lots of samples (>100), need more memory, run as job

my $config=$ARGV[0];
(my $cohort = $config) =~s/\.config//;
open(CONFIG,"$config")||die "Could not open $config: $!\n";
my $samplehash;
my @wantlist;
my @datasets;
while(<CONFIG>){
        chomp;
        if($_!~/^#/){
                my($fastq, $sampleid, $dataset, $ref, $seqcentre, $platform, $runtype, $lib) = split(" ",$_);
		#print "Looking for $sampleid from $dataset\n";
		$samplehash->{$dataset}->{$sampleid}->{1}=1;
		push @wantlist, $sampleid;
		push @datasets, $dataset;
        }
}

# Get list of gene ids

# Loop through datasets, find count data
my $genecountshash;
my @allgenes;
my @empty;
my $i=1;
my @uniq_datasets=uniq(@datasets);
foreach my $dataset(@uniq_datasets){
	my $countpath="../$dataset\_htseq-count";
	my @setids=`(find $countpath -not -empty -printf "%f\n" -name "*counts" | sed 's/.counts//')`;
        @empty=`(find $countpath -empty -printf "%f\n" -name "*counts" | sed 's/.counts//')`;
	#chomp @setids;
       	#print @setids;
        #print scalar @setids;
	
	## NEED TO INCLUDE IF STATEMENT TO CHECK IF @WANTLIST EXISTS IN @SETIDS
        foreach my $setids (@setids){
		chomp $setids;

                if($samplehash->{$dataset}->{$setids}->{1}){
                        my $file="$countpath\/$setids\.counts";
			#print "Including counts from $file\n";
                        #my $geneids=`awk '{print \$2}' $file`;
                        #my $samplecount=`awk '{print \$2}' $file`;
                        open(FILE,$file)||die "Could not open $file: $!\n";
                        while(<FILE>){
                                chomp;
                                my($gene,$count)=split(" ",$_);
				if($i==1){
                                        #Save gene ids only once
                                        push @allgenes, $gene;
                                }
                                $genecountshash->{$setids}->{$gene}->{count}=$count;
                        }
                        $i++;
                }
        }
        $i++;
}
chomp(@empty);

print "Empty count files: @empty\n";
# Remove samples with empty files from want list
my %seen = ();
my @wantlist_filtered;

foreach my $empty (@empty) { $seen{$empty} = 1};

# Print count file matrix
# Print headers. Note R doesn't like rows starting with #
# Whilst printing headers, check and include only count files that are not empty
my $output="$cohort\_counts.txt";
print "Writing counts to $output\n";
open(O,">$output")||die "Could not write to $output: $!\n";
print O "GeneName\t";

# Find elements in @wantlist that are not in @empty
foreach my $want (@wantlist){
	if($seen{$want}){
        	#print "$want is empty\n";
	}
	else{
		#print "$want is not empty\n";
		push(@wantlist_filtered, $want);
		print O "$want\t";
	}
}
print O "\n";

# Print count data
foreach my $ensembl(@allgenes){
        print O "$ensembl\t";
        foreach my $want (@wantlist_filtered){
                print O "$genecountshash->{$want}->{$ensembl}->{count}\t";
        }
        print O "\n";
}
exit;

## Subroutines
#
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
