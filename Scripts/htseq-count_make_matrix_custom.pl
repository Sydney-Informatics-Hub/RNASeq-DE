#! /usr/bin/perl

use strict;
use warnings;
use File::Basename;

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
# All files ../$cohort\_htseq-count/*.counts are included
# For lots of samples (>100), need more memory, run in a job
#
# Collect sampleids from all files in a directory
# Files must be named <sampleid>.counts

my $config=basename($ARGV[0]);
(my $cohort = $config) =~s/\.config//;
my $countdir="../$cohort\_htseq-count";
my $out="$countdir\/$cohort\.counts";

# Loop through datasets, find count data
my @sampleids=`(find $countdir -not -empty -type f -printf "%f\n" -name "*counts" | sed 's/.counts//')`;
chomp(@sampleids);
my @empty=`(find $countdir -empty -type f -printf "%f\n" -name "*counts" | sed 's/.counts//')`;
chomp(@empty);

if(@empty){
	print "WARNING: Empty htseq-count files found: @empty\n";
}

my $genecountshash = {};
my $i=0;
my @allgenes; my $gene; my $count;
foreach my $sampleid (@sampleids){
        my $file="$countdir\/$sampleid\.counts";
        open(FILE,$file)||die "Could not open $file: $!\n";
        
	while(<FILE>){
		chomp;
                my($gene,$count)=split(" ",$_);
		if($i==1){
			# Save gene IDs just once in the same order of the first sample (should be same for all samples)
			push @allgenes, $gene;
		}
		$genecountshash->{$sampleid}->{$gene}->{count}=$count;
	}
        $i++;
}

# Print to file
open(O,">$out")|| die "Could not print to $out: $!\n";

# Print headers
print O "GeneName\t";
foreach my $sampleid (sort keys %{$genecountshash} ){
	print O "$sampleid\t";
}
print O "\n";

# Print counts
foreach my $gene (@allgenes){
	print O "$gene\t";
	foreach my $sampleid (sort keys %{$genecountshash} ){
		print O "$genecountshash->{$sampleid}->{$gene}->{count}\t";
	}
	print O "\n";
}
exit;
