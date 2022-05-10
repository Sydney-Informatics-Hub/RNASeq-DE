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

# Create matrix of TPM counts from TPMCalculator output
# Compatible with 0.0.4 release - output file formats have changed over releases
# Output will contain gene level TPMs
# TPMCalculator does also give TPM at the transcript, exon/intron level
# rows = genes, columns = sample
# 
# Takes all files from a directory named <sampleid>.final_genes.out
# For lots of samples (>100), need more memory, run in job

# Collect sampleids from all files in a directory
# Files must be named <sampleid>.counts

my $tpmdir=$ARGV[0];

if($tpmdir){
	$tpmdir=~s/\/$//;
}
else{
	print "Please provide the path to the TPMCalculator output directory, e.g. perl tpmcalculator_make_matrix.pl ../Cohort_TPMCalculator\n";
	exit;
}
my $out="$tpmdir\/TPM_GeneLevel.counts";

# Loop through datasets, find count data
my @sampleids=`(find $tpmdir -not -empty -type f -name "*final_genes.out" -printf "%f\n" | sed 's/.final_genes.out//')`;
chomp(@sampleids);
my @empty=`(find $tpmdir -empty -type f -name "*final_genes.out" -printf "%f\n" | sed 's/.final_genes.out//')`;
chomp(@empty);

if(@empty){
	print "WARNING: Empty files found: @empty\n";
}

print "Saving gene level TPMs into memory...\n";
my $genecountshash = {};
my $i=0;
my @allgenes; my $gene; my $count;
foreach my $sampleid (@sampleids){
        my $file="$tpmdir\/$sampleid\.final_genes.out";
        open(FILE,$file)||die "Could not open $file: $!\n";
        my $header=<FILE>;
	while(<FILE>){
		chomp;
                my(@col)=split(" ",$_);
		my $gene=$col[0];
		my $tpm=$col[6];
		if($i==1){
			# Save gene IDs just once in the same order of the first sample (should be same for all samples)
			push @allgenes, $gene;
		}
		$genecountshash->{$sampleid}->{$gene}->{tpm}=$tpm;
	}
        $i++;
}
my $num_samples=scalar(@sampleids);
print "Printing TPM counts for $num_samples samples to $out\n";

# Print to file
open(O,">$out")|| die "Could not print to $out: $!\n";

# Print headers
print O "Gene\t";
foreach my $sampleid (sort keys %{$genecountshash} ){
	print O "$sampleid\t";
}
print O "\n";

# Print counts
foreach my $gene (@allgenes){
	print O "$gene\t";
	foreach my $sampleid (sort keys %{$genecountshash} ){
		if($genecountshash->{$sampleid}->{$gene}->{tpm}){
			print O "$genecountshash->{$sampleid}->{$gene}->{tpm}\t";
		}
		else{
			# No record for the gene & sample
			print O "0\t";
		}
	}
	print O "\n";
}
exit;

