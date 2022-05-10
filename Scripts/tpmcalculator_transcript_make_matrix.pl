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
# This version contains ensemble gene ID, ensemble transcript ID followed by
# per sample TPMs
# Compatible with 0.0.4 release - output file formats have changed over releases
# TPMCalculator does also give TPM at the transcript, exon/intron level
# rows = genes, columns = sample
# 
# Takes all files from a directory named <sampleid>.final_transcripts.out
# For lots of samples (>100), need more memory, run in job

my $tpmdir=$ARGV[0];

if($tpmdir){
	$tpmdir=~s/\/$//;
}
else{
	print "Please provide the path to the TPMCalculator output directory, e.g. perl tpmcalculator_transcript_make_matrix.pl ../Cohort_TPMCalculator\n";
	exit;
}
my $out="$tpmdir\/TPM_TranscriptLevel.counts";

# Loop through datasets, find count data
my @sampleids=`(find $tpmdir -not -empty -type f -name "*final_transcripts.out" -printf "%f\n" | sed 's/.final_transcripts.out//')`;
chomp(@sampleids);
my @empty=`(find $tpmdir -empty -type f -name "*final_transcripts.out" -printf "%f\n" | sed 's/.final_transcripts.out//')`;
chomp(@empty);

if(@empty){
	print "WARNING: Empty files found: @empty\n";
}

print "Saving transcript level TPMs into memory...\n";
my $transcriptcountshash = {};
my $geneidhash = {};
my $i=0;
my @alltranscripts; my $gene; my $transcript;  my $count;
foreach my $sampleid (@sampleids){
        my $file="$tpmdir\/$sampleid\.final_transcripts.out";
	open(FILE,$file)||die "Could not open $file: $!\n";
        my $header=<FILE>;
	while(<FILE>){
		chomp;
                my(@col)=split(" ",$_);
		my $gene=$col[0];
		my $transcript=$col[1];
		my $tpm=$col[7];
		if($i==1){
			# Save transcript IDs just once in the same order of the first sample (should be same for all samples)
			push @alltranscripts, $transcript;
		}
		$transcriptcountshash->{$sampleid}->{$transcript}->{tpm}=$tpm;
		# Save ensembl gene ID for each transcript ID
		$geneidhash->{$transcript}->{geneid}=$gene;
	}
        $i++;
}

my $num_samples=scalar(@sampleids);
print "Printing TPM counts for $num_samples samples to $out\n";

# Print to file
open(O,">$out")|| die "Could not print to $out: $!\n";

# Print headers
print O "Gene\tTranscript\t";
foreach my $sampleid (sort keys %{$transcriptcountshash} ){
	print O "$sampleid\t";
}
print O "\n";

# Print counts
foreach my $transcript (@alltranscripts){
	my $gene = $geneidhash->{$transcript}->{geneid};
	print O "$gene\t$transcript\t";
	foreach my $sampleid (sort keys %{$transcriptcountshash} ){
		if($transcriptcountshash->{$sampleid}->{$transcript}->{tpm}){
			print O "$transcriptcountshash->{$sampleid}->{$transcript}->{tpm}\t";
		}
		else{
			# No record for transcript & sample
			print O "0\t";
		}
	}
	print O "\n";
}
exit;

