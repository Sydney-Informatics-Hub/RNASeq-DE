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

# Collate all mapping information from STAR output for each aligned file 
# perl summarise_STAR_alignment_stats.pl ../cohort.config
# Output: ../QC_reports/cohort_STAR_metrics.txt

my $config=basename($ARGV[0]);
(my $cohort = $config) =~s/\.config//;
my @logs;
# Collect sampleIDs and find STAR Log.final.out files
# ../${dataset}_STAR/${sampleid}_${lane}_Log.final.out
open(CONFIG,$ARGV[0])||die "$!\n";
while(<CONFIG>){
	chomp;
	if($_!~/^#/){
		my($fastq, $sampleid, $dataset, @other)=split(" ",$_);
		my $indir="../$dataset\_STAR";
		my @datasetlogs=`ls $indir\/*Log.final.out`;
		push(@logs, @datasetlogs);
	}
}

# Most inefficient way ever not not worth fixing now!
# Still runs pretty fast
my @starlogs = do { my %seen; grep { !$seen{$_}++ } @logs };
my $num_logs = scalar(@starlogs);
print "Number of STAR logs found for: $num_logs\n";
my $output="../QC_reports\/$cohort\_STAR_metrics.txt";

open(O,">$output")||die "Could not write to $output: $!\n";
print O "#File\tNumber_of_input_reads\tAverage_read_length\t";
print O "Uniquely_mapped_reads\tUniquely_mapped_reads_%\tAverage_mapped_length\t";
print O "Number_of_splices:_Total\tNumber_of_splices:_Annotated_(sjdb)\tNumber_of_splices:_GT/AG\tNumber_of_splices:_GC/AG\tNumber_of_splices:_AT/AC\tNumber_of_splices:_Non-canonical\t";
print O "Mismatch_rate_per_base_%\tDeletion_rate_per_base\tDeletion_average_length\tInsertion_rate_per_base\tInsertion_average_length\t";
print O "Reads_mapped_to_multiple_loci\tReads_mapped_to_multiple_loci_%\tReads_mapped_to_too_many_loci\tReads_mapped_to_too_many_loci_%\tReads_unmapped_%:_too_many_mismatches\tReads_unmapped_%:_too_short\tReads_unmapped_%:_other\t";
print O "Chimeric_reads\tChimeric_reads_%\n";

foreach my $file (@starlogs){
        chomp $file;
        print O "$file\t";
        open(F,$file)||die "Could not open $file: $!\n";
        while(<F>){
                chomp;
                if($_=~/(\s)+(.*)(\s)\|(\s)+(.*)/){
                        my $header=$2;
                        my $value=$5;
                        if($header=~/Number of input reads/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Average input read length/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Uniquely mapped reads number/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Uniquely mapped reads %/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Average mapped length/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of splices: Total/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of splices: Annotated \(sjdb\)/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of splices: GT\/AG/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of splices: GC\/AG/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of splices: AT\/AC/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of splices: Non-canonical/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Mismatch rate per base\, %/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Deletion rate per base/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Deletion average length/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Insertion rate per base/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Insertion average length/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of reads mapped to multiple loci/i){
                                print O "$value\t";
                        }
                        elsif($header=~/% of reads mapped to multiple loci/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of reads mapped to too many loci/i){
                                print O "$value\t";
                        }
                        elsif($header=~/% of reads mapped to too many loci/i){
                                print O "$value\t";
                        }
                        elsif($header=~/% of reads unmapped: too many mismatches/i){
                                print O "$value\t";
                        }
                        elsif($header=~/% of reads unmapped: too short/i){
                                print O "$value\t";
                        }
                        elsif($header=~/% of reads unmapped: other/i){
                                print O "$value\t";
                        }
                        elsif($header=~/Number of chimeric reads/i){
                                print O "$value\t";
                        }
                        elsif($header=~/% of chimeric reads/i){
                                print O "$value\n";
                        }
                }
        }
}

print "Printed output to: $output\n";
