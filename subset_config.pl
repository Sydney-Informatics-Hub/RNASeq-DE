#! /usr/bin/env perl
#
use strict;
use warnings;

use List::Util 'shuffle';
#
# Create subsets of samples.config containing different number of samples
# Samples is the typical unit of processing (rather than files)
# Resulting config files can be used for benchmarking and demonstrating scalability of jobs
#
my $config="../samples.config_benchmarking";

open(CONFIG,"$config")||die "Could not open $config: $!\n";
my $samplehash = {};
while(<CONFIG>){
	chomp;
	if($_!~/^#/){
		my($fastq, $sampleid, $dataset, $ref, $seqcentre, $platform, $runtype, $lib) = split(" ",$_);
		$samplehash->{$sampleid}->{1}=1;
	}
}

my @samplehash_keys = keys %{$samplehash};
my $size = keys %{$samplehash};
print "Number of samples in $config: $size\n";

# Create 5 config files with: 20, 40, 60, 80, 100 samples each
for(my $i=20; $i<=100; $i+=20){
	my $out="../$i\_samples.config";
	`rm -rf $out`;
	
	# Shuffle numbers 0-number of samples to get a random selection of indexes
	# Take first i values and use as indexes to extract random i samples
	my @shuffle= shuffle( 0..$#samplehash_keys);
	my @indexes=@shuffle[ 0 .. $i-1 ];

	# Check there are the correct number of samples
	my $length=scalar(@indexes);
	print "Writing to $out for $length samples\n";

	open(OUT,">$out")||die "Could not write to $out: $!\n";
	# Print config headers
	my $header=`head -1 $config`;
	print OUT "$header";
	foreach my $index (@indexes){
		my $line=`awk -v pattern="$samplehash_keys[$index]" '\$2 ~ pattern{print \$0}' $config`;
		chomp($line);
		print OUT "$line\n"
	}
	close OUT;
}
exit;

