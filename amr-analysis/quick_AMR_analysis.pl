# Original location: /gpfs/project/dilthey/projects/KI/AMR_analysis/quick_AMR_analysis.pl

use strict;
use Data::Dumper;

my @directories = grep {-d $_} glob('*');

my %all_elements;
my %elements_per_sample;
my %reads_and_bases_per_sample;
foreach my $directory (@directories)
{
	die unless($directory =~ /_MinION/);
	(my $sampleID = $directory) =~ s/_MinION//;
	my $AMR_file = $directory . '/AMR_output_identifiedElements.txt';

	{
		die unless(-d '../combined_fastq');
		
		my $FASTQ_file;
		
		(my $sampleID_d = $sampleID) =~ s/KI/KID/;
		(my $sampleID_d_no11 = $sampleID_d) =~ s/A11/A/;
		
		if(-e '../combined_fastq/' . $sampleID . '_MinION_combined_fixed.fastq')
		{
			$FASTQ_file = '../combined_fastq/' . $sampleID . '_MinION_combined_fixed.fastq';
		}
		elsif(-e '../combined_fastq/' . $sampleID . '_MinION_combined.fastq')
		{
			$FASTQ_file = '../combined_fastq/' . $sampleID . '_MinION_combined.fastq';
		}
		elsif(-e '../combined_fastq/' . $sampleID_d . '_MinION_combined.fastq')
		{
			$FASTQ_file = '../combined_fastq/' . $sampleID_d . '_MinION_combined.fastq';
		}
		elsif(-e '../combined_fastq/' . $sampleID_d_no11 . '_MinION_combined.fastq')
		{
			$FASTQ_file = '../combined_fastq/' . $sampleID_d_no11 . '_MinION_combined.fastq';
		}	
		else
		{
			die "Can't find FASTQ for $sampleID";
		}	

		die "File $FASTQ_file not existing" unless(-e $FASTQ_file);
		$reads_and_bases_per_sample{$sampleID} = count_reads_data_in_fastq($FASTQ_file);
	}
	
	if(-e $AMR_file)
	{
		open(AMR, '<', $AMR_file) or die "Cannot open $AMR_file";
		while(<AMR>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			die unless($line_fields[2]);
			my $AMR_element = $line_fields[2];
			$all_elements{$AMR_element}++;
			$elements_per_sample{$sampleID}{$AMR_element}++;
			
			$all_elements{'Total'}++;
			$elements_per_sample{$sampleID}{'Total'}++;			
		}
		close(AMR);
	}
	else
	{
		warn "No AMR file for $sampleID";
	}
}

my @AMR_elements = sort {$a cmp $b} keys %all_elements;
my $output_fn = 'combined_output.txt';
open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
print OUTPUT join("\t", 'sampleID', 'nReads', 'nBases', @AMR_elements, (map {$_ . '_propReads'} @AMR_elements)), "\n";
foreach my $sampleID (keys %elements_per_sample)
{
	my @output_AMR;
	my @output_AMR_rescaled;
	if(exists $elements_per_sample{$sampleID})
	{
		foreach my $AMRid (@AMR_elements)
		{
			if(exists $elements_per_sample{$sampleID}{$AMRid})
			{
				my $nReads = $elements_per_sample{$sampleID}{$AMRid};
				my $nAllReads = $reads_and_bases_per_sample{$sampleID}[0];
				die unless($nAllReads);
				push(@output_AMR, $nReads);
				push(@output_AMR_rescaled, $nReads / $nAllReads);
			}
			else
			{
				push(@output_AMR, 0);
				push(@output_AMR_rescaled, 0);
			}
		}
	}
	else
	{
		@output_AMR = ('NA' x scalar(@AMR_elements));
		@output_AMR_rescaled = ('NA' x scalar(@AMR_elements));
	}
	die unless(scalar(@output_AMR) == scalar(@AMR_elements));
	print OUTPUT join("\t", $sampleID, $reads_and_bases_per_sample{$sampleID}[0], $reads_and_bases_per_sample{$sampleID}[1], @output_AMR, @output_AMR_rescaled), "\n";
}
close(OUTPUT);

sub count_reads_data_in_fastq
{
	my $in = shift;
	
	my $n_reads = 0;
	my $n_bases = 0;
	die unless(defined $in);
	open(FASTQ, '<', $in) or die "Cannot open $in";
	while(<FASTQ>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die unless(substr($line, 0, 1) eq '@');
		my $sequence = <FASTQ>; chomp($sequence);
		my $plus = <FASTQ>; die unless(substr($plus, 0, 1) eq '+');
		my $qual = <FASTQ>; chomp($qual);
		unless(length($sequence) == length($qual))
		{
			warn "Observed corruption in file $in";
		}
		
		$n_reads++;
		$n_bases += length($sequence);
	}
	close(FASTQ);
	
	return [$n_reads, $n_bases];
}

