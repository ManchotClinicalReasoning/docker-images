use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;
use Getopt::Long;
use File::Basename;

my $FASTQ;
my $sampleID;
my $report; #PDF containing some info
my $histogram; #Histogram
GetOptions (
	'FASTQ:s' => \$FASTQ, 
	'sampleID:s' => \$sampleID, 
	'report:s' => \$report, 
	'histogram:s' => \$histogram, 
);

unless($FASTQ)
{
	die "Please provide parameter --FASTQ";
}

unless($sampleID)
{
	$sampleID = basename($FASTQ);
	if($sampleID =~ /\./)
	{
		$sampleID =~ s/^(.+)\.(.+?)$/\1/;
	}
	$sampleID =~ s/\s//g;
	
	print STDERR "Using '$sampleID' for sample ID. If this is not desired, call me with --sampleID\n";
}

my $rl_hist_href = getReadLengthHistogram($FASTQ);


open(HIST, '>', $histogram) or die "Cannot open $histogram";
foreach my $readLength (keys %$rl_hist_href)
{
	print HIST join("\t", $readLength, $rl_hist_href->{$readLength}), "\n";
	# $np_total_bases += ($readLength * $rl_hist_href->{$readLength});
}
close(HIST);

my $cmd = qq(Rscript /scripts/plotRunStats_oneSample.R $histogram $report);
system($cmd) and die "Command $cmd failed";

print "\nDone. Produced $report\n\n";

sub getReadLengthHistogram
{
	my $file = shift;	
	
	my %hist;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{		
		my $readID = $_;
		chomp($readID);
		$readID =~ s/[\n\r]//g;
		die unless(substr($readID, 0, 1) eq '@');
		
		my $sequence = <F>;
		chomp($sequence);
		
		my $intermediate = <F>;
		die unless(substr($intermediate, 0, 1) eq '+');
		
		my $qualities = <F>;
		chomp($qualities);
		warn "FASTQ error, line $. of $file - qualities string has length " . length($qualities) . " v/s " . length($sequence) . " for sequence." unless(length($qualities) == length($sequence));
		
		my $readLength = length($sequence);
		$hist{$readLength}++;
	}	
	close(F);
		
	return \%hist;
}
