# original location: /gpfs/project/dilthey/software/MetaMaps_fresh/removeHumanReads.pl

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use Getopt::Long;

my $humanGenomeFASTA;
my $inputFASTQ;
my $outputFASTQ;
my $maxmemory = '';
my $threads = 1;
my $minimap2_bin = 'minimap2';
GetOptions (
	'humanGenomeFASTA:s' => \$humanGenomeFASTA, 
	'inputFASTQ:s' => \$inputFASTQ, 
	'outputFASTQ:s' => \$outputFASTQ, 
	'maxmemory:s' => \$maxmemory, 
	'threads:s' => \$threads, 
);

unless($humanGenomeFASTA and $inputFASTQ and $outputFASTQ)
{
	print_help();
}

print "File $humanGenomeFASTA not found" unless(-e $humanGenomeFASTA);
print "File $inputFASTQ not found" unless(-e $inputFASTQ);

my $minimap2_output_file = $outputFASTQ . '.mappings';
$maxmemory = ($maxmemory) ? "--maxmemory $maxmemory " : '';
$threads = ($threads) ? "-t $threads " : '';

my $minimap2_cmd = qq($minimap2_bin $threads -c -x map-ont $humanGenomeFASTA $inputFASTQ > $minimap2_output_file);
system($minimap2_cmd) and die "Cannot execute: $minimap2_cmd";
 
my %removeReadIDs; 
open(PAF, '<', $minimap2_output_file) or die "Could not opne $minimap2_output_file";
while(<PAF>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	die unless(scalar(@line_fields) > 2);	
	my $alignment_length = $line_fields[3] - $line_fields[2];
	die unless($alignment_length);
	my $alignment_prop = $alignment_length / $line_fields[1];
	die unless(($alignment_prop >= 0) and ($alignment_prop <= 1));
	my $identity = $line_fields[9] / $line_fields[10];
	die unless(($identity >= 0) and ($identity <= 1));
	
	if(($alignment_prop >= 0.8) and ($identity >= 0.8))
	{
		my $readID = $line_fields[0];
		$removeReadIDs{$readID}++;		
	}
	
}
close(PAF);

print "Now remove " . scalar(keys %removeReadIDs) . " reads identified as human.\n";

my $reads_notRemoved = 0;
open(FASTQIN, '<', $inputFASTQ) or die "Cannot open $inputFASTQ";
open(FASTQOUT, '>', $outputFASTQ) or die "Cannot open $outputFASTQ";
while(<FASTQIN>)
{
	my $readID_line = $_;
	chomp($readID_line);
	if($readID_line)
	{
		die "Weird - read ID should begin with '\@' - line $. of $inputFASTQ - '$readID_line'" unless(substr($readID_line, 0, 1) eq '@');
		my $readID = substr($readID_line, 1);
		$readID =~ s/\s.+//;
		my $read_sequence = <FASTQIN>;
		my $plus = <FASTQIN>;
		my $read_qualities = <FASTQIN>;
		die "FASTQ corruption around line $. of $inputFASTQ - expect a plus character" unless(substr($plus, 0, 1) eq '+');
		die "FASTQ corruption around line $. of $inputFASTQ - expect length(read) == length(qualities)" unless(length($read_sequence) == length($read_qualities));
		if(exists $removeReadIDs{$readID})
		{
			delete($removeReadIDs{$readID});
		}
		else
		{
			$reads_notRemoved++;
			print FASTQOUT $readID_line, "\n", $read_sequence, $plus, $read_qualities;
		}
	}
}
close(FASTQIN);
close(FASTQOUT);
print "Done. $reads_notRemoved remaining reads.\n";
unless(scalar(keys %removeReadIDs) == 0)
{
	die "Could not remove " . scalar(keys %removeReadIDs) . " reads. Abort and remove generated filtered FASTQ.";
	unlink($outputFASTQ);
}

sub print_help
{
	print qq(
removeHumanReads.pl

  Remove reads that map to the human genome.
  
Usage:

  perl removeHumanReads.pl --humanGenomeFASTA myFASTA --inputFASTQ inputFASTQ --outputFASTQ outputFASTQ
  
Example:

 perl removeHumanReads.pl --humanGenomeFASTA humanGenomeFASTA inputFASTQ --outputFASTQ outputFASTQ
  
  
);
exit;
}
