# location: /gpfs/project/dilthey/projects/AMR/analyzeAMR.pl

use strict;
use Data::Dumper;
use Getopt::Long;
use Cwd qw/abs_path getcwd/;
use File::Basename;
use File::Spec qw/abs2rel/;
use lib "/gpfs/project/dilthey/software/MetaMaps_fresh/perlLib";
use taxTree;
use validation;

$| = 1;

my $reference;
my $FASTQ;
my $alignment;
my $identifiedElements;
my $surroundingSequences;
my $threads = 8;

GetOptions (
	'FASTQ:s' => \$FASTQ,
	'reference:s' => \$reference,
	'alignment:s' => \$alignment,
	'identifiedElements:s' => \$identifiedElements,
	'surroundingSequences:s' => \$reference,
	'threads:s' => \$threads
);

die "Please specify --FASTQ" unless($FASTQ);
die "Specified file --FASTQ $FASTQ not existing" unless(-e $FASTQ);


my $fn_output_identifiedElements = abs_path($identifiedElements);
my $fn_output_identifiedElements_withSpecies = abs_path($identifiedElements . '.with_species');
my $fn_output_surroundingSequences = abs_path($surroundingSequences);

my $AMR_href = readFASTA($reference, 1);
my $reads_href = readFASTQ_sequencesOnly($FASTQ, 1);

open(OUTPUT, '>', $fn_output_identifiedElements) or die "Cannot open $fn_output_identifiedElements";
open(SURROUNDINGFASTQ, '>', $fn_output_surroundingSequences) or die "Cannot open $fn_output_surroundingSequences";

my $n_printed_entries = 0;
my %n_alignments_histogram;
my @alignments;
my $currentReadID;
my $processAlignments = sub {
	my $n_alignments = scalar(@alignments);
	$n_alignments_histogram{$n_alignments}++;
	my %read_covered_allAlignmentsCombined;
	my $readLength_outer;
	
	my @primary_alignments = grep {$_->{isPrimary} and not $_->{isSupplementary}} @alignments;
	die "Not exactly one primary aligment for read $primary_alignments[0]{readID}?" unless(scalar(@primary_alignments) == 1);
	my $refID_primary = $primary_alignments[0]->{refID};
	my @alignments_onPrimaryContig = grep {$_->{refID} eq $refID_primary} @alignments;
	unless(scalar(@alignments_onPrimaryContig) == 1)
	{
		warn "Read $primary_alignments[0]{readID} has more than one alignment on $refID_primary - these are ignored.";
	}
	foreach my $alignment (@primary_alignments)
	{
		unless($reads_href->{$alignment->{readID}})
		{
			die "Have no FASTQ entry for $alignment->{readID}";
		}
		my $readSequence = $reads_href->{$alignment->{readID}};
		
		my $readSequence_length = length($readSequence);
		my $alignedReadSequence_asInAlignment = $alignment->{read_alignment};
		(my $alignedReadSequence_asInAlignment_noGaps = $alignedReadSequence_asInAlignment) =~ s/-//g;
		
		my $read_start_smaller = $alignment->{read_start};
		my $read_start_larger = $alignment->{read_stop};
		my $read_covered_bases = $alignment->{read_stop} - $alignment->{read_start} + 1;
		my $isReverseComplement = 0;
		my $read_start_orientationLikeAlignment_0based = $alignment->{read_start} - 1;
		my $read_stop_orientationLikeAlignment_0based  = $alignment->{read_stop} - 1;
		
		die Dumper("Read lengths mismatch", $alignment->{readID}, $alignment->{readLength}, $readSequence_length, $alignment) unless($alignment->{readLength} == $readSequence_length);
		if($read_covered_bases < 0)
		{	
			$read_covered_bases = $alignment->{read_start} - $alignment->{read_stop} + 1;
			$read_start_smaller = $alignment->{read_stop};
			$read_start_larger =  $alignment->{read_start};			
			$isReverseComplement = 1;
			$read_start_orientationLikeAlignment_0based = $readSequence_length - $alignment->{read_start};
			$read_stop_orientationLikeAlignment_0based  = $readSequence_length - $alignment->{read_stop};			
		}
		
		die unless($read_start_orientationLikeAlignment_0based >= 0);
		die unless($read_stop_orientationLikeAlignment_0based >= 0);
		
		# print join("\t", $read_start_orientationLikeAlignment_0based, $read_stop_orientationLikeAlignment_0based), "\n";
		
		my $rawReadSequence_asInAlignment = ($isReverseComplement) ? reverseComplement($readSequence) : $readSequence;		
		my $readSequence_asInAlignment_noGaps_supposed = substr($rawReadSequence_asInAlignment, $read_start_orientationLikeAlignment_0based, $read_stop_orientationLikeAlignment_0based - $read_start_orientationLikeAlignment_0based + 1);

		die Dumper("Read sequence mismatch", length($alignedReadSequence_asInAlignment), length($readSequence_asInAlignment_noGaps_supposed), $alignedReadSequence_asInAlignment_noGaps, $readSequence_asInAlignment_noGaps_supposed) unless($alignedReadSequence_asInAlignment_noGaps eq $readSequence_asInAlignment_noGaps_supposed, $rawReadSequence_asInAlignment);

		die unless($read_covered_bases > 0);

		my $fraction_covered_read = $read_covered_bases / $alignment->{readLength};
		$readLength_outer = $alignment->{readLength};		
		
		my $AMR_ref_length = length($AMR_href->{$alignment->{refID}});
		my $fraction_covered_ref = ($alignment->{ref_stop} - $alignment->{ref_start} + 1) / $AMR_ref_length;
		die unless($fraction_covered_ref > 0);
		
		die "Can't parse ref ID $alignment->{refID}" unless($alignment->{refID} =~ /ARO:\d+\|(\S+)/);
		my $resistanceType = $1;
		
		my $accuracy = sprintf("%.2f", 1 - $alignment->{nm} / $alignment->{alignmentLength});
		
		print $n_alignments, "\t", sprintf("%.2f", $fraction_covered_read * 100) . '% read_coverage [', $read_start_smaller , '-', $read_start_larger, ' of ', $alignment->{readLength}, ']', "  ", sprintf("%.2f", $fraction_covered_ref * 100) . '% ref_coverage [', $alignment->{ref_start} , '-', $alignment->{ref_stop}, ' of ', $AMR_ref_length, ']',, " ", "mapQ=" . $alignment->{mapQ}, " ", "idy=" . $accuracy, " ", $alignment->{refID}, " ", $alignment->{refID}, "\n";
		if(not exists $read_covered_allAlignmentsCombined{$alignment->{refID}})
		{
			$read_covered_allAlignmentsCombined{$alignment->{refID}} = [];
			$#{$read_covered_allAlignmentsCombined{$alignment->{refID}}} = ($alignment->{readLength} - 1);
		}
		die unless($read_start_smaller >= 1);
		die unless($read_start_larger <= $alignment->{readLength});		
		die unless($read_start_smaller < $read_start_larger);
		
		for(my $posI = ($read_start_smaller - 0); $posI <= ($read_start_larger - 1); $posI++)
		{
			$read_covered_allAlignmentsCombined{$alignment->{refID}}[$posI]++;
		}
		
		my $missing_reference_left = $alignment->{ref_start}; die unless($missing_reference_left >= 0);
		my $missing_reference_right = $AMR_ref_length - $alignment->{ref_stop}; die unless($missing_reference_right >= 0);
		my $readSequence_asInAlignment_leftProtruding = substr($rawReadSequence_asInAlignment, 0, $read_start_orientationLikeAlignment_0based);
		my $readSequence_asInAlignment_rightProtruding = substr($rawReadSequence_asInAlignment, $read_stop_orientationLikeAlignment_0based + 1);
		
		my $reconstructed_read_length = (length($alignedReadSequence_asInAlignment_noGaps) + length($readSequence_asInAlignment_leftProtruding) + length($readSequence_asInAlignment_rightProtruding));
		die Dumper("Weird sequence length mismatch", $readSequence_length, $reconstructed_read_length, [length($alignedReadSequence_asInAlignment_noGaps), length($readSequence_asInAlignment_leftProtruding), length($readSequence_asInAlignment_rightProtruding)], [$read_start_orientationLikeAlignment_0based, $read_stop_orientationLikeAlignment_0based], $alignment) unless($readSequence_length == (length($alignedReadSequence_asInAlignment_noGaps) + length($readSequence_asInAlignment_leftProtruding) + length($readSequence_asInAlignment_rightProtruding))); 
		
		my $missing_reference_left_prop = $missing_reference_left / $AMR_ref_length;
		my $missing_reference_right_prop = $missing_reference_right / $AMR_ref_length;
		
		my $verbose = 0;
		if($verbose)
		{
			print "Examining alignment for $alignment->{readID} - along AMR element of length $AMR_ref_length from $alignment->{ref_start} to $alignment->{ref_stop} - along read of length $readSequence_length (aligned) from $read_start_orientationLikeAlignment_0based to $read_stop_orientationLikeAlignment_0based\n";
		}

		my $minAlignmentLength = 500;
		my $max_missingness_bases = 50;
		my $minLength_checkSpecies = 100;

		my $leftOK = 0;
		my $leftSequence_forSpeciesAssignment;
		if($missing_reference_left < $max_missingness_bases)
		{
			print "\tAssume left-hand side of AMR completely covered\n" if($verbose);
			$leftOK = 1;
			if(length($readSequence_asInAlignment_leftProtruding) > $missing_reference_left)
			{
				# we have some protruding sequence, which we correct by the start coordinate in the AMR sequence
				$leftSequence_forSpeciesAssignment = substr($readSequence_asInAlignment_leftProtruding, 0, length($readSequence_asInAlignment_leftProtruding) - $missing_reference_left);
				die unless(length($leftSequence_forSpeciesAssignment) > 0);
			}
			else
			{
				# no protruding sequence
			}
		}
		else
		{
			print "\tLeft-hand side of AMR is not completely covered\n" if($verbose);
			if(length($readSequence_asInAlignment_leftProtruding) < $max_missingness_bases)
			{
				print "\t\t... but the read seems to end here, so OK.\n" if($verbose);
				$leftOK = 1;
			}
			else
			{
				print "\t\t... but the read continues for another " . length($readSequence_asInAlignment_leftProtruding) . " bases, so assume this is not a good alignment.\n" if($verbose);
			}
		}
		
		my $rightOK = 0;
		my $rightSequence_forSpeciesAssignment;
		if($missing_reference_right < $max_missingness_bases)
		{
			print "\tAssume right-hand side of AMR completely covered\n" if($verbose);
			$rightOK = 1;
			if(length($readSequence_asInAlignment_rightProtruding) > $missing_reference_right)
			{
				# we have some protruding sequence, which we correct by the start coordinate in the AMR sequence
				$rightSequence_forSpeciesAssignment = substr($readSequence_asInAlignment_rightProtruding, $missing_reference_right);
				die unless(length($rightSequence_forSpeciesAssignment) > 0);
			}
			else
			{
				# no protruding sequence
			}
		}
		else
		{
			print "\tRight-hand side of AMR is not completely covered\n" if($verbose);
			if(length($readSequence_asInAlignment_rightProtruding) < $max_missingness_bases)
			{
				print "\t\t... but the read seems to end here, so OK.\n" if($verbose);
				$rightOK = 1;
			}
			else
			{
				print "\t\t... but the read continues for another " . length($readSequence_asInAlignment_rightProtruding) . " bases, so assume this is not a good alignment.\n" if($verbose);
			}
		}	

		if(($read_covered_bases >= $minAlignmentLength) and $leftOK and $rightOK)
		{
			# this is a good alignment
			
			$n_printed_entries++;
			
			my @outputFields = ($n_printed_entries, $alignment->{readID}, $resistanceType);
			if($leftSequence_forSpeciesAssignment and (length($leftSequence_forSpeciesAssignment) >= $minLength_checkSpecies))
			{
				push(@outputFields, $leftSequence_forSpeciesAssignment);
				print SURROUNDINGFASTQ '@', $n_printed_entries . '_left', "\n",	$leftSequence_forSpeciesAssignment, "\n", '+', "\n", getPseudoQualities($leftSequence_forSpeciesAssignment), "\n";		
			}
			else
			{
				push(@outputFields, undef);
			}
			
			if($rightSequence_forSpeciesAssignment and (length($rightSequence_forSpeciesAssignment) >= $minLength_checkSpecies))
			{
				push(@outputFields, $rightSequence_forSpeciesAssignment);
				print SURROUNDINGFASTQ '@', $n_printed_entries . '_right', "\n", $rightSequence_forSpeciesAssignment, "\n", '+', "\n", getPseudoQualities($rightSequence_forSpeciesAssignment), "\n";		
			} 
			else
			{
				push(@outputFields, undef);
			}
			
			print OUTPUT join("\t", @outputFields), "\n";
		}
	}
	
	foreach my $refID (sort keys %read_covered_allAlignmentsCombined)
	{
		my $n_bases_covered_combined = scalar(grep {$_> 0} @{$read_covered_allAlignmentsCombined{$refID}});
		my $fraction_covered_combined = $n_bases_covered_combined / $readLength_outer;
		# print "COMBINED ", $n_alignments, "\t", sprintf("%.2f", $fraction_covered_combined * 100) . '%', "\n";
	}
};

open(ALIGNMENTS, '<', $alignments) or die "Cannot open $alignments";
while(<ALIGNMENTS>)
{
	my $meta = $_; chomp($meta);
	my $reference_alignment = <ALIGNMENTS>; chomp($reference_alignment);
	my $read_alignment = <ALIGNMENTS>; chomp($read_alignment);
	<ALIGNMENTS>;
	die "Weird line format in line " . ($. - 2) . " of $alignments:\n$meta" unless($meta =~ /^(\S+) (\S+?):(\d+)-(\d+) read:(\d+)-(\d+) \((\d+)\/(\d+)\) mapQ=(\d+) .+? secondary=(\d+) supplementary=(\d+) readLength=(\d+)/);
	
	
	my %alignment = (
		readID => $1,
		refID => $2,
		ref_start => $3,
		ref_stop => $4,
		read_start => $5,
		read_stop => $6,
		nm => $7,
		alignmentLength => $8,
		mapQ => $9,
		isSecondary => $10,
		isSupplementary => $11,
		readLength => $12,
		reference_alignment => $reference_alignment,
		read_alignment => $read_alignment,
	);
	$alignment{isPrimary} = (! $alignment{isSecondary}) ? 1 : 0;
	
	die "Unknown reference ID $alignment{refID}" unless(exists $AMR_href->{$alignment{refID}});
	
	if($alignment{readID} ne $currentReadID)
	{
		if(@alignments)
		{
			$processAlignments->();
		}
		@alignments = ();
		$currentReadID = $alignment{readID};
	}
	
	push(@alignments, \%alignment);
	
}
close(ALIGNMENTS);

if(@alignments)
{
	$processAlignments->();
}

close(OUTPUT);
close(SURROUNDINGFASTQ);

my $cmd_Kraken = qq(cd /gpfs/project/dilthey/software/MetaMaps_fresh && module load kraken/1.1 && module load Perl && perl callKrakenOnConvertedDB.pl databases/miniSeq+H $fn_output_surroundingSequences $outputDir_kraken noBracken);
system($cmd_Kraken) and die "Could not execute kraken command: $cmd_Kraken";

# warn "Kraken skipped, work with existing clasisfications!";

my $taxonomy = taxTree::readTaxonomy('/gpfs/project/dilthey/software/MetaMaps/databases/miniSeq+H/taxonomy');
my %seqID_2_classification;
my $kraken_output_file = $outputDir_kraken . '/results_kraken.txt.reads2Taxon';
open(CLASSIFICATIONS, '<', $kraken_output_file) or die "Cannot open $kraken_output_file";
while(<CLASSIFICATIONS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	die unless(scalar(@line_fields) == 2);
	$seqID_2_classification{$line_fields[0]} = $line_fields[1];
	
}
close(CLASSIFICATIONS);

open(OUTPUT2, '>', $fn_output_identifiedElements_withSpecies) or die "Cannot open $fn_output_identifiedElements_withSpecies";
open(GENERATEDOUTPUT, '<', $fn_output_identifiedElements) or die "Cannot open $fn_output_identifiedElements";
while(<GENERATEDOUTPUT>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	my $idx = $line_fields[0];
	my $readID = $line_fields[1];
	my $element = $line_fields[2];
	my $seq_left = $line_fields[3];
	my $seq_right = $line_fields[4];
	
	my @output_fields = ($readID, $element);
	my @lca_fields;
	if($seq_left)
	{
		my $seqID = $idx . '_left';
		die unless(exists $seqID_2_classification{$seqID});
		push(@output_fields, $seqID_2_classification{$seqID});
		push(@lca_fields, $seqID_2_classification{$seqID});
	}
	else
	{
		push(@output_fields, 'NA');
	}
	
	if($seq_right)
	{
		my $seqID = $idx . '_right';
		die unless(exists $seqID_2_classification{$seqID});
		push(@output_fields, $seqID_2_classification{$seqID});
		push(@lca_fields, $seqID_2_classification{$seqID});
	}
	else
	{
		push(@output_fields, 'NA');
	}
	
	@lca_fields = grep {$_ != 0} @lca_fields;
	
	if(@lca_fields)
	{
		my $lca = taxTree::lowestCommonAncestor($taxonomy, \@lca_fields);
		push(@output_fields, $lca);
		my $lightning = validation::getAllRanksForTaxon_withUnclassified($taxonomy, $lca, {}, [qw/species genus family order class phylum superkingdom/]);
		push(@output_fields, map {($_ != 0) ? taxTree::taxon_id_get_name($_, $taxonomy) : 'Unclassified'} ($lightning->{species}, $lightning->{genus}, $lightning->{family}, $lightning->{order}, $lightning->{class}, $lightning->{phylum}, $lightning->{superkingdom}));
	}
	else 
	{
		push(@output_fields, 'NA');
		push(@output_fields, 'NA');
		push(@output_fields, 'NA');
		push(@output_fields, 'NA');
	}
	
	push(@output_fields, $seq_left);
	push(@output_fields, $seq_right);

	print OUTPUT2 join("\t", @output_fields), "\n";
}
close(GENERATEDOUTPUT);
close(OUTPUT2);

print "Alignments per read:\n";
foreach my $nAlignments (sort {$a <=> $b} keys %n_alignments_histogram)
{
	print "\t", $nAlignments, " => ", $n_alignments_histogram{$nAlignments}, " reads\n";
}

print "\n\nDone. Generated output in $outputDir\n\n";

foreach my $f (@filesToDelete)
{
	unlink($f) or warn "Could not delete file $f";
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGTacgt/TGCAtgca/;
	return scalar(reverse($kMer));
	return $kMer;
}

sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}

sub readFASTQ_sequencesOnly
{
	my $file = shift;
	my $cut_sequence_ID_after_whitespace = shift;
		
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000000) == 0)
		{
		# 	print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		die unless(substr($line, 0, 1) eq '@');
		my $readID = substr($line, 1);
		
		my $sequence = <F>;
		chomp($sequence);
		$sequence =~ s/[\n\r]//g;
		
		my $plus = <F>;
		chomp($plus);
		$plus =~ s/[\n\r]//g;
		die "No plus character in line '$plus' file $file" unless($plus eq '+');
		
		my $qualities = <F>;
		
		if($cut_sequence_ID_after_whitespace)
		{
			$readID =~ s/\s+.+//;
		}		
		$R{$readID} = uc($sequence);
	}	
	close(F);
		
	return \%R;
}

sub getPseudoQualities
{
	my $FASTA = shift;
	die unless($FASTA);
	my $qualities = ('A' x length($FASTA));
	die unless(length($qualities) == length($FASTA));
	return $qualities;
}
