use strict;
use Data::Dumper;

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Find;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use List::Util qw/shuffle/;
use List::MoreUtils qw/all/;
use Cwd qw/abs_path getcwd/;

my $inputFASTQ;
my $outputDir;
my $database = 'databases/miniSeq+H';
my $maxmemory = '';
my $threads = '';

GetOptions (
	'inputFASTQ:s' => \$inputFASTQ, 
	'outputDir:s' => \$outputDir, 
	'database:s' => \$database,
	'maxmemory:s' => \$maxmemory, 
	'threads:s' => \$threads
);

unless((defined $inputFASTQ) and (defined $outputDir) and (defined $database))
{
	die "Please specify --inputFASTQ --outputDir and --database";
}

unless(-e $inputFASTQ)
{
	die "--inputFASTQ $inputFASTQ not existing";
}

unless(-e $database)
{
	die "--database $database not existing";
}

unless(-d $outputDir)
{
	mkdir($outputDir) or die "Cannot mkdir $outputDir";
}


$inputFASTQ = abs_path($inputFASTQ);
my $metamaps_output_dir = abs_path($outputDir);

$maxmemory = ($maxmemory) ? "--maxmemory $maxmemory " : '';
$threads = ($threads) ? "--threads $threads " : '';

my $metamaps_cmds = qq(./metamaps mapDirectly --all $maxmemory $threads -r ${database}/DB.fa -q $inputFASTQ -o ${metamaps_output_dir}/analysis && ./metamaps classify $threads --mappings ${metamaps_output_dir}/analysis --DB ${database});
system($metamaps_cmds) and die "Command failed: $metamaps_cmds";

