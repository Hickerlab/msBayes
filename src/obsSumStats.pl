#! /usr/bin/perl -w

# obsSumStats.pl
#
# Copyright (C) 2006   Naoki Takebayashi and Michael Hickerson
#
# This program is distributed with msBayes.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

# version: 20060615
# Originally written by ????, and completely rewritten by 
# Mike Hickerson <mhickerson@ccny.cuny.edu> and 
# Naoki Takebayashi <ntakebayashi@alaska.edu>
# Licensing: GPL

## dummy strings (DUM_...) will be replaced.
# simpler version.  The most numbers used in this equation doesn't
# influence the output.  Thre of them (-t, 12.0 and 3.1) influence the
# 3rd column of summary stat vector, but this value is not used in later 
# analysis.

my $headerTmpl = <<HEADER_TEMPLATE;
./msDQH 1 DUM_TotSampleSize 1 -t 4.0 -Q 2.0 0.25 0.25 0.25 0.25 -H 999.0 -r 0 DUM_SeqLen -D 5 2 DUM_SampleSize1 DUM_SampleSize2 0 I 0.0 1.0 1.0 0.9 0.9 3.1 2 1 0 0 1 0 I 0.0 Nc 0.4 0.09 12.0 1 Nc 0.03 DUM_SeqLen 1 Nc 0.03 DUM_TaxonPairID 1 Nc 0.03 DUM_NumTaxonPairs

777
11
//
segsites: DUM_NumSeqSites
freqACGT: 0.25 0.25 0.25 0.25
positions:DUM_positions
HEADER_TEMPLATE
###### END of HERE doc ######

my $usage = "Usage: $0 [-h] [-T tblOutFileName] [-t headerTmplFile] [-s sortPattern] SampleSize_MuModel_Vector\n".
    "  -h: help, print this message\n" .
#    "  -g: Do not remove the sites with any gaps. This probably cause problems.\n".
#    "      So do not use this option." .
    "  -s: Specify the sorting pattern and moments of summary statistics.\n".
    "      By default, the columns are sortedby pi.b (pi between the pairs).\n".
    "      So the column 1 corresponds to the pair with the lowest pi.b\n".
    "      column 2 is the pair with the 2nd lowest\n" .
    "      pi.b etc, and all other summary statistics are ordered accordingly\n" .
    "      Integer value can be used to specify different sorting.\n".
    "       Arguments:\n".
    "        0: do not sort the column at all\n" .
    "        1: simple sort\n" .
    "        2: group by locus then sort by pi.b within each locus(?)\n" .
    "        3: group by locus then sort by the average of pi.b, use 4 moments\n" .
    "        4-7: group by taxon then sort by the average of pi.b\n" .
    "         Moments of summary statistics across loci:\n".
    "          4: use first 4 moments for each summary statistics\n".
    "          5: use first 3 moments for each summary statistics\n".
    "          6: use first 2 moments for each summary statistics\n".
    "          7: use first 1 moment (mean) for each summary statistics (default)\n".
    "      To create the input file for acceptRej.pl, this option\n" .
    "      SHOULD BE LEFT AS THE DEFAULT\n".
    "  -T: Print out an additional human readable table of summary statistics\n" .
    "  -t: an example output file of msDQH (which is used in msbayes.pl)\n" .
    "      You probably do not need to use this option\n" ;

# Note that the aligned sequence data are represented by a simple 1 dimensional
# array.  Each element is a tab-separated (or a character specified by $sep)
# string, which is in the form of "seqName\tseqData".
# It's more efficient to use array of hash or array of array to represent them,
# but it works.

use strict;
use IO::File;
use File::Copy;
use Getopt::Std;
use IPC::Open2;

# Adding the following paths to @INC, so we can find the R scripts.
# The R scripts should be in the same directory as this perl script,
# or inside of ../lib/msbayes/ relative to this perl script.
# i.e. If this script is inside of ~/bin/, the 3 r-scripts should be inside of
# ~/lib/msbayes/.
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../lib/msbayes";
use lib ".";  # looks for required files in the current directory at first

my $sumStatsBinName  = "sumstatsvector";
my $sep = "\t";  # used as the internal separater between seq and seq name
my $filename = "SampleSize_MuModel_Vector";  # default filename

our($opt_h, $opt_t, $opt_g, $opt_T, $opt_s);
my ($verbose);  # assign 1 for verbose (for debug).

getopts('hgt:s:T:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $sumStatsBin = FindFile($sumStatsBinName);  # locate the binary
if ($sumStatsBin eq '-1') {
    die "Can't find $sumStatsBinName in directories:\n", join(":", @INC), "\n";
}

if (@ARGV > 0) {
    if (@ARGV == 1) {
	$filename = shift;
    } else {
	die "Please give only 1 argument\n";
    }
}

my $sumstatsOptions = "";
if (defined($opt_s)) {
    if ($opt_s != 7) {
	warn "INFO: Using the sorting pattern of $opt_s, DO NOT USE THIS as the input to acceptRej.pl\n";
    }
    if ($opt_s < 0 || $opt_s > 7) {
	die "ERROR: argument of -s has to be between 0 and 1\n$usage\n";
    }
    $sumstatsOptions = "-s $opt_s";
} else {
    $sumstatsOptions = "-s 7";
}

my @obsSumStats = CreateObsSumStats($filename, $headerTmpl);

print join "", @obsSumStats;

if (defined($opt_T))
{
    my $dataSummary = "dataSummary.pl";
    my $dsObj = FindExec($dataSummary);
    # `chmod a+x $dsObj`;
    system("$dataSummary $filename > $opt_T");
}

exit;

### This is the actual main function
sub CreateObsSumStats {
    my ($fileName, $header) = @_;

    my @master_matrix = ReadInMaster($fileName);

    my ($numTaxonPairs, $numLoci, $taxonIDHashRef, $locusIDHashRef) = 
	TaxonLocusInfo(@master_matrix);
    my %taxonID = %$taxonIDHashRef;
    my %locusID = %$locusIDHashRef;
    my $numTaxonLocusPairs = @master_matrix;
    # create a matrix from the masterinfile (=$fileName) data.
    # In this master matrix following info will be used later:
    #	v DUM_TotSampleSize = column 0 # total sample
    #	w DUM_SampleSize1 = column 1   # taxon A sample
    #	x DUM_SampleSize2 = column 2   # taxon B sample
    #	y DUM_SeqLen = colomn 4        seqLen
    #	u = column 8                   taxon name
    # column numbers here are 0-offset
    
    # Create an array of filenames in the current working directory
    my @directory_files = glob("*");

    # get the header from a file
    # probably we don't need this, but keep it for now
    if (defined($opt_t)) {
	$header = GetFileContentsAsScalar($opt_t);	
    }
        
    # Each row in @master_matrix corresponds to an alignment file For
    # each row, modify the header template (fake msDQH command line)
    # with appropriate values of @master_matrix, and save it to
    # $sumStatInput.  After the header, extracted variable sites data
    # (and their fake positions) are attached. Basically,
    # $sumStatInput is a long string variable, cotaining fake outputs
    # of msDQH.

    # This info get fed to sumstats program, and receive the results
    # with @sumStatsResultArr.
    my $sumStatInput = "# BEGIN MSBAYES\n" .
	"# numTaxonLocusPairs $numTaxonLocusPairs ".
	"numTaxonPairs $numTaxonPairs numLoci $numLoci\n";
    my $updateSeqLen = 0;
    my @newSeqLen = (); # stores seqLens after removing gaps
    for (my $i = 0; $i < @master_matrix; ++ $i) {
	## Extracting relevant columns
	my ($taxonName, $locusName, $NScaler, $mutScaler) =  @{$master_matrix[$i]}[0..3];
	my ($totSampleSize, $sampleSize1, $sampleSize2) = 
	    @{$master_matrix[$i]}[4..6];
	my ($seqLen, $fastaFile) = @{$master_matrix[$i]}[8,12];
	
	$sumStatInput .= "# taxonID $taxonID{$taxonName} locusID " .
	    "$locusID{$locusName}\n";

	### read in the aligned sequence file, and process it.
	# column 13 (index 12) contains a file name for a taxon-pair
	my $fileName = "";
	if ( -e $fastaFile) {
	    die "ERROR: $fastaFile has 0 size\n" if ( -z $fastaFile);
	    die "ERROR: $fastaFile is not readable\n" if ( ! (-r $fastaFile));
	    $fileName = $fastaFile;
	} else {
	    $fileName = FindSeqFile($fastaFile, \@directory_files);
	}
	if ($fileName eq "") {
	    die "ERROR: Couldn't open a file for taxon, $fastaFile\n";
	}
	
	my @alignedSeq =  (IsFASTA($fileName)) ? ReadInFASTA($fileName) : 
	    ReadInTxt($fileName);
	
	# '?' attached to shorter seqs
	@alignedSeq = AdjustSeqLength(@alignedSeq);
	
	# converting all degenerate characters to '?'
	@alignedSeq = SubstDegenerateCode(\@alignedSeq, '?');

	unless  (defined($opt_g)) {  # sites with any gaps are removed
	    @alignedSeq = RemoveSitesWithGaps(\@alignedSeq);
	}

	# check the seqLen after removal of ambiguous (degenerate) and gaps
	push (@newSeqLen, MaxSeqLen(@alignedSeq));
	if ($seqLen != $newSeqLen[$i]) {
	    $updateSeqLen = 1; # master.batchIn need to be updated
	    $seqLen = $newSeqLen[$i];
	}

	@alignedSeq = ExtractVarSites(\@alignedSeq);
	@alignedSeq = GetSeqDat(@alignedSeq); # get rid of the sequence names
	
	### prepare the header file for this taxon pair
	my $serialNum = $i + 1;
	my $numSeqSites = length($alignedSeq[0]);  # number of variable sites
	
	# printing this to make sure correct files are read in.
	warn "INFO: taxon:locus = $taxonName:$locusName\tfile= $fileName\t".
	    "# variable sites = $numSeqSites\n";

	my $positionString = "";  # making fake positions
	foreach my $srNum (1..$numSeqSites) {
	    $positionString .= sprintf("%7d", $srNum);
	}
	
	my $header_interim = $header;
	$header_interim =~ s/DUM_TotSampleSize/$totSampleSize/g;
	$header_interim =~ s/DUM_SampleSize1/$sampleSize1/g;
	$header_interim =~ s/DUM_SampleSize2/$sampleSize2/g;
	$header_interim =~ s/DUM_SeqLen/$seqLen/g;
	$header_interim =~ s/DUM_NumSeqSites/$numSeqSites/g;
	$header_interim =~ s/DUM_TaxonPairID/$serialNum/g;
	$header_interim =~ s/DUM_NumTaxonPairs/$numTaxonLocusPairs/g;
	# Mike said the next one doesn't matter, but I'm doing it anyway.
	# It is supposed to contain site positions of mutations
	$header_interim =~ s/DUM_positions/$positionString/g;
	
	### some consistency check
	if (@alignedSeq != $totSampleSize) {
	    die "ERROR: taxon, $fastaFile, should have " .
		"$totSampleSize samples.\nBut ", scalar(@alignedSeq),
		" samples are in the file $fileName\n";
	}
	if ($totSampleSize != $sampleSize1+$sampleSize2){
	    die "ERROR: Total sample size ($totSampleSize) for taxon $fastaFile".
		"should be \n       the sum of sampleSizes for the pairs: ".
		"$sampleSize1 + $sampleSize2.\n".
		"Check the sampleSize/mutModel File\n";
	}
	if ($seqLen < $numSeqSites) {
	    die "ERROR: For taxon, $fastaFile, lengths of sequences ($seqLen)".
		" should be\n       longer than number of variable sites " .
		"($numSeqSites)\n";
	}
	
	$sumStatInput .= "$header_interim\n" . join("\n", @alignedSeq) . "\n";
    }  # done with processing all fasta, and making fake msDQH output

    ### run sumstat.
    my $ssPid = open2(\*READ_SS, \*WRITE_SS, "$sumStatsBin $sumstatsOptions -H ");
    print WRITE_SS "$sumStatInput";;
    close(WRITE_SS);  # need to close this to prevent dead-lock
    my @sumStatsResultArr = <READ_SS>;
    close (READ_SS);
    waitpid ($ssPid, 0);

    if (@sumStatsResultArr != 2) {
	die "ERROR: $sumStatsBin returned " . scalar(@sumStatsResultArr) .
	    " lines. It should be 2 lines.\n\nOutput is:\n".
	    join("", @sumStatsResultArr) . "\n";
    }
    
    ### Attaching completely fake prior columns, is this needed? NT Feb 6, 2009
    $sumStatsResultArr[0] = 
	join("\t", ('PRI.Psi', 'PRI.var.t', 'PRI.E.t', 'PRI.omega')) . 
	"\t$sumStatsResultArr[0]";
    $sumStatsResultArr[1] = 
	join("\t", ('1', '0', '1', '0')) . "\t$sumStatsResultArr[1]";
    
    warn "INFO: Total number of (taxon pairs):locus in the data set = $numTaxonLocusPairs\n";

    ### updating the file
    if($updateSeqLen) {
	UpdateSeqLen($fileName, \@master_matrix, \@newSeqLen);
    }
    return @sumStatsResultArr;
}

## this will update the config file $fileName with the corrected sequence length
sub UpdateSeqLen {
    my ($fileName, $master_matrixRef, $newSeqLenRef) = @_;
    my @master_matrix = @$master_matrixRef;
    my @newSeqLen = @$newSeqLenRef;
    my $newSampleTbl = "";
    print STDERR 
	"INFO: seqLen of $fileName is adjusted (to the length after\n".
	"INFO: removals of sites with any gaps or ambiguous/degenerate\n".
	"INFO: characters). The original file is saved as $fileName.orig\n".
	"INFO: The fasta files are NOT modified at all.\n";
    
    for my $i (0..$#master_matrix) {
	my @thisRow = @{$master_matrix[$i]};
	if ($thisRow[8] != $newSeqLen[$i]) {
	    # update 9-th column (seqLen).
	    my $nRmChar = $thisRow[8] - $newSeqLen[$i];
	    my $action = "removed";
	    if ($nRmChar < 0) {
		$action = "added";
		$nRmChar *= -1; # convert to positive number
	    }
	    print STDERR "INFO: $thisRow[0]:$thisRow[1], $nRmChar sites ".
		"$action\n";
	    $thisRow[8] = $newSeqLen[$i];
	}
	
	# remove 5-th column totalSampleSize
	splice(@thisRow, 4, 1);
	$newSampleTbl .= join("\t", @thisRow) . "\n";
    }
    
    CheckNBackupFile("$fileName.orig", 'file');	
    move($fileName, "$fileName.orig") || 
	die "Can't move $fileName to $fileName.orig";
    open NEWCONF, ">$fileName" || die "Can't reopen $fileName\n";
    open OLDCONF, "<$fileName.orig" || die "Can't open $fileName.orig\n";
    my $state = 0;
    while(<OLDCONF>) {
	my $saveLine = $_;
	s/#.*$//;
	chomp;
	if (/BEGIN\s+SAMPLE_TBL/) {
	    if ($state == 0) {
		print NEWCONF "BEGIN SAMPLE_TBL\n$newSampleTbl\n".
		    "END SAMPLE_TBL\n";
		$state++;
		next;
	    } else {
		warn "WARN: there are more than 1 SAMPLE_TBL section.\n" .
		    "WARN: Check $fileName to make sure it is correct.\n";
	    }
	}
	
	if ($state == 1) { # inside of SAMPLE_TBL
	    if (/END\s+SAMPLE_TBL/) {
		$state++;
		next;
	    } else {
		next;
	    }
	}
	# when reached here, it is outside of SAMPLE_TBL
	print NEWCONF "$saveLine";
    }
    
    close OLDCONF;
    close NEWCONF;

    return;
}

# subroutine to get header template in a file and put it into a scalar
# Takes 1 argument: the file name
sub GetFileContentsAsScalar {    
    my $filename = shift;
    my $contents = "";  # scalar with data
    
    open(INPUT, $filename) || die "Cannot open file \"$filename\"!\n";
    while(<INPUT>) {
	$contents .= $_;
    }
    close INPUT;
    return $contents;
}

sub TaxonLocusInfo {
    my @matrix = @_;
    my %nameHash = ();
    my %locusNameHash = ();
    my $cntr = 1;
    for my $i (0..$#matrix) {
	next if (defined ($nameHash{$matrix[$i][0]}));
	$nameHash{$matrix[$i][0]} = $cntr;
	$cntr++;
    }
    my $numUniqTaxon = $cntr - 1;
    $cntr = 1;
    for my $i (0..$#matrix) {
	next if (defined ($locusNameHash{$matrix[$i][1]}));
	$locusNameHash{$matrix[$i][1]} = $cntr;
	$cntr++;
    }
    my $numUniqLoci = $cntr - 1;

    return ($numUniqTaxon, $numUniqLoci, \%nameHash, \%locusNameHash);
}

# Takes a filename as an argument
# Parse the file, and returun a 2-dim matrix
# The file contains information about sample sizes and mutational models.
# It should be tab delimited text file with following columns.
#   1: TaxonPairName
#   2: locusName
#   3: NScaler (0.25 for chloroploast in dioecy, 1 for nuclear genes, etc.)
#   4: mutScaler (e.g. 10.0, mtDNA has a higher mutation rate)
#   5: TotalSampleSize
#   6: SampleSize1 
#   7: SampleSize2
#   8: transition/transversion Ratio
#   9: baseTotalpairs (length of sequences)
#  10: Afreq 
#  11: Cfreq
#  12: Gfreq
#  13: Filename
# Each line contains a data for 1 taxon pair (sp.1 and sp.2).
# The returned 2-dim matrix contain these information, each line = each row.
# '#' is used to indicate comments, and ignored.
# Empty lines are ignored.
# The configuration info for msprior can be included before this data matrix.
# They include '=', so the first non-empty line which does not contain the
# '=' is considered as the beginning of the samplesize/mut. model data.
sub ReadInMaster {
    my $filename = shift;
    open(INPUT, $filename) || die "Cannot open file \"$filename\"!\n";
    my @master_matrix;
    my $numCol = -1;
    my $paraConfSect = 1;
    my $warnFlag = 1;
    while (<INPUT>) {
	s/#.*$//;  # remove any thing after "#", which indicates comments
	next if (/^\s*$/);
	if ($warnFlag) {   # Maybe we should handle line ending better at some point
	    if (/\r\n$/) {
		warn "WARN: the file $filename appears to have DOS style line ending, it may not work correctly if it's not converted to unix style newline characters\n";
	    }
	    $warnFlag = 0;
	}

	# There could be parameter config lines for msprior in the
	# beginning of this file, which should be ignored.  These are
	# in the format of "parameterName = value", so they are
	# indicated by existence of '=' character. The first line
	# which doesn't contain '=' is considered to be the 1st line
	# with mutation parameter data.
	if ($paraConfSect) {
	    if (/BEGIN\s+SAMPLE_TBL/i) {
		$paraConfSect = 0;
		next;
	    } elsif (/=/) {
		next;
	    } else { # this is for old format without the constraints
		$paraConfSect = 0;
		# continue with analysis
	    }
	}
	
	# when reached here, we are getting into sample sizes/mutational para
	if (/END\s+SAMPLE_TBL/i) {
	    last;
        }

	my @master = split /\t/; 
	chomp ($master[$#master]); # remove newline in the taxon name cell
	# check all rows have good column numbers
	if ($numCol < 0) {
	    $numCol = @master;
	    if ($numCol >13 || $numCol < 12) {
		die "ERROR: reading $filename, the 1st line of sample " .
		    "sizes/mutation parameter lines should have 12 or 13 " .
		    "columns.  But, it has $numCol:  $_\n";
	    }
	} elsif ($numCol != @master) {
	    die "ERROR: the number of columns in $filename should be $numCol, " .
		"but the following line has ", scalar(@master), " columns: $_\n";
	}
	
	# total sample numbers should be calculated from the sample sizes of
        # taxon pairs.
	if ($numCol == 12) {
	    splice @master, 4, 0, $master[4] + $master[5];
	}
	push @master_matrix, [ @master ];
    }
    close INPUT;
    
    return @master_matrix;
}

# This have a problem (original method)
# E.g., if "taxa1.txt" and "taxa12.txt" exist "taxa1" matches with both file
sub FindSeqFileVague {
    my ($name, $arrRef) = @_;
    
    my @matching = grep { $_ =~ /$name/ } @$arrRef;

    if (@matching == 1) {
	return $matching[0];
    } else {
	return "";
    }
}

# Takes the name for the taxon pair as the first arg. and array ref to
# an array containing the name of files in the directory as the 2nd arg.
# It returns the candidate file for the taxon pair.
# If it can't identify a single file, it returns "".
#
# I'm assuming that the filename is taxa1.txt, taxa12.fasta, taxa3,
# etc.  It takes the part of filename before the 1st '.' and check if
# it matches with the name for the taxon pair.  Also exact match is
# ok.  case-insensitive match
sub FindSeqFile {
    my ($name, $arrRef) = @_;
    
    my @matching = grep { ($_ =~ /^$name\./ || $_ =~ /^$name$/) 
			      && $_ !~ /~$/ } @$arrRef;

    if (@matching == 1) {
	return $matching[0];
    } elsif (@matching > 1) {
	warn "ERROR: trying to find the sequence file for taxon pair, $name.\n".
	    "       Multiple matches: ", join(" ", @matching), "\n";
	return "";
    }

    # no match, so ignore the case and see if some file matches
    @matching = grep { $_ =~ /^$name\./i || $_ =~ /^$name$/i} @$arrRef;
    if (@matching == 1) {
	if ($verbose) {
	    warn "Trying to find filename containing $name, similar to " .
		"$matching[0].  This is used\n";
	}
	return $matching[0];
    } elsif (@matching > 1) {
	warn "ERROR: trying to find the sequence file for taxon pair, $name.\n".
	    "       Multiple matches: ", join(" ", @matching), "\n";
	return "";
    }
    # we could use the original relaxed match here, but I'm not doing it
    return "";
}

# each line contains one sequence (sample), no names are give.
sub ReadInTxt {
    my $infile = shift;
    open(INFILE, "<$infile")  || die "Can't open $infile\n";
    
    my @result = ();
    my $baseName = "seq";
    my $cntr = 1;
    while(<INFILE>) {
	chomp;
	next if (/^\s*$/);
	push @result, "$baseName$cntr\t$_";
	$cntr++;
    }
    return @result;
}

# Take a filename as the arg.
# Check if it appear to be fasta file or not.
# Returns:
#   1: the file is FASTA
#   0: the file is not FASTA
# criterion:
#  First non empty line starts with \s*>someName  &&
#  The non empty line after '>' does not conatin '>'
sub IsFASTA {
    my $name = shift;
    open(INFILE, "<$name") || die "Can't open $name\n";
    my $state = 0;
    while(<INFILE>) {
	chomp;
	next if (/^\s*$/);
	
	if ($state == 0) {
	    if (/^\s*>.+$/) {
		$state = 1;
	    } else {
		close INFILE;
		return 0;
	    }
	} elsif ($state == 1) {
	    if (/^\s*>/) {
		close INFILE;
		return 0;
	    } else {
		close INFILE;
		return 1;
	    }
	}
    }
    close INFILE;
    return 0;
}

# takes an arg; name of a file from which data are read Then read in
# the data and make an array.  Each element of this array corresponds
# to a sequence, name tab data.
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";
    
    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
            s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            tr/uU/tT/;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . $_;
        }

        # checking no occurence of internal separator $sep.
        die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
             "the input FASTA file contains this charcter. Make sure this " .
             "separator character is not used in your data file or modify " .
             "variable \$sep in this script to some other character.\n")
            if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
        $result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}

# take std seq data array (name\tseq), and attach "?" at the end for
# shorter sequences
sub AdjustSeqLength {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);
    my $maxLen = MaxSeqLen(@_);

    foreach my $i (0 .. $#seqDat) {
        my $thisLen = length($seqDat[$i]);
        if ($thisLen == $maxLen)  {
            ; # do nothing
        } elsif ($thisLen < $maxLen) {
            my $diff = $maxLen - $thisLen;
            warn "WARN: $seqName[$i] shorter.  " .
                "$diff '?' (missing character) were added at the end\n";
            for (my $j=0; $j < $diff; $j++) {
                $data[$i] = $data[$i] . "?";
            }
        } else {
            die "ERROR: the length of sequence $seqName[$i] is $thisLen, " .
                "longer than \$maxLen = $maxLen.  Weird!!";
        }
    }
    return (@data);
}

# convert any degenerate code to the character specified by $substChar
sub SubstDegenerateCode {
    my ($datArrRef, $substChar) = @_;
    my @seqDat = GetSeqDat(@$datArrRef);
    my @seqName = GetSeqName(@$datArrRef);
    my $maxLen = MaxSeqLen(@$datArrRef);

    my @result = ();
    foreach my $i (0..$#seqDat) {
	$seqDat[$i] =~ s/[RYKMSWBDHVN]/$substChar/gi;
	push @result, "$seqName[$i]\t$seqDat[$i]";
    }
    return @result;
}

# '-' or '?' are considered as gaps
sub RemoveSitesWithGaps {
    my $datArrRef = shift;
    my @seqDat = GetSeqDat(@$datArrRef);
    my @seqName = GetSeqName(@$datArrRef);
    my $maxLen = MaxSeqLen(@$datArrRef);
    my @gapSites = ();
    my @notGapSites = ();
    my ($posi, $seqNumber);
    my @seqMat = ();

    # make 2 dimensional matrix
    foreach $seqNumber (0..$#seqDat) {
        my @tmpArray = split(//, $seqDat[$seqNumber]);
        # Check the length
        if (@tmpArray != $maxLen)  {
            die "ERROR: the sequence $seqName[$seqNumber] is not same length ".
                "as \$maxLen = $maxLen.  Weird!!";
        }
        push @seqMat, [ @tmpArray ];
    }
    
    # now identify the sites with any gap.
    for $posi (0 .. ($maxLen-1)) {
        my $gap = 0;
        for $seqNumber (0 .. $#seqMat){
            if ($seqMat[$seqNumber][$posi] =~ /^[-\?]$/) {
                $gap = 1;
                last;
            }
        }
        if ($gap == 1) {  #  a gap at these sites
            push (@gapSites, $posi+1); # now unit-offset
        } else {          # there are some non-gap character at these sites
            push (@notGapSites, $posi);
        }
    }
    
    # select sites and make results
    my @result = ();
    for $seqNumber (0 .. $#seqMat) {
        my @thisSeq = SelectSites($seqMat[$seqNumber], \@notGapSites);
        my $line = $seqName[$seqNumber] . $sep . (join("", @thisSeq));
        push (@result, $line);
    }
    
    if ($verbose && @gapSites > 0) {
        warn ("Following sites contains gap(s), removed from analysis\n");
        print STDERR join(" ", @gapSites);
        print STDERR "\n";
    }
    return (@result);
}

## select sites with any variation.  This include varaition due to gaps or 
# degenerate (ambiguous) sites
sub ExtractVarSites {
    my $datArrRef = shift;
    my @seqDat = GetSeqDat(@$datArrRef);
    my @seqName = GetSeqName(@$datArrRef);
    my $maxLen = MaxSeqLen(@$datArrRef);

    my ($posi, $seqNumber);
    my @seqMat = ();

    # make 2 dimensional matrix
    foreach $seqNumber (0..$#seqDat) {
        my @tmpArray = split(//, $seqDat[$seqNumber]);
        # Check the length
        if (@tmpArray != $maxLen)  {
            die "ERROR: the sequence $seqName[$seqNumber] is not same length ".
                "as \$maxLen = $maxLen.  Weird!!";
        }
        push @seqMat, [ @tmpArray ];
    }
    
    # now identify the sites with variable sites
    my @varSites = ();
    for $posi (0 .. ($maxLen-1)) {
        my $variable = 0;
	my $char1stSeq = $seqMat[0][$posi];
        for $seqNumber (1 .. $#seqMat){
            if ($seqMat[$seqNumber][$posi] !~ /$char1stSeq/i) {
                $variable = 1;
                last;
            }
        }
        if ($variable == 1) {  #  variable site
            push (@varSites, $posi);
        }
    }
    
    # select sites and make results
    my @result = ();
    for $seqNumber (0 .. $#seqMat) {
        my @thisSeq = SelectSites($seqMat[$seqNumber], \@varSites);
        my $line = $seqName[$seqNumber] . $sep . (join("", @thisSeq));
        push (@result, $line);
    }
    
    return (@result);
}

sub SelectSites {
    my ($arrayRef, $indexRef) = @_;
    unless (@_ == 2 && ref($arrayRef) eq 'ARRAY' && ref($indexRef) eq 'ARRAY'){
        die "args to SelectSites() should be ARRAY REF, ARRAY REF\n";
    }

    my $maxIndex = @$arrayRef -1;
    my @result = ();
    foreach my $posi (@$indexRef) {
        if ($maxIndex < $posi) {
            push @result, "?";
        } else {
            push @result, $$arrayRef[$posi];
        }
    }
    return @result;
}

sub MaxSeqLen {
    my @data = GetSeqDat(@_);
    my $maxLen = 0;
    foreach my $i (@data) {
        my $len = length($i);
        $maxLen = $len if ($len > $maxLen);
    }
    return ($maxLen);
}

sub GetSeqDat {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
        @line = split (/$sep/, $i);
	# Note that @line == 1 when the sequence data is empty: ""
	# Instead of returning undef, it's returning "" in this case.
	push @result, (@line == 1) ? "" : $line[1];
    }

    return (@result)
}

sub GetSeqName {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
        @line = split (/$sep/, $i);
        push @result, $line[0];
    }
    return (@result)
}


# Check @INC library locations, and find a file
sub FindFile {
    my $name = shift;

    foreach my $dir (@INC) {
        if ( -e "$dir/$name" ) {
	    print STDERR "FILEINFO: using $dir/$name\n";
            return "$dir/$name";
        }
    }
    return -1;
}

# This fucntion check if the argument (fileName) exists. If it exists,
# it get renamed to fileName.oldN, where N is a digit.
# In this way, no files will be overwritten.
# 2nd argument $type is either 'file' or 'dir'.
# It will create file or directory with the name $fileName.
# rmtree() requires File::Path
use File::Path;
sub CheckNBackupFile {
    my ($fileName, $type) = @_;
    my $maxSave = 3;

    if (-e $fileName) {
	my $i = 1;
	while (-e "$fileName.old$i") {  # checking if the file exists
	    $i++;
	}
	
	if ($i > $maxSave) {
	    $i = $maxSave;
	    rmtree("$fileName.old$i") || die "Can't delete $fileName.old$i";
	}
	
	# oldest file has .old5, newest file has .old1
	while ($i > 1) {
	    move("$fileName.old" . ($i - 1), "$fileName.old$i") ||
		die "Can't rename $fileName.old" . ($i-1) . 
		" to $fileName.old$i";
	    $i--;
	}
	
	move("$fileName", "$fileName.old$i") ||
	    die "Can't rename $fileName to $fileName.old$i";
    }
    if ($type eq 'file') {
	# create the empty outfile, so other processes don't use the name.
	open(OUT,">$fileName");
	close(OUT);
    } elsif ($type eq 'dir') {
	mkdir $fileName;
    }
}

### Supply the name of program, and it will try to find the executable.
### In addition to regular path, it search for several other places
sub FindExec {
    my $prog = shift;
    # I'm making it to find the binary in the current directory (.) at first.
    # I do not personally like this.  But since Mike doesn't like to
    # install the binaries in the appropriate directories, we need to
    # force this behavior to reduce confusion. 
    # When this program become more matured, we should reevaluate this.
    # Similar behavior in acceptRej.pl introduced  Naoki Feb 8, 2008
    $ENV{'PATH'} = ".:" . $ENV{'PATH'} . 
	":/bin:/usr/bin:/usr/local/bin:$ENV{'HOME'}/bin";
    my $bin = `which $prog 2>/dev/null`;
    chomp $bin;

    if ($bin eq "") {
	die "ERROR: $prog not found in PATH $ENV{'PATH'}\n";
    }

    print STDERR "INFO: using $bin\n";
    return $bin;
}

