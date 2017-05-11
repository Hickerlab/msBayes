#!/usr/bin/perl -w

# acceptRej.pl
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


my $pdfOut = "figs.pdf";  # default pdf output file name
my $defaultAcceptCnt = 500;  # used to set default tolerance
my $statString = "'pi','wattTheta','pi.net','tajD.denom'";  # default stat strings
my $defaultAcceptedFile = "posterior_table";

my $usage="Usage: $0 [-hncr] [-p outPDF] [-a acceptedOutFileName] [-s summary_stats]\n" .
          "                       [-t tolerance] obsData simData \n".
    "  -h: help\n".
    "  -n: print out the names of all available summary stats\n".
    "  -r: Simple rejection method without regression will be used\n".
    "  -p: output pdf filename (default: $pdfOut)\n".
    "  -a: specify a filename of accepted parameter values and summary stats\n".
    "      Without this flag, a default filename \"$defaultAcceptedFile\" is used.\n".
    "  -s: statString (e.g. -s '$statString', <=default)\n".
    "      The summary statistics listed here will be used\n".
    "  -t: tolerance (a value between 0 an 1, default: set the tolerance so that\n".
    "      $defaultAcceptCnt simulations will be accepted)\n".
    "  -i: Leave intermediate data files used for R and print out the\n".
    "      filenames of these temporary files.  By default, these files get\n".
    "      erased after the analysis.  These files can be useful for\n".
    "      debugging.  Or these files can be used to conduct more detailed\n".
    "      analysis or to create custom figures in R.  With the manual R\n".
    "      analysis, you can source() the file pointed by \$tmpR to set up\n".
    "      the environment.  A list object \"res\" contains relevant data.\n" 
#    "  -a: old analysis with all R processing (nobody needs it)\n"
    ;

use Getopt::Std;
getopts('a:hindp:rt:s:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

our($opt_a, $opt_h, $opt_i, $opt_n, $opt_d, $opt_p, $opt_r, $opt_t, $opt_s);

use File::Copy;
use IO::File;
use POSIX qw(tmpnam);

# if you set this to 1, it will not use msReject, and use old way of
# doing analysis with just R.  But R only analysis is a memory hog,
# and this option is not needed any more.
my $ROnlyAnalysis = 0;

# R command
my $Rbin = "R";

# The names (not paths) of 3 R scripts and a C prog required for this program.
my $mainRscript = "acceptRej.r";
my $make_pdRscript = "make_pd2005.r";
my $loc2plotRscript = "loc2plot.r";
my $calmodRscript = "calmod.r";
my $rejectionExe = "msReject";         # rejection program
my $sumStatExe = "sumstatsvector";

# Adding the following paths to @INC, so we can find the R scripts.
# The R scripts should be in the same directory as this perl script,
# or inside of ../lib/msbayes/ relative to this perl script.
# i.e. If this script is inside of ~/bin/, the 3 r-scripts should be inside of
# ~/lib/msbayes/.
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../lib/msbayes";
use lib ".";  # looks for Required files in the current directory at first
              # -d will print out the search path
if (defined($opt_d)) {
    print STDERR "FILEINFO: searching path is ";
    print STDERR join(":", @INC);
    print STDERR "\n";
}

if(defined($opt_n)) {
    PrintSumStatNames();
    exit(0);
}

###  create a bunch of temp files
my $rmTempFiles = 1;  # set this to 0 for debugging
if (defined($opt_i)) {
    $rmTempFiles = 0;
}

# open a temporary file to store the dynamically created R script
do {$tmpR = tmpnam()} until $tmpRfh = 
    IO::File->new($tmpR, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if ($rmTempFiles) {
	if (defined($tmpR) && -e $tmpR) {
	    unlink($tmpR) || die "Couldn't unlink $tmpR : $!"
	}
    } else {
	print STDERR "FILE: \$tmpR = $tmpR\n" if(defined($tmpR));
	
    }
};

# open a temp file to preprocess the observed data
do {$tmpObs = tmpnam()} until $tmpObsfh = 
    IO::File->new($tmpObs, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if ($rmTempFiles) {
	if (defined($tmpObs) && -e $tmpObs) {
	    unlink($tmpObs) || die "Couldn't unlink $tmpObs : $!";
	}
    } else {
	print STDERR "FILE: \$tmpObs = $tmpObs\n" if (defined($tmpObs));
    }
};
$tmpObsfh->close();

# open a temp file to preprocess the data with msrejection
my $tmpSimDatfh;
do {$tmpSimDat = tmpnam()} until $tmpSimDatfh = 
    IO::File->new($tmpSimDat, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if ($rmTempFiles) {
	if (defined($tmpSimDat) && -e $tmpSimDat) {
	    unlink($tmpSimDat) || die "Couldn't unlink $tmpSimDat : $!";
	}
    } else {
	print STDERR "FILE: \$tmpSimDat = $tmpSimDat\n" if (defined($tmpSimDat));
    }
};

# open a temp file to remove the header from sim dat
my $tmpSimDat2fh;
do {$tmpSimDat2 = tmpnam()} until $tmpSimDat2fh = 
    IO::File->new($tmpSimDat2, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if ($rmTempFiles) {
	if (defined($tmpSimDat2) && -e $tmpSimDat2) {
	    unlink($tmpSimDat2) || die "Couldn't unlink $tmpSimDat2 : $!";
	}
    } else {
	print STDERR "FILE: \$tmpSimDat2 = $tmpSimDat2\n" if (defined($tmpSimDat2));
    }
};

# open a temp file to extract the prior columns.
do {$tmpPrior = tmpnam()} until $tmpPriorfh = 
    IO::File->new($tmpPrior, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if ($rmTempFiles) {
	if (defined($tmpPrior) && -e $tmpPrior) {
	    unlink($tmpPrior) || die "Couldn't unlink $tmpPrior : $!";
	}
    } else {
	print STDERR "FILE: \$tmpPrior = $tmpPrior\n" if (defined($tmpPrior));
    }
};

#### end of setting-up temp files

if (defined($opt_s)) {
    $statString = MkStatString($opt_s);
}

if (@ARGV != 2) {
    warn "ERROR: This script requires two arguments";
    die $usage;
}

my ($obsDat, $simDat) = @ARGV;

if (! (-r $simDat && -T $simDat )) {
    warn "ERROR: Problem with $simDat. Please give a readable, non-empty ".
	"file name\n";
    die $usage;
}
if (! (-r $obsDat && -T $obsDat )) {
    warn "ERROR: Problem with $obsDat. Please give a readable, non-empty ".
	"file name\n";
    die $usage;
}

if(defined($opt_p)) {
    $pdfOut=$opt_p;
}
CheckNBackupFile($pdfOut, 'file');

## save the accepted sims (+ transfomed prior values + sum stats) to a file
my $acceptedFile = (defined ($opt_a)) ? $opt_a : $defaultAcceptedFile;
CheckNBackupFile($acceptedFile, 'file');
copy($tmpSimDat, $acceptedFile) || 
    warn "WARN: copying $tmpSimDat to $acceptedFile failed: $!";

{   # Making sure R is installed
    # Mac may have R64 for 64 bit operation
    my $checkR64 = system("which R64 > /dev/null");
    if ($checkR64 == 0) {
	$Rbin = "R64";
    } else {
	my $checkR = system("which R > /dev/null");
	die "\nERROR: Can't find R.  Make sure R is installed and in your PATH.\n" 
	    unless ($checkR == 0);
    }
}

## preprocess the simDat

## find the number of columns
my $numColInfile = ColNumTabDelimFile($simDat);

# getting column names/structures
my ($arrRef1, $arrRef2, $numTaxonPairs) = GetPriorSumStatNames($simDat);
my @priorNames = @$arrRef1;
my @sumStatNames = @$arrRef2;
my $numPriorCols = scalar(@priorNames);

# processing the header of obsData
my $numTaxonPairs2;
($arrRef1, $arrRef2, $numTaxonPairs2) = GetPriorSumStatNames($obsDat);
if (@sumStatNames != scalar (@$arrRef2)) {
    my $sn1 = scalar (@sumStatNames);
    my $sn2 = scalar (@$arrRef2);
    die "ERROR: mismatch in summary stats of obsData ($sn2 columns) ".
	"and simData ($sn1 columns)\n";
} elsif ($numTaxonPairs != $numTaxonPairs2) {
    die "ERROR: obsData contains $numTaxonPairs2 taxon pairs and simData".
	" contains $numTaxonPairs taxon pairs.\n";
}

# making sure summary stats headers are matching bet simDat and obsDat
for my $hCnt (0..$#sumStatNames) {
    if($sumStatNames[$hCnt] ne $$arrRef2[$hCnt]) {
	die "ERROR: simData has $sumStatNames[$hCnt] but obsData ".
	    "has $$arrRef2[$hCnt]\n";
    }
}

# read in obsData
open OBS, "<$obsDat" || die "Can't open $obsDat\n";
my @obsDataArray = <OBS>;
close OBS;
if (@obsDataArray != 2) {
    die "ERROR: Observed data should have two lines: 1 header line and ".
	"another data line\n";
}

# adjust the column length of obsData
my @tmpObsVectLine = split /\t/, $obsDataArray[1];
my @tmpObsHeader = split /\t/, $obsDataArray[0];
my %obsHash = ();
foreach my $index (0..$#tmpObsHeader) {
    if (defined($obsHash{$tmpObsHeader[$index]})) {
	die "ERROR: Processing observed summary stat data, and the  ".
	    "header doesn't contain unique column names: $tmpObsHeader[$index]";
    }
    $obsHash{$tmpObsHeader[$index]} = $tmpObsVectLine[$index];
}


# removing the prior columns in obsData, and replacing it with 0 with the
# length matching with the simulated data set.
#splice @tmpObsVectLine, 0, scalar(@$arrRef1), Rep(0, scalar(@priorNames));

# Making sure the columns in obsData are consistent with simulated data set.
# When fake observed data are created with different model, PRI.* columns
# may not match.  If obsData are missing columns, dummy 0 are inserted.
# This should not influence the subsequent calculation
@tmpObsVectLine = ();
foreach my $cname (@priorNames, @sumStatNames) {
    if (defined($obsHash{$cname})) {
	push @tmpObsVectLine, $obsHash{$cname};
    } else {
	if ($cname !~ /^PRI\./) {  # it should never come here
	    die "ERROR: summary statistics for observerd and simulated data ".
		"do not match\n";
	}
	push @tmpObsVectLine, 0; # pushing in dummy 0
    }
}
$obsDataArray[1] = join "\t", @tmpObsVectLine;
$obsDataArray[0] = join "\t", (@priorNames, @sumStatNames);
$obsDataArray[0] .= "\n";

## finding the length of SIMDAT
my $numSimCount = `grep -v PRI.Psi $simDat | wc -l`;
die "ERROR: grep -v PRI.Psi $simDat | wc -l failed: $?" if $?;
chomp $numSimCount;

## Set tolerance
my $tol;
if (defined($opt_t)) {
    $tol = $opt_t;
} else {
    my $defaultTolerance = $defaultAcceptCnt / $numSimCount;
    $tol = (defined($opt_t)) ?  $opt_t :  $defaultTolerance;

    if ($defaultTolerance > 1) {
	$defaultTolerance = 1;
	print STDERR 
	    "WARN: Hmmm, you ran only $numSimCount simulations.  You'll need about".
	    "500 accepted simulations for reasonable estimation. Tolerance is set".
	    "to 1, but this analysis is probably not meaningfull.\n";
    }
    print STDERR "INFO: tolerance set to $defaultTolerance.\n" ;
}

#### Making the R script
MkStdAnalysisRScript($tmpRfh);
close($tmpRfh);  # the tmp R script is ready to use

if(defined($opt_d)) {
    open RSCRIPT, "<$tmpR" || die ("ERROR: Can't open $tmpR\n");
    warn "### debug: contents of R script file\n";
    while(<RSCRIPT>) {
	print STDERR $_;
    }
    close RSCRIPT;
    warn "### debug: end of R script file\n";
}

if (! $ROnlyAnalysis) {  # use the external acceptRejection C program
    
    ## create the prior columns only file, and a file without header
    open SIMDAT, "<$simDat" || die "Can't open $simDat\n";

    ## Setting up thinning of prior for Bayes Factor.  Sampling
    ## approximately 100k prior to make its distn
    my $priorSampleIntvl = ($numSimCount>100000) ? int($numSimCount/100000) : 1;
    my $savedHeader = "";
    while (<SIMDAT>) {
	if ($. != 1) {  # creating a file without the header for rejection
	    if (/PRI\.numTauClass/) {  # multiple concatenated sims
		if ($_ eq $savedHeader) { # ignore the line
		    next;
		} else {
		    die "ERROR: You concatenated the multiple simulation ".
			"results.  But the headers do not match (indicating ".
			"that the simulations weren't run with the same ".
			"parameters). 1st header was\n\n$savedHeader\n".
			"The header at line $. was\n\n$_\n";
		}
	    }
	    print $tmpSimDat2fh "$_";
	} else {  # 1st line, prepare the header for the output of rejection 
	    print $tmpSimDatfh "$_";
	    $tmpSimDatfh->close();
	    $savedHeader= $_;  # used to check multiple headers when
	                       # multipe simulations are concatenated
	}
	
	# shrink the prior for R
	if (($. -  1) % $priorSampleIntvl == 0) {
	    chomp;
	    my @line = split /\t/;
	    my @priors = splice @line, 0, $numPriorCols;
	    print $tmpPriorfh join("\t", @priors), "\n";
	}
    }
    close SIMDAT;
    $tmpSimDat2fh->close();
    
    # finding the column numbers to use as the summary statistics
    my @usedSS = split /\s*,\s*/, $statString;
    print STDERR "INFO: Using following summary statistics: ";
    print STDERR join(",", @usedSS), "\n";
    my @index=();
    foreach my $ss (@usedSS) {
	my @tmp = FindMatchingIndex($ss, @sumStatNames); # 0-offset
	if (@tmp > $numTaxonPairs) {
	    die "ERROR: More than 1 summary stat matches with $ss.  Please ".
		"inform this bug to the developper\n";
	} elsif (@tmp != $numTaxonPairs) {
	    die "ERROR: Could not find specified summary stats $ss\n";
	}
	my @tmp2 = map { $_ + 1 + scalar(@priorNames)} @tmp;
	push @index, @tmp2;  # @index is 1-offset
    }

    # run rejection program
    my $columns = join " ", @index;
    my $rejExe = FindFile($rejectionExe);
    if ($rejExe eq '-1') {
	die  "ERROR: Can't find $rejectionExe, is it installed in your PATH?\n";
    }
    
    # Making a temporary obsDat for rejection
    open OBS, ">$tmpObs" || die "Can't open temporary obs file $tmpObs\n";
#    for (my $extra = 0; $extra < $numExtraPrior; $extra++) { # padding
#	print OBS "0\t";
#    }
    print OBS $obsDataArray[1];  # only print the data part
    close OBS;
    # the headers are stripped from obsDat and simDat below
    print STDERR "INFO: running $rejExe $tmpObs $tmpSimDat2 $tol $columns >> $tmpSimDat\n";
    my $rc = system ("$rejExe $tmpObs $tmpSimDat2 $tol $columns >> $tmpSimDat");
    unless ($rc == 0) {
	print STDERR "ERROR: $rejExe ran funny: $?\n";
	print STDERR "ERROR: Did you make sure that tolerance ($tol) * ".
	    "(number of sim runs) is\nERROR: greater than 100 or so?\n";
	die;
    }
    # checking that there are decent number of accepted points
    $rc = `wc -l < $tmpSimDat`;
    die "wc $tmpSimDat failed: $?\n" if $?;
    chomp $rc;
    if ($rc < 100) {
	warn "WARN: With tolerance of $tol, there are only $rc sampling " .
	    "points.\nWARN: If you encouter a problem, you may need to run ".
	    "more simulations\nWARN: or use a higher tolerance (-t option).\n";
    }
} else {  # Not using preprocessing by rejection
    open SIMDAT, "<$simDat" || die "Can't open $simDat\n";
    my $savedHeader = "";
    while (<SIMDAT>) {
	if ($. != 1) {
	    if (/PRI\.numTauClass/) {  # multiple concatenated sims
		if ($_ eq $savedHeader) { # ignore the line
		    next;
		} else {
		    die "ERROR: You concatenated the multiple simulation ".
			"results.  But the headers do not match (indicating ".
			"that the simulations weren't run with the same ".
			"parameters). 1st header was\n\n$savedHeader\n".
			"The header at line $. was\n\n$_\n";
		}
	    }
	} else {  # 1st line, prepare the header for the output of rejection 
	    $savedHeader= $_;  # used to check multiple headers when
	                       # multipe simulations are concatenated
	}
	print $tmpSimDatfh "$_";
	$tmpSimDatfh->close();
    }
    close SIMDAT;
}

# prepare the obsDat for the acceptRej.r
open OBS, ">$tmpObs" || die "Can't open temporary obs file $tmpObs\n";
#for (my $extra = 0; $extra < $numExtraPrior; $extra++) { # header padding
#    print OBS "$priorNames[$extra]\t";
#}
print OBS $obsDataArray[0];  # only print the data part
#for (my $extra = 0; $extra < $numExtraPrior; $extra++) { # data padding
#    print OBS "0\t";
#}
print OBS $obsDataArray[1];  # only print the data part
close OBS;

# The following is an attempt to implement Bayes Factor in perl
# But it's now done in R.
# my @critVals = (0.01, 0.05, 0.1);
# my @tmpCritVals =  ();
# foreach my $cc (@critVals) {
#     push  @tmpCritVals, ($cc) x $numPriorCols;
# }

# return the proportion of values below the threshold
# used to calculate bayes factor
#my @priorLT = 
#    FreqOfValuesLessThan([1..$numPriorCols],\@tmpCritVals, $simDat, 1);

# run R
my @output = `$Rbin --quiet --no-save --no-restore --slave < $tmpR`;

print @output;

exit (0);

# find all required R scripts, the main R script has to source several
# R scripts within.  This function will return a text string of the R main script after
# a modification so that source() points to the correct file
sub ProperMainRScript {
    my $result = "";
    my $mainR = FindFile($mainRscript);
    if ($mainR eq '-1') {
	die "Can't find $mainRscript in directories:\n", join(":", @INC), "\n";
    }

    my $make_pd = FindFile($make_pdRscript);
    if ($make_pd eq '-1') {
	die "Can't find $make_pdRscript in directories:\n",join(":", @INC),"\n";
    }

    my $loc2plot = FindFile($loc2plotRscript);
    if ($loc2plot eq '-1') {
	die "Can't find $loc2plotRscript in directories:\n",join(":", @INC),"\n";
    }

    my $calmod = FindFile($calmodRscript);
    if ($calmod eq '-1') {
	die "Can't find $calmodRscript in directories:\n",join(":", @INC),"\n";
    }

    $result = "source(\"$make_pd\")\nsource(\"$loc2plot\")\nsource(\"$calmod\")\n";

    open MAIN_R_TMPL, "<$mainR";
    while(<MAIN_R_TMPL>) {
	s/^\s*source\s*\(\s*["']$make_pdRscript['"]\s*\)\s*$//;
	s/^\s*source\s*\(\s*["']$loc2plotRscript['"]\s*\)\s*$//;
	s/^\s*source\s*\(\s*["']$calmodRscript['"]\s*\)\s*$//;
	$result .= $_;
    }
    close MAIN_R_TMPL;

    return $result;
}

# Making R script to print out stat names
sub PrintSumStatNames {
    my $sumStatExe = FindFile($sumStatExe);
    if ($sumStatExe eq '-1') {
	die  "ERROR: Can't find $sumStatExe, is it installed in your PATH?\n";
    }

    system("$sumStatExe -n");
    return;
}

# making the R script
sub MkStdAnalysisRScript {
    my $fh = shift;
    
    my $mainRScript = ProperMainRScript();
    print $fh "$mainRScript\n";
    
    if($ROnlyAnalysis) {
      print $fh "res <- stdAnalysis(\"$tmpObs\", \"$tmpSimDat\", pdf.outfile=\"$pdfOut\",posterior.tbl.file=\"$acceptedFile\",pre.rejected=F";
      print $fh ", tol=$tol";
    } else {
      print $fh "res <- stdAnalysis(\"$tmpObs\", \"$tmpSimDat\", \"$tmpPrior\",pdf.outfile=\"$pdfOut\",posterior.tbl.file=\"$acceptedFile\",pre.rejected=T";
      print $fh ", tol=1";
    }
    
    print $fh ", used.stats=c($statString)";

    if (defined ($opt_i)) {
	print $fh ", return.res=T";
    }

    if (defined($opt_r)) {
	print $fh ", rejmethod=T";  # no regression
    } else {
	print $fh ", rejmethod=F";  # regression
    }

    print $tmpRfh ")\n";

    return;
}

sub MkStatString {
    my $ss = shift;
    $ss =~ s/^\s+//; $ss =~ s/\s+$//;
    my @arr = split /\s*,\s*/, $ss;
    my @result =();
    foreach my $ssName (@arr) {
	$ssName =~ s/^['"]//;
	$ssName =~ s/['"]$//;
	if ($ssName =~ /['"]/) {
	    die "ERROR: -s $ss was given, make sure the summary stats names ".
		"are correct, and comma delimited (e.g -s 'pi,pi.b,wattTheta')\n";
	}
	$ssName =~ s/^\s+//; $ssName =~ s/\s+$//;
	next if ($ssName =~ /^$/);
	push @result, "\"$ssName\"";
    }
    @result = Unique(@result);
    $ss = join ",", @result;
    return $ss;
}

# Receive an array, and extract unique elements
sub Unique {
    my %seen = ();
    my @uniq = grep { ! $seen{$_} ++ } @_;
    return @uniq;
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

# The following function is not used, it was an attempt to implement
# Bayes Factor.  But BF is done in acceptRej.r now.

## Read in a tab-delimited text file.
## Then it returns the array whose elements are the frequencies
## of values in columns less than the critical values.
## The column numbers are specified by an array, and the reference to
## the array should be passed as the first argument.
## 2nd argument is the reference to the array of critical values.
## If the two arrays of column numbers and critical values have different
## length, the shorter ones are recycled in the way similar to R.
## Example:
##   FreqOfValuesLessThan([1,3],[0.5,0.5,1,1], "infile.txt")
## returns an array of
##   (freqs of values less than 0.5 in 1st column,
##    freqs of values less than 0.5 in 2nd column,
##    freqs of values less than 1   in 1st column,
##    freqs of values less than 1   in 2nd column)
sub FreqOfValuesLessThan {
    my($colNumArrRef, $valueArrRef, $filename, $header) = @_;

    open IN, "<$filename" || die "Can't open $filename\n";

    my $numCol = scalar(@$colNumArrRef);
    my $numCritVal = scalar(@$valueArrRef);

    # implement R style index recycling
    my @colNums = ();
    my @critVals = ();
    if ($numCol < $numCritVal) {
        for (my $i = 0; $i < $numCritVal; $i++) {
            push @colNums, $$colNumArrRef[$i % $numCol];
        }
        @critVals = @$valueArrRef;
        $numCol = $numCritVal;
    } elsif ($numCritVal < $numCol) {
        $numCol = $numCol;
        @colNums = @$colNumArrRef;
        for (my $i = 0; $i < $numCol; $i++) {
            push @critVals, $$valueArrRef[$i % $numCritVal];
        }
    } else {
        @colNums = @$colNumArrRef;
        @critVals = @$valueArrRef;
    }

    my @result = map {0} 1..$numCol;

    if ($header) {  # get rid of the first line
	my $firstLine = <IN>;
    }
    my $cntr = 0;
    while(<IN>) {
        chomp;
        my @line = split /\t/;
        for( my $i = 0; $i < $numCol; $i++) {
            my $index = $colNums[$i] - 1;  # change to 0-offset
            if ($line[$index] < $critVals[$i]) {
                $result[$i] ++;
            }
        }
        $cntr ++;
    }
    close (IN);
    return map {$_ / $cntr} @result;
}

# find the names of Prior columns and sumstats
# Read in the 1st line of simDat. The first several items contains  prior columns "PRI."
# Then the summary statistics are given.  This function separate these,
# and return two arrays of column names for prior and summary stats
sub GetPriorSumStatNames {
    my $simDatFileName = shift;
    my @priorNames = ();

    # get the header
    open IN, "<$simDatFileName" || die "Can't open $simDatFileName\n";
    my $line1 = <IN>;
    close IN;
    
    $line1 =~ s/^\s+//;  $line1 =~ s/\s+$//;
    my @header = split /\t/, $line1;

    if ($header[0] =~ /^[.\d]+$/) {
	die "ERROR: The first line of the input file should contain header.\n".
	    "You may have created the files with older version of msBayes.pl\n".
	    "and/or obsSumStats.pl, which is incompatible with this version\n".
	    "of acceptRej.pl.  Please redo the analysis with the most recent\n".
	    "version of programs\n"
    }
    while (my $name = shift @header) {
	if ($name =~ /^PRI\./) {
	    push @priorNames, $name;
	} else {  # the rest should be summary statistics
	    unshift @header, $name;
	    last;
	}
    }

    # Find out number of taxon pairs, and some sanity check.
    my @sumStatNames = ();  # extracting only the summary statistics name part, not used now
    my $prevSS = "";
    my $prevCnt = 0;
    my $maxTaxonID = -1;
    for my $index (0..$#header) {
	my $ss  = $header[$index];
	if ($ss =~ /^(.+)\.(\d+)$/) {
	    if ($prevSS ne $1) {
		push @sumStatNames, $1;
		$prevSS = $1;
		if ($2 != 1) {
		    die "ERROR: Header of summary stats in simulated file weird.  " .
			"$ss should end with '.1'\n";
		} elsif (@sumStatNames > 1 && $prevCnt != $maxTaxonID) {
		    die "ERROR: Previous to $ss had $prevCnt, but should be $maxTaxonID\n";
		} elsif ($index == $#header) {  # the very last element
		    if ($prevCnt != -1 && $prevCnt != 1) {
			die "ERROR: Check the last element of $line1\n";
		    }
		} else {
		    $prevCnt = 1;  # everything seems to be ok, so reset
		    $maxTaxonID = $2 if ($maxTaxonID < $2);
		}
	    } else {
		if ($maxTaxonID < $2){
		    # This should happen only in the first kind of summary stats.  
		    # If not, number of taxon pairs differ for different summayr stats.
		    if (@sumStatNames != 1) {
			die "ERROR: Header of summary stats in simulated file weird ($ss)\n";
		    }
		    $maxTaxonID = $2;
		}
		if ($prevCnt + 1 != $2) {
		    die "ERROR: Header of summary stats in simulated file weird.  ".
			"the number in $ss should be " . ($prevCnt + 1);
		}
		if ($index == $#header && $2 != $maxTaxonID) {  # final element
			die "ERROR: check the final element of $line1, ".
			    "should be $maxTaxonID\n";
		}
		$prevCnt++;
	    }
	} else {  # should never come here
	    die "ERROR: $ss doesn't have the correct naming scheme: ss.number\n";
	}
    }

    return (\@priorNames, \@header, $maxTaxonID);
}


# Find the elements in @array which matches $target.digits and return the array of index
# Returns 0-offset index array
#  Argument: ($target, @array)
sub FindMatchingIndex {
    my $target = shift;
    $target =~ s/^["']//;
    $target =~ s/["']$//;
    my @result = ();
    for (my $i = 0; $i < @_; $i++) {
	push (@result, $i) if ($_[$i] =~ /^$target\.\d+$/);
    }
    return (@result);
}

# Find the number of columns in tab delimited text file.
sub ColNumTabDelimFile {
  my $simDat = shift;
  my ($line, $numColInFile);
  open IN, "<$simDat" || die "Can't open $simDat\n";
  while (defined ($line = <IN>)) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @a = split /\t/, $line;
    $numColInfile = @a;
    last;
  }
  close IN;
  return $numColInfile;
}


# Return an array with $what repeated $times times.
sub Rep {
    my ($what, $times) = @_;

    my @result = ();
    foreach my $i (0..($times-1)) {
	push @result, $what;
    }
    return @result;
}
