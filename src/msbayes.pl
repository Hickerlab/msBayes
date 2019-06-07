#!/usr/bin/perl -w

# msbayes.pl
#
# Copyright (C) 2009  Wen Huang, Naoki Takebayashi and Michael Hickerson
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

my $usage="Usage: $0 [-hd] [-b observedSummaryStatistics] [-r numSims] [-c config] [-i IMconfig] [-o outputFileName] [-s sumStatsSortPattern] [-S seed]\n".
    "  -h: help\n".
    "  -r: number of repetitions\n".
    "  -c: configuration file for msprior.  Parameters setup interactively,\n".
    "      if this option is not specified\n" .
    "  -i: IM format configuration file for msprior,\n" .
    "      it will be converted so that msprior can read \n" .      
    "  -o: output file name.  If not specified, output is STDOUT\n" .
    "  -b: observed summary statistics file name.\n".
    "      If not specified, it won't calculate summary statisctics \n".
    "  -s: Specify how to sort the summary statistics.\n".
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
    "  -S: set the initial seed (but not verbose like -d)\n" .
    "      By default (without -s), unique seed is automaically set from time\n".
    "  -d: debug (msprior and msDQH uses the same initial seed = 1)\n";

use strict;

use File::Copy;
use IO::File;
use File::Temp qw(tempfile);
use IPC::Open2;

use Getopt::Std;

our ($opt_h, $opt_d, $opt_o, $opt_b, $opt_c, $opt_i, $opt_r, $opt_S, $opt_s);
getopts('hdo:b:c:i:r:S:s:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $defaultOutFile = "Prior_SumStat_Outfile";

my $batchFile;

my $debug=0;
my $rmTempFiles = 1; # if this is set to 0, temp files will NOT be deleted
if (defined($opt_d)) {
    $debug = 1;
    $rmTempFiles = 0;
}

my $options = "";

if($debug) {  # force msprior to use the same seed
    $options = "-d 1 ";
}

if (defined($opt_r)) {
    if ($opt_r < 0) {
	die "The number of repetitions should be positive integer\n$usage\n";
    }
    $options = $options . " --reps $opt_r ";
}

if (defined($opt_i) && defined($opt_c)) {
    die "ERROR: Configuaration file can be in only one format. " .
	"Specify EITHER -i OR -c\n";
}

if (defined($opt_i)) {   
    die "ERROR: $opt_i is not readable\n" unless (-r $opt_i);
    die "ERROR: $opt_i is empty\n" if (-z $opt_i);
    die "ERROR: $opt_i is not a text file\n" unless (-T $opt_i);
    
    # `chmod a+x convertIM.pl`;
    my $convertIM = FindExec("convertIM.pl");
    
    $batchFile = `$convertIM $opt_i`;
    chomp $batchFile;

    print STDERR "INFO: batchFile : $batchFile \n";
    $options = $options . " --config $batchFile";
}

if (defined($opt_c)) {
    die "ERROR: $opt_c is not readable\n" unless (-r $opt_c);
    die "ERROR: $opt_c is empty\n" if (-z $opt_c);
    die "ERROR: $opt_c is not a text file\n" unless (-T $opt_c);
    $batchFile = $opt_c;
    if ($debug) {
	print STDERR "INFO: msprior options are: $options\n";
    }
}

if (defined($opt_b)) {
    my $obsSS = FindExec("obsSumStats.pl");
    CheckNBackupFile($opt_b, 'file');
    system("$obsSS $batchFile > $opt_b");
}
 
if (defined ($opt_S)) {
    $options = $options . " --seed $opt_S ";
}

my $ssvSortOpt = "-s 1";
if (defined($opt_s)) {
    die "$usage\nERROR: sorting pattern (-s) should be 0 t0 7"
	if ($opt_s != 0 && $opt_s != 1 && $opt_s != 2 && $opt_s != 3 && $opt_s != 4 && $opt_s != 5 && $opt_s != 6 && $opt_s != 7);
    $ssvSortOpt = "-s $opt_s";
} else  {
    $ssvSortOpt = "-s 7";  # default sorting option set to 7
}

if (defined($opt_S) || defined($opt_d)) {  # set the msDQH use the same seeds
    if (defined($opt_S)) {
	srand($opt_S);
    } else {
	srand(1);
    }
}

#### Find programs
my $msprior = FindExec("msprior");
my $msDQH = FindExec("msDQH");
my $sumstatsvector = FindExec("sumstatsvector");

#### setting up final output filename
#### This will contain summary stats from simulations
my $outFile = '-';  # by default, print out to STDOUT
if(defined($opt_o)) {
    $outFile = $opt_o;
    CheckNBackupFile($outFile, 'file');
} 
open FINAL_OUT, ">$outFile" || die "Can't open $outFile";

# obtain configuration from msprior, which read in the config file
# and return the configuration.
# This may invoke the Interactive parameter set up mode.
my $mspriorConfOut = (defined ($opt_c)) ? 
    `$msprior $options --config $opt_c --info` :
    `$msprior $options --info` ;  # getting the config info
my %mspriorConf = ExtractMspriorConf($mspriorConfOut);

# If parameters are interactively set by the user, we need to incorporate
# this info to make an updated config file to do the real msprior run.
# This is a work-around, so users don't have to type in parameters twice.
# Create a new temp file with this updated parameter info here.
my $newMspriorConf = MkNewMspriorBatchConf($mspriorConf{configFile}, 
					   \%mspriorConf);

# open and close a temp file
# This is used to store the new msprior conf file
my ($tmpMspriorConffh, $tmpMspriorConf);
if ($rmTempFiles) {
    ($tmpMspriorConffh, $tmpMspriorConf) = tempfile(UNLINK => 1) ;    
} else {
    ($tmpMspriorConffh, $tmpMspriorConf) = tempfile() ;
    print STDERR "FILE: \$tmpMspriorConf = $tmpMspriorConf\n";
}
print $tmpMspriorConffh "$newMspriorConf";
$tmpMspriorConffh->close();

$options = $options . " --config $tmpMspriorConf ";

open (RAND, "$msprior $options |") || 
    die "Hey buddy, where the hell is \"msprior\"?\n";

my @msOutCache = ();
my @priorCache = ();

my $prepPriorHeader = 1;

my $msCacheSize = $mspriorConf{'numTaxonLocusPairs'} * 50; # ADJUST MULTIPLIER to reduce I/O
# $numTaxonLocusPairs is the total number of taxa:locus pairs.  If
# taxon pair 1 have 3 loci, taxon pair 2 have 4 loci, and taxon pair 3
# have 1 locus, the 2nd term = 8 So $totalNumSims are the number of
# times msDQH will be invoked.

my $counter  = 0;
my $totalNumSims = $mspriorConf{reps} * $mspriorConf{numTaxonLocusPairs};
# This is used to print out the last simulations.

my $headerOpt = " -H ";

#### Overview of what's going on in the following while loop
# - msprior spits out a line of the prior parameters drawn from the prior distn.
# - Each line contains all of the required parameter for 1 msDQH run (= 1 gene 
#   of 1 taxon pair)
# - msDQH will run with the specified parameter.  The output of msDQH is 
#   stored as a single string, and it will be added to @msOutCache array
# - After every n parameter lines (n = sum n_i, where n_i are the number of 
#   genes used for i-th taxon pair), msprior print out some information in the
#   following format:
#  '# TAU_PSI_TBL setting: 0 realizedNumTauClasses: 2 tauTbl:,1.48,7.92 psiTbl:,1,2'
# - This line indicates that one set of simulations are run for total n genes
# - When certain number of runs (specified by $msCacheSize above), sumstats
#   will be run to process the accumulated outputs of msDQH
# - Then the results of sumStats (and prior informations) will be 
#   outputed to the final file
while (<RAND>) {
    s/^\s+//; s/\s+$//; 

    unless (/^# TAU_PSI_TBL/) {
        # When reached here, it is regular parameter lines, so need to run msDQH
	my ($taxonLocusPairID, $taxonID, $locusID, $theta, $gaussTime, $mig, $rec, 
	    $BottleTime, $BottStr1, $BottStr2, 
	    $totSampleNum, $sampleNum1, $sampleNum2, $tstv1, $tstv2, $gamma,
	    $seqLen, $N1, $N2, $Nanc, 
	    $freqA, $freqC, $freqG, $freqT) = split /\s+/;

## the output line of msprior contains following parameters in this order	
#  0 $taxonLocusPairID, $taxonID, $locusID, $theta, $gaussTime, 
#  6 $mig, $rec, $BottleTime, $BottStr1, $BottStr2,
# 11 $totSampleNum, $sampleNum1, $sampleNum2, $tstv1, $tstv2, 
# 16 $gamma, $seqLen, $N1, $N2, $Nanc, 
# 21 $freqA, $freqC, $freqG, $freqT
	
	# The duration of bottleneck after divergence before the pop. growth
	my $durationOfBottleneck = $gaussTime - $BottleTime;
	
	# option for -r was fixed to 0, so changed to $rec, then forcing
	# it to be 0 here
	# $rec = 0;
	
	my $SEED = int(rand(2**32));  # msDQH expect unsigned long, the max val (2**32-1) is chosen here

	# At the bottom of this script, msDQH options are explained.
	my $ms1run = `$msDQH $SEED $totSampleNum 1 -t $theta -Q $tstv1 $freqA $freqC $freqG $freqT -H $gamma -r $rec $seqLen -D 5 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $BottStr1 $N2 $BottStr2 $BottleTime 2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $durationOfBottleneck 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonLocusPairID 1 Nc $Nanc $mspriorConf{numTaxonLocusPairs}`;

	if ($? != 0) { # Returned weird exit status.
	    warn "WARN: following command returned exit status of $?. Contact the developers\n";
	    warn "$msDQH $SEED $totSampleNum 1 -t $theta -Q $tstv1 $freqA $freqC $freqG $freqT -H $gamma -r $rec $seqLen -D 5 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $BottStr1 $N2 $BottStr2 $BottleTime 2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $durationOfBottleneck 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonLocusPairID 1 Nc $Nanc $mspriorConf{numTaxonLocusPairs}\n";
	    next;
	}
	$ms1run = "# taxonID $taxonID locusID $locusID\n" . $ms1run;
	push @msOutCache, $ms1run;
	
	$counter++;
	next;
    }
    
    # When reached here, msprior printed a TAU_PSI_TBL line which looks like:
    # # TAU_PSI_TBL setting: 0 realizedNumTauClass: 3 tauTbl:,8.71,4.92,4.01 psiTbl:,1,1,1
    # This means one set of simulations are done.
    # Processing this line to prepare 'prior columns' for the final output
    my ($tauClassSetting, $realizedNumTauClass);
    if (/setting:\s+(\d+)\s+realizedNumTauClasses:\s+(\d+)\s+tauTbl:,([\d\.,]+)\s+psiTbl:,([\d\.,]+)/) {
	$tauClassSetting = $1;
	$realizedNumTauClass = $2; # this is not actually used
	my @tauTbl = split /\s*,\s*/, $3;
	my @psiTbl = split /\s*,\s*/, $4;
	
	# prep header
	if ($prepPriorHeader){ 
	    my $headString = "PRI.numTauClass";
	    if ($mspriorConf{numTauClasses} > 0) {
		for my $suffix (1..$mspriorConf{numTauClasses}) {
		    $headString .= "\tPRI.Tau.$suffix";
		}
		for my $suffix (1..$mspriorConf{numTauClasses}) {
		    $headString .= "\tPRI.Psi.$suffix";
		}
	    }
	    $headString .= "\tPRI.Psi\tPRI.var.t\tPRI.E.t\tPRI.omega";
	    push @priorCache, $headString;
	    $prepPriorHeader = 0;  # print this only 1 time
	}

	##### PREPARING PRIOR COLUMNS #####
	# prepare the prior columns
	# PRI.numTauClass
	my @tmpPrior = ($mspriorConf{numTauClasses});

	# conversion of tau
	if ($mspriorConf{numTauClasses} > 0) {
	    #PRI.Tau.1 PRi.Tau.2 ... PRI.Tau.numTauClasses
	    #PRI.Psi.1 PRi.Psi.2 ... PRI.Psi.numTauClasses
	    push @tmpPrior, @tauTbl, @psiTbl;
	}
	
	# PRI.Psi PRI.var.t PRI.E.t PRI.omega 
	#  (= #tauClasses, Var, Mean, weirdCV of tau)
	push @tmpPrior, SummarizeTau(\@tauTbl, \@psiTbl);

	push @priorCache, join("\t", @tmpPrior);

    } else {
	die "ERROR: TAU_PSI_TBL line is weird.\n$_\n";
    }

    # Check if it is time to run sumstats
    if (@msOutCache % $msCacheSize == 0 || $counter == $totalNumSims) {
	# getting read write access to sumstatsvector
	my $pid = open2(\*READ_SS, \*WRITE_SS, "$sumstatsvector $headerOpt $ssvSortOpt"); 
	
	$headerOpt = "";  # remove -H, so header is printed only once

	print WRITE_SS "# BEGIN MSBAYES\n".
	    "# numTaxonLocusPairs $mspriorConf{numTaxonLocusPairs} ".
	    "numTaxonPairs $mspriorConf{numTaxonPairs} ".
	    "numLoci $mspriorConf{numLoci}\n";

	for my $index (0..$#msOutCache) {
	    print WRITE_SS "$msOutCache[$index]";
	}

	close(WRITE_SS);  # need to close this to prevent dead-lock
	
	my @ssOut = <READ_SS>;
	close(READ_SS);
	waitpid($pid,0);

	if (@priorCache != @ssOut) {
	    die "ERROR: size of priorCache (". scalar(@priorCache).
		") differ from sumStatsCache(" .scalar(@ssOut)."\n";
	}
	
	# print out prior etc.
	for my $index (0..$#priorCache) {
	    print FINAL_OUT "$priorCache[$index]\t$ssOut[$index]";
	}
	
	@msOutCache = ();  # clear the cache
	@priorCache = ();
    }
}

close RAND;
close FINAL_OUT;

exit(0);

# interactively setting up. Not used anymore.
sub InteractiveSetup {
    my $outFileName;

    # output filname
    print "Output Filename? [Return] to use default of " .
	"\"$defaultOutFile\"\n";
    chomp($outFileName = <STDIN>);
    if($outFileName eq "") {
	$outFileName = $defaultOutFile;
    }
    CheckNBackupFile($outFileName, 'file');

    return ($outFileName)
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

# This is not used any more

# Take filenames of two files, and concatenate them side by side to
# produce the outputfile given as the 3rd argument.  The tab will be
# inserted after each line of the first file.
# If two files do not have the same number of lines, the output file
# have the same length as the files with less lines.  The rest of the
# longer file is ignored.

# Return value
# 1 two files same number of lines, and successfully concatenated
# 2 two files have different number of lines
# -1 Either one file or bot files were empty.  non-empty file is copied as the 
#    output file
sub ColCatFiles {
    my ($infilename1, $infilename2, $outfilename) = @_;

    # check empty files, if empty, copy is enough
    if ( -s $infilename1 && -z $infilename2) { # file 2 empty
	copy($infilename1, $outfilename) ||
	    warn "WARN: copy $infilename1, $outfilename failed";
	return -1;
    } elsif (-z $infilename1 ) { # file 1 empty or both empty
	copy($infilename2, $outfilename) ||
	    warn "WARN: copy $infilename2, $outfilename failed";
	return -1;
    }

    # both infiles are not empty
    my $retval = 1;

    open FILE1, "<$infilename1" || die "Can't open $infilename1\n";
    open FILE2, "<$infilename2" || die "Can't open $infilename2\n";
    open OUT, ">$outfilename" || die "Can't open $outfilename\n";

    my $numLines1 = `wc -l < $infilename1`;
    die "wc failed: $?" if $?;
    chomp $numLines1;
    my $numLines2 = `wc -l < $infilename2`;
    die "wc failed: $?" if $?;
    chomp $numLines2;

    if ($numLines1 != $numLines2) {
	warn "WARN: number of lines differ between $infilename1 ($numLines1 lines) and $infilename2 ($numLines2 lines)\n";
	$retval = 2;
    }
    
    my $maxLines =  ($numLines1 > $numLines2) ? $numLines1 : $numLines2;
    
    for(my $i = 0; $i < $maxLines; $i++) {
	if ($i < $numLines1) {
	    my $line = <FILE1>;
	    chomp $line;
	    print OUT $line;
	}
	print OUT "\t";

	if ($i < $numLines2) {
	    my $line = <FILE2>;
	    chomp $line;
	    print OUT $line;
	}
	print OUT "\n";
    }
    
    close OUT;
    close FILE1;
    close FILE2;
    return $retval;
}


sub ExtractMspriorConf {
    my $mspriorConfOut = shift;

    my %result =();

    my @generalKwdArr = qw(lowerTheta upperTheta upperTau upperMig upperRec upperAncPopSize reps numTauClasses bufferTauClasses concentrationShape concentrationScale constrain subParamConstrain numTaxonLocusPairs numTaxonPairs numLoci prngSeed configFile);

    for my $kkk (@generalKwdArr) {
	if ($mspriorConfOut =~ /\s*$kkk\s*=\s*([^\s\n]+)\s*\n/) {
	    $result{$kkk} = $1;
	} else {
	    die "In:\n $mspriorConfOut\nCan't find $kkk\n";
	}
    }

    my $mutPara;
    if ($mspriorConfOut =~ /## gMutParam ##\s*\n(.+)## gConParam/s) {
	# note s at the end of /regex/ will let . to match \n
	$mutPara = $1;
    } else {
	warn "Couldn't find mutation parameter table";
    }    
    # I'm not using this, but following info can be extrcted
# ### taxon:locus pair ID 1 taxonID 1 (lamarckii) locusID 1 (mt) NScaler 1 mutScaler 1 ###
# numPerTaxa =    15
# sample =        10 5
# tstv =  11.600000  0.000000
# gamma = 999.000000
# seqLen =        614
# freq:A, C, G, T = 0.323000, 0.268000 0.212000 0.197000
# fileName =      lamarckii.fasta
# ### taxon:locus pair ID 2 taxonID 2 (erosa) locusID 1 (mt) NScaler 0.25 mutScaler 10 ###
# numPerTaxa =    16
# sample =        10 6
# tstv =  13.030000  0.000000
# gamma = 999.000000
# seqLen =        614
# freq:A, C, G, T = 0.266000, 0.215000 0.265000 0.254000
# fileName =      erosa.fasta
# ### taxon:locus pair ID 3 taxonID 3 (clandestina) locusID 2 (adh) NScaler 1 mutScaler 1 ###

    return %result;
}

# This takes a filename of a msprior config file and a hash containing
# the msprior config parameters.  From the file, SAMPLE_TBL (and
# CONSTRAIN_TBL if it exists) are extracted.  Then the tables and the
# parameters in the hash table are combined to create a new string
# corresponding to the config file.  parameters in the config file are
# ignored.
sub MkNewMspriorBatchConf {
    my ($oldConfFileName, $mspriorConfHashRef) = @_;
    my %confHash = %$mspriorConfHashRef;

    # the following kwd should match with SetupParams() in setup.c
    my @generalKwdArr = qw(lowerTheta upperTheta upperTau upperMig upperRec upperAncPopSize reps numTauClasses bufferTauClasses concentrationShape concentrationScale constrain subParamConstrain prngSeed);

    open CONFIN, "<$oldConfFileName" || die "Can't open $oldConfFileName\n";
    my $conf = "";
    while (<CONFIN>) {
	s/#.*$//;  # remove any comments
	next if (/^\s*$/);
	$conf = $conf . $_;
    }
    close CONFIN;

    my $newConf = "";
    if ($conf =~ /\s*(BEGIN\s+SAMPLE_TBL\s*\n.+END\s+SAMPLE_TBL)\s*\n/s) {
	# /s means . match newline
	$newConf = $newConf . "$1\n";
    } else  {
	die "ERROR: BEGIN SAMPLE_TBL and END SAMPLE_TBL not found in " .
	    "$oldConfFileName\n";
    }

    if ($conf =~ /\s*(BEGIN\s+CONSTRAIN\s*\n.+END\s+CONSTRAIN)\s*\n/s) {
	$newConf = $newConf . "$1\n";
    }

    foreach my $kwd (@generalKwdArr) {
	if (defined ($$mspriorConfHashRef{$kwd})) {
	    $newConf = "$kwd = $$mspriorConfHashRef{$kwd}\n" . $newConf;
	} else {
	    warn "WEIRD: in MkNewMspriorBatchConf(), confHash doesn't have ".
		"$kwd\n";
	}
    }

    return $newConf;
}

sub SummarizeTau {
    my ($tauArrRef, $cntArrRef)  = @_;
    
    my $numTauClasses = scalar(@$tauArrRef); # num elements = Psi
    
    my ($sum, $ss, $n) = (0,0,0);
    foreach my $index (0..($numTauClasses-1)) {
	$n += $$cntArrRef[$index];
	$sum += $$tauArrRef[$index] * $$cntArrRef[$index];
	$ss += ($$tauArrRef[$index] ** 2) * $$cntArrRef[$index];	
    }
    
    my $mean = $sum / $n;
    my $var = ($n==1) ? 'NA': ($ss -  $n * ($mean ** 2)) / ($n-1); # estimated, or sample var
    my $dispersionIndex = ($n==1 || $mean == 0) ? 'NA': $var/$mean;
    
    return ($numTauClasses, $var, $mean, $dispersionIndex);
}


#### Explanation of msDQH commandline options by Eli
#    system("$msDQH $SEED $totSampleNum 1 -t $theta -Q $tstv1 $freqA $freqC $freqG $freqT -H $gamma -r $rec $seqLen -D 6 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $BottStr1 $N2 $BottStr2 $BottleTime 2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $durationOfBottleneck 1 Nc $Nanc $numTauClasses 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonLocusPairID 1 Nc $Nanc $mspriorConf{numTaxonLocusPairs} | $sumstatsvector -T $mspriorConf{upperTheta} --tempFile $tmpSumStatVectScratch $headerOpt >> $tmpMainOut");

# minor change 9/8/06; $N1 $N1 $N2 $N2 to $N1 $BottStr1 $N2 $BottStr2

# The command line format is the same as Dick Hudson's ms except for
# everything after -D. Everything after -D specifies what happens
# during X number of time intervals. -D 6 2 means 6 intervals starting
# with 2 populations at time 0 (going backwards in time).

# The rest -D is explained using the following template:
#
# -D 6 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $N1 $N2 $N2 $BottleTime \
#  -2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $durationOfBottleneck 1 Nc $Nanc \
#   -$numTauClasses 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonLocusPairID 1 Nc
#    -$Nanc $mspriorConf{numTaxonLocusPairs}
#
# $sampleNum1 $sampleNum2; the 2 sample sizes of the 2 populations

# 0 I $mig;  $mig is migration rate under an island model of migration

# $N1 $N1 $N2 $N2; Relative size pop1 at begening of timestep,
#   Relative size pop1 at end of timestep, Relative size pop2 at
#   begening of timestep, Relative size pop2 at end of timestep

# $BottleTime; length of first time step (begining of bottleneck going
#    backwards in time)

# 2 1 0 0 1; this is the admixture matrix a[][]; this allows
#   population divergence or admixture, or can specify neither occuring
#   during the #time step. In this case nothing happens (2 populations
#   become two #populations)
#   1 0
#   0 1
#   in a[][] this means all of pop1 goes into pop1, and all of pop2 goes
#   into pop 2 (2 populations remain isolated)
#   If we had:
#   2 1 0 1 0 0 1
#   1 0
#   1 0
#   0 1
#   this would mean that at first we have 3 populaations and then all of
#   pop2 fuses with pop1 (going back in time and hence divergence), and
#   pop3 would remain intact

# 0 I $mig;  the migration region of the next time step

# Nc; Nc specifies that all populations have constant size in this
#   next time step

# $BottStr1 $BottStr2; these are the two constant relative sizes of
#   the two populations during this next time step.

# $durationOfBottleneck; this is the length of this next time step (in this case it
#   ends at the divergence time)

# 1 Nc $Nanc $numTauClasses; specifies that the next time step has
#   only one population (divergence) and the population is constant in
#   size through the time step
#     $Nanc; relative size of this ancestral population
#     $numTauClasses; this is the length of the time step, but has a 2nd
#        meaning unrelated to msDQH. the actual value gets passed on to the
#        summary stats program for parameter estimation purposes. The actual
#        value is somewhat arbitray for msDQH because there is only one
#        population remaining going back in time.  The length of the period
#        can be infinite.

# Three more time-steps use the same "1 Nc $Nanc $length" pattern,
# where "length" has a 2nd use

# If one wants to use msDQH independently on the command line, one can
# add "-P" to see what population model is being used.  Example below
#
# ./msDQH 35 1 -t 20.0 -Q 5.25 0.25 0.25 0.25 0.25 -H 999.000000 -r 0 1000 -D 5 2 20 15 0 I 0.000000 0.8 0.05 0.9 0.05 6.03 2 1 0 0 1 0 I 0.000000 Nc 0.05 0.05 0.001 1 Nc 0.42 6 1 Nc 0.42 1000 1 Nc 0.42 1 -P
 
# Output example using "-P"

# In the below example 2 populations (20 and 15 individuals) diverged
# from a common ancestor of relative size 0.42 at the third time step

#./msDQH 35 1 -t 20.0 -Q 5.25 0.25 0.25 0.25 0.25 -H 999.000000 -r 0 1000 -D 5 2 20 15 0 I 0.000000 0.8 0.05 0.9 0.05 6.03 2 1 0 0 1 0 I 0.000000 Nc 0.05 0.05 0.001 1 Nc 0.42 6 1 Nc 0.42 1000 1 Nc 0.42 1 -P 
#./msDQH nsam 35 howmany 1
#  theta 20.00 segsites 0
#seQmut 1 output 0
#tstvAG CT 5.25 5.25, freqACGT 0.25 0.25 0.25 0.25
#gammaHet alpha 999.00
#  r 0.00 f 0.00 tr_len 0.00 nsites 1000
#  Dintn 5 
#    Dint0 npops 2
#      config[] 20 15 
#      Mpattern 0
#      M[][] 0.00 0.00 0.00 0.00 
#      (Nrec_Npast)[] 0.80 0.05 0.90 0.05 
#       tpast 6.03
#    Dint1 npops 2
#      a[][] 1.00 0.00 0.00 1.00 
#      Mpattern 0
#      M[][] 0.00 0.00 0.00 0.00 
#      (Nrec_Npast)[] 0.05 0.05 0.05 0.05 
#       tpast 6.03
#    Dint2 npops 1
#      (Nrec_Npast)[] 0.42 0.42 
#       tpast 12.03
#    Dint3 npops 1
#      (Nrec_Npast)[] 0.42 0.42 
#       tpast 1012.03
#    Dint4 npops 1
#      (Nrec_Npast)[] 0.42 0.42 
#       tpast 1013.03
