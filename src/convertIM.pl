#!/usr/bin/perl -w

# convertIM.pl
#
# Copyright (C) 2009 Wen Huang, Naoki Takebayashi and Michael Hickerson
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

# This is the default output filename
my $batchFileName = "batch.masterIn.fromIM";
# to avoid mess, all fasta created from IM will be stored in this directory
my $fastaDir = "fastaFromIM";  

my $defaultMutScaler = 20;  # with shrimp pi_b it is about 19.71.

my $usage="Usage: $0 [-h] [-o batchFileName] [-m mutationRateMultiplier] [-s substModel] imfileListFile\n".
    "  -h: help\n".
    "  -o: specify the output file (configFile for msbayes.pl)\n".
    "      default filename: $batchFileName\n" .
    "  -m: If the theta scaler of IM file is set to 0.25, it assumes animal\n".
    "      mtDNA.  The mutation rate is higher, and the value given to this\n".
    "      option is used as the mutation rate scaler (4-th column of\n".
    "      SAMPLE_TBL) for mtDNA. default value is $defaultMutScaler.\n".
    "  -s: specify mutational model\n".
    "        Argument:\n".
    "          JC: Jukes Cantor model\n".
    "          F81: equal input model\n".
    "      By default, it creastes the configFile which uses equal-input\n".
    "      model (= F81).  To use HKY model, you need to estimate transition\n".
    "      transversion rate ratio for each locus of each taxon pair, and change \n".
    "      the values of the 7-th column (tstv) of configFile for msbayes.pl (output\n".
    "      of this script).\n";

#use strict;  # need to check the exact syntax
use File::Copy;
use IO::File;

use Getopt::Std;

our ($opt_h, $opt_o);

getopts('ho:m:s:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

if(defined($opt_o)) {
    $batchFileName = $opt_o;
}

$defaultMutScaler = $opt_m if (defined($opt_m));

my $mutModel="F81";
if (defined($opt_s)) {
    $mutModel = "JC" if ($opt_s =~ /^JC$/i);
}

my $obsSS = "obsSumStats.pl";
my $tstv = 1;
my $verbose = 1;  # prints extra info in stderr
my $sep = "\t";

open LISTFILE, "< $ARGV[0]" ||
   die "Can't open file $ARGV[0] : $! \n";    

my @fileList = <LISTFILE>;
close LISTFILE;
chomp(@fileList);

# This will create an empty directory to store the fasta files
CheckNBackupFile("$fastaDir", 'dir');

my $warnFlag = 1;

# create batch.masterIn file
CheckNBackupFile($batchFileName, 'file');
open BATCH, "> $batchFileName" || 
    die "Can't open file $batchFileName : $! \n";    

print BATCH "BEGIN SAMPLE_TBL\n";

foreach my $currentFile (@fileList)
{ 
   if ($warnFlag) {   # Maybe we should handle line ending better at some point
       if ($currentFile =~ /\r\n$/) {
	   warn "WARN: In convertIM.pl, the file $ARGV[0] appears to have DOS style line ending, it may not work correctly if it's not converted to unix style newline characters\n";
	   $warnFlag = 0;
       }
   }

   # trim the line
   $currentFile =~ s/#.*//;    # remove comments
   $currentFile =~ s/^\s+//;   # remove leading whitespace
   $currentFile =~ s/\s+$//;   # remove trailing whitespace
   
   next if ($currentFile =~ /^\s*$/); # skip empty line   

   # check whether currentFile is OK
   die "ERROR: $currentFile is not readable\n" unless (-r $currentFile);
   die "ERROR: $currentFile is empty\n" if (-z $currentFile);
   die "ERROR: $currentFile is not a text file\n" unless (-T $currentFile);
      
   my ($locusName, 
       $number_pop1,$number_pop2, $sample_length, 
       $mutation_model,$NScaler, $mutScaler, $total);
   my $lineRead = 0;
   
   # open a .im file, grab the data, then close it (to save resource)
   open IMFILE, "< $currentFile" ||
       die "ERROR: Can't open file $currentFile : $! \n";    
   
   my @imData =  <IMFILE>;
   
   close IMFILE;

   # Maybe we should handle line ending better at some point
   if ($imData[0] =~ /\r\n$/) {
       warn "WARN: In convertIM.pl, the file $currentFile appears to have DOS style line ending, it may not work correctly if it's not converted to unix style newline characters\n";
   }

   shift(@imData); # 1st line contains arbitary text
   
   # remove leading/trailing spaces and newline char from each line.
   map {s/^\s+//; s/\s+$//} @imData;

   # remove empty lines
   @imData = grep { ! /^\s*$/ } @imData;

   my $thisLine = shift(@imData);
   # comment lines starting with # can be between line 1 and 2.
   while ($thisLine =~ /^\s*#.*$/) {
	  $thisLine = shift(@imData);
   }

   # pop names
   my @tmpArr = split /\s+/, $thisLine;
   unless (@tmpArr == 2) {
       RemovePartialFiles();
       die "ERROR: 2nd line of $currentFile should be TWO population names separated by spaces\n" .
	   "$thisLine\n";
   }
   my ($pop1Name, $pop2Name) = @tmpArr;
   
   # number of loci
   $thisLine = shift (@imData);
   if ($thisLine !~ /^\d+$/) {
       RemovePartialFiles();
       die "ERROR: 3rd line of $currentFile should be integer (number of Loci)\n";
   }
   my $numLoci = $thisLine;

   # Overall headers are done processing
   # So start to deal with each locus
   
   for my $ll (1..$numLoci) {
       my $thisLocusHeader = shift(@imData);
       my @info = split /\s+/, $thisLocusHeader;
       if (@info < 6) {
	   RemovePartialFiles();
	   die "ERROR: Reading the locus info of $currentFile, it".
	       " should contain at least 6 elements:\n$thisLocusHeader\n";
       }
       
       ($locusName, $number_pop1, $number_pop2, 
	$sample_length, $mutation_model, $NScaler) = @info;

       warn "\n\n### Processing file=$currentFile, $locusName ###\n";

       my $numSamples = $number_pop1 + $number_pop2;
       
       # Make sure enough lines left for sequence data
       if (@imData < $numSamples) {
	   RemovePartialFiles();
	   die "ERROR: In $currentFile, there should be $numSamples samples ".
	       "(= $number_pop1 + $number_pop2)\n".
	       "for locus= $locusName.  But only " . scalar(@imData) .
	       " samples are left in the file.\n";
       }

       # Checking if there are extra lines left after the last locus
       if ($ll == $numLoci) {
	   if (@imData > $numSamples) {
	       warn "WARN: After processing all $ll loci in $currentFile, ".
		   "there are still\n".
		   "WARN: some lines left in the file (see below).\n".
		   "WARN: These extra lines are ignored:\n".
		   join ("\n", @imData) . "\n";
	   }
       }
       
       # a set of sequence is extracted
       my @samples = splice @imData, 0, $numSamples;

       # first 10 characters of each line are the sequence name
       my @alleleNames = map {substr $_, 0, 10} @samples;
       my @seqDat = map {uc(substr $_, 10)} @samples;

       # remove leading/trailing spaces
       # For the name, all spaces are converted to '_'
       # U -> T
       map {s/^\s+//; s/\s+$//; s/\s+/_/g} @alleleNames;
       map {s/^\s+//; s/\s+$//; s/U/T/g} @seqDat;

       # check non-DNA char, and deal with degenarate code
       my @degenArr = ();
       foreach my $thisSeqIndex (0..$#seqDat) {
	   my $degenChar = $seqDat[$thisSeqIndex]; 
	   $degenChar =~ s/[^RYKMSWBDHVN]//g;
	   if (length($degenChar) > 0) {
	       push @degenArr, "$alleleNames[$thisSeqIndex]: $degenChar";
	       # convert degenerate char to '?'
	       $seqDat[$thisSeqIndex] =~ s/[RYKMSWBDHVN]/\?/g;
	   }
	   
	   if ($seqDat[$thisSeqIndex] !~ /^[ATGC\?\-]*$/) {
			$seq = $seqDat[$thisSeqIndex];
			print "$seq\n";
			while($seq =~/(.)/g){
    			if( $1 !~ /^[ATGC\?\-]/ ) {
					print "$1\n";
				}
			}
	       RemovePartialFiles();
	       $seqDat[$thisSeqIndex] =~ s/[ATGC\?\-]//g;
	       die "ERROR: a sequence ($alleleNames[$thisSeqIndex]) in ".
		   "file=$currentFile, locus=$locusName \n".
		   "ERROR: contains non-DNA characters: $seqDat[$thisSeqIndex]\n";
	   }
       }

       if (@degenArr > 0) {
	   warn "INFO: In $currentFile, locus=$locusName,\n".
	       "INFO: following degenerate characters were replaced by '?'\n".
	       "INFO: in the fasta file created by this script:\n".
	       join ("\n", @degenArr) . "\n\n";
       }
       
       # remove gaps
       # make an array with each element = name tab DNA-sequence
       my @gapRmDat = ();
       for my $i (0..$#seqDat) {
	   push @gapRmDat, "$alleleNames[$i]\t$seqDat[$i]";
       }
       # '?' attached to shorter seqs
       @gapRmDat = AdjustSeqLength(@gapRmDat);

       # update the original data with the adjusted sequences
       @seqDat = GetSeqDat(@gapRmDat);

       # removing sites with any gaps
       @gapRmDat = RemoveSitesWithGaps(\@gapRmDat);

       # get rid of seq names before tab
       @gapRmDat = GetSeqDat(@gapRmDat);

       # NOTE: these data (sites with any gaps are removed) are used
       # to calc base composition and seqLen, but fasta files created
       # by this script will retain gaps.  However the data in fasta
       # are modified: (1) all ambig. chars are converted to '?',
       # (2), the sequence length are adjusted by adding ? to the ends
       # if some seqs are shorter.

       # calculate the proportion of A,T,C and G
       my %baseCnter = ("A"=>0,"T"=>0,"G"=>0,"C"=>0);
       foreach my $thisSeq (@gapRmDat) {
	   my %thisBaseComposi = CntBase($thisSeq);
	   foreach my $bb ('A', 'T', 'G', 'C') {
	       $baseCnter{$bb} += $thisBaseComposi{$bb};
	   }
       }
       
       # convert to base frequency
       my $totalSampleLength = $baseCnter{'A'} + $baseCnter{'T'} +
	   $baseCnter{'G'}+$baseCnter{'C'};
       if ($totalSampleLength != 0)
       {
	   foreach my $bb ('A', 'T', 'G', 'C') {
	       $baseCnter{$bb} /= $totalSampleLength;
	   }
       } else {
	   RemovePartialFiles();
	   die "ERROR: in $currentFile, sequence length of 0 encoutered\n";
       }
       
       # construct the header for batch.masterIn file
       # this is the seqLen after removing gaps
       my $newSeqLen = 
	   MyRound($totalSampleLength / ($number_pop1 + $number_pop2));
#       if ($newSeqLen != $sample_length) {
	   # this is already printed in RemoveSitesWithGaps().
#	   warn "INFO: For $currentFile, sites containing an ambiguous char ".
#	       "or a gap was found. " . ($sample_length-$newSeqLen) . 
#	       " sites are removed\n";
#       }

       # removing .im and removing directory names
       my $currentFileBasename;
       if ($currentFile =~ /([^\/]+).(im|phy)$/) {
	   $currentFileBasename = $1;
       } else {
	   $currentFileBasename = $currentFile;
       }
       my $fastaFileName = $currentFileBasename ."_${locusName}.fasta";

       # note the sequence length is after sites with any gaps are removed
       if ($NScaler == 0.25) {
	   $mutScaler = $defaultMutScaler;
	   
	   warn "WARN: in convertIM.pl, For $currentFileBasename:$locusName,".
	       " theta scaler is 0.25, so assuming mutation rate ".
	       "is $mutScaler times higher than nuclear loci.  ".
	       "If this is not the case, correct the 4-th column ".
	       "of $batchFileName.\n" unless (defined ($opt_m));
       } else {
	   if ($NScaler != 1) {
	       warn "WARN: in convertIM.pl, For ".
		   "$currentFileBasename:$locusName, theta scaler is ".
		   "$NScaler, but assuming mutation rate ".
		   "is same as nuclear loci.  ".
		   "If this is not the case, correct the 4-th column ".
		   "of $batchFileName.\n";

	   }
	   $mutScaler = 1;
       }

       if ($mutModel eq "JC") {
	   $baseCnter{A} = $baseCnter{C} = $baseCnter{G} = $baseCnter{T} = 0.25;
       }
       print BATCH
	   join("\t",($currentFileBasename, $locusName, $NScaler, $mutScaler, 
		      $number_pop1, $number_pop2, $tstv, $newSeqLen)) .
		      sprintf ("\t%.3f\t%.3f\t%.3f\t", $baseCnter{A}, 
			       $baseCnter{C}, $baseCnter{G}). 
			       "$fastaDir/$fastaFileName\n";
       
       # create fasta file
       CheckNBackupFile("$fastaDir/$fastaFileName", 'file');
       open FASTA, ">$fastaDir/$fastaFileName" ||
	    die "Can't open file $fastaFileName : $! \n";
       
       for (my $i = 0; $i < @seqDat; $i ++)
       {
	   my $popName = ($i < $number_pop1) ? $pop1Name : $pop2Name;
	   
	   print FASTA ">${currentFileBasename}_${popName}_${locusName}_$alleleNames[$i]\n$seqDat[$i]\n";
       }
       
       close FASTA;
       
   } # for my $ll (1..$numLoci)
} # while (still files in the file list)


print BATCH "END SAMPLE_TBL\n\n";

close BATCH;

my %ssHash = GetSumStats($batchFileName);
my @sampleMat = GetSampleTable($batchFileName);
my @thetaRangeFromObs = EstimateRangeOfThetaPerSite(\%ssHash, \@sampleMat);

# print mutation table at the begining of batch file
my @settings = (
    "# bounds for theta per site (guessed from observed pi within subpops)",
    "upperTheta = $thetaRangeFromObs[1]",
    "lowerTheta = $thetaRangeFromObs[0]",
    "# upper limit of tau (divergence time)",
    "upperTau = 1.0",
    "# number of tau classes (Psi): 0 means Psi are drawn from [1,#taxonPairs]",
    "numTauClasses = 0",
    "# upper bound of migration rate (0 disables migration)",
    "upperMig = 0.0",
    "upperRec = 0.0",
    "# Ancestral theta multiplier:",
    "#  product of this and upperTheta is the upper bound of ancestral theta",
    "upperAncPopSize = 0.25",
    "reps = 1000000",
    "# Most users don't want to constrain the subparameters",
    "constrain = 0",
    "subParamConstrain = 111111111");


#print BATCH "BEGIN CONSTRAIN\n";

my @csData = ("1.0	0.9	0.1	0.5	0.0	10.1	1.5	0.1	0.0",
 	      "1.1	0.8	0.2	0.6	0.0	20.1	1.4	0.2	0.0",
	      "1.2	0.7	0.3	0.7	0.0	30.1	1.3	0.3	0.0",
	      "1.0	0.3	0.7	0.8	0.0	40.1	1.2	0.4	0.0",
	      "1.0	0.3	0.8	0.9	0.0	5.1	1.1	0.5	0.0",
	      "1.0	0.3	0.9	0.3	0.0	25.1	1.0	0.5	0.0");

# opening with update mode to produce the final config file
open (UPDATE, "+<$batchFileName") || die "Can't open $batchFileName\n";

my @sampleTblString = <UPDATE>;
my $output = join("\n", @settings) . "\n\n" . 
    "# taxonName\tlocusName\tNe_Scalar\tMut_Scalar\tsampleSizeA\tsampleSizeB\ttstv\tseqLen\tAfreq\tCfreq\tGfreq\tfastaFileName\n" .
    join("",@sampleTblString) . "\n" .
    "# Most users can ignore the following table\n" .
    "BEGIN CONSTRAIN\n" . join ("\n", @csData) . "\nEND CONSTRAIN\n";

seek(UPDATE, 0, 0) || die "Can't seek to start of $batchFileName: $!";
print UPDATE $output || die "Can't print to $batchFileName: $!";
truncate(UPDATE, tell(UPDATE)) || die "Can't truncate $batchFileName: $!";
close(UPDATE);

print "$batchFileName\n";

exit(0);

sub EstimateRangeOfThetaPerSite {
    my ($sumStatHashRef, $smplTblRef) = @_;

    my $NScalerColumnIndex = 2;
    my $mutScalerColumnIndex = 3;
    
    my $sumStatNameForTheta1 = 'pi.wPop1';
    my $sumStatNameForTheta2 = 'pi.wPop2';
    my $numCol = @{${$smplTblRef}[0]};
        
    my $numTaxonLocus = @$smplTblRef;

    my @thetaPerSite = ();
    for my $i (1..($numTaxonLocus)) {
	my $piw1 = ${$sumStatHashRef}{$sumStatNameForTheta1 . "." . $i};
	my $piw2 = ${$sumStatHashRef}{$sumStatNameForTheta2 . "." . $i};

	my $thetaScaler = ${${$smplTblRef}[$i-1]}[$NScalerColumnIndex] *
	    ${${$smplTblRef}[$i-1]}[$mutScalerColumnIndex];
	push @thetaPerSite, ($piw1 / $thetaScaler, $piw2 / $thetaScaler);
    }

    my $minT = Min(@thetaPerSite);
    $minT = Max($minT * 0.0001, 0.00000000004);
    # larger of the two values are taken. 4*10^-11 is from mu = 10^(-11), Ne=1

    return ($minT, Max(@thetaPerSite) * 4);
}

# The argument of this function is the filename of msbayes
# configuration file (SAMPLE_TBL is all it needs).
# It reads in SAMPLE_TBL section, and create two dimensional array (matrix)
# and return this array.
# So $result[0][0] is the taxon-pair name of the 1st pair, $result[0][1]
# is the gene name, $result[0][2] is the NScaler.
# $result[1][] is the taxon-pair name of the 2nd pair.
sub GetSampleTable {
    my $confFile = shift;

    open CONFIN, "<$confFile" || 
	die "ERROR: Can't open $confFile in GetSampleTable\n";

    my $conf = "";
    while (<CONFIN>) {
	s/#.*$//;  # remove any comments
	next if (/^\s*$/);
	$conf = $conf . $_;
    }
    close CONFIN;

    my $sampleTbl = "";
    if ($conf =~ /\s*BEGIN\s+SAMPLE_TBL\s*\n(.+)END\s+SAMPLE_TBL\s*\n/s) {
	# /s means . match newline
	$sampleTbl = $1;
    } else  {
	die "ERROR: BEGIN SAMPLE_TBL and END SAMPLE_TBL not found in " .
	    "$$confFile in convertIM.pl\n";
    }

    my @tmpTbl = split "\n", $sampleTbl;
    
    # make two dimentional array
    my @result = ();
    foreach my $i (@tmpTbl) {
	$i =~ s/^\s+//; $i =~ s/^\s+$//;
	my @tmp1row = split /\t/, $i;
	push @result, \@tmp1row;
    }

    return (@result);
}

# The argument of this function is the filename of msbayes
# configuration file (SAMPLE_TBL is all it needs).
# Run obsSumStats.pl with -s 0 option, and return a hash which contains
# all of the info.
# Returned hash has the keys like pi.b.1, pi.b.2, ..., pi.w.1, pi.w.2, ...
# and the corresponding summary statistics values.
sub GetSumStats {
    my $confFile = shift;

    die "ERROR: $confFile is not readable\n" unless (-r $confFile);
    die "ERROR: $confFile is empty\n" if (-z $confFile);
    die "ERROR: $confFile is not a text file\n" unless (-T $confFile);
    
    my $obsSSbin = FindExec($obsSS);

    my @sumStats = `$obsSSbin -s 0 $confFile 2> /dev/null`;
    
    if (@sumStats != 2) {
	die "ERROR: hmmm, in GetSumStats of convertIM.pl, obsSumStats.pl ".
	    "returned weird results\n" . join ("", @sumStats);
    }

    my @statNames = split /\t/, $sumStats[0];
    my @statValues = split /\t/, $sumStats[1];
    if (@statNames != @statValues || @statNames < 1) {
	die "ERROR: hmmm, in GetSumStats of convertIM.pl, obsSumStats.pl ".
	    "returned weird results\n" . join ("", @sumStats);
    }
    
    my %result = ();
    foreach my $i (0..$#statNames) {
	if (defined($result{$statNames[$i]})) {
	    die "ERROR: $result{$statNames[$i]} occurs multiple times.  ".
		"This is a programmer's error, please notify the developper\n";
	}
	$result{$statNames[$i]} = $statValues[$i];
    }
    return %result;
}

####### following functions are common with obsSumStats.pl ###############
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
        warn ("INFO: Following " . scalar(@gapSites) . 
	      " site(s) contain gap(s), removed from analysis:\n");
        print STDERR join(" ", @gapSites);
        print STDERR "\n";
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
        push @result, $line[1];
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

####### end of functions common with obsSumStats.pl ###############


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

# takes an array, and find the max value
sub Max {
    my $max = shift;
    foreach my $i (@_) {
	$max = $i if ($i > $max);
    }
    return $max;
}

sub Min {
    my $min = shift;
    foreach my $i (@_) {
	$min = $i if ($i < $min);
    }
    return $min;
}

sub CntBase {
  my $seq = shift;
  my @baseArray = split (//, $seq);
  my %cntResult = (
    "A" => 0,
    "T" => 0,
    "G" => 0,
    "C" => 0
  );

  foreach my $b (@baseArray) {
    $cntResult{$b}++;
  }
  return %cntResult;
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



# rounding function x.y becomes x when y < 5, x+1 otherwise
sub MyRound {
    my $num = shift;

    return int($num + 0.5 * ($num <=> 0));
}

# When something fails, clean-up the partial files
sub RemovePartialFiles {
    warn "REMOVING PARTIAL FILES: $batchFileName & $fastaDir\n";
    system ("rm -rf $fastaDir");
    unlink $batchFileName || die "Can't delete $batchFileName: $!\n";
}
