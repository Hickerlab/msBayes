#!/usr/bin/perl -w

use IO::File;
use locale;

my @ssNameVect=
  ( "pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom",
  "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1",
  "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between",
  "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2",
  "S_x", "S_y","S_xy","JWakeley_Psi");

my $obsSS = "obsSumStats.pl";
# `chmod a+x $obsSS`;

if (@ARGV != 1) {
    die "ERROR: This script requires only one argument";
}
my $batchFile = $ARGV[0];

#store taxa_gene order in the config file in %indexHash
#also all available taxa and loci
my %indexHash;
my @smplTbl = GetSampleTable($batchFile);

my @taxa = ();
my @loci = ();
for $i(0..$#smplTbl){
    push @taxa, uc($smplTbl[$i][0]);
    push @loci, uc($smplTbl[$i][1]);
    my $taxaLoci = uc($smplTbl[$i][0])."+".uc($smplTbl[$i][1]);
    $indexHash{$taxaLoci} = $i+1;
}

my %seenTaxa = ();
my @uniqueTaxa = grep { ! $seenTaxa{$_} ++ } @taxa;
my %seenLoci = ();
my @uniqueLoci = grep { ! $seenLoci{$_} ++ } @loci;

#store summary statistics in %ssHash
my %ssHash = GetSumStats($batchFile);

#print suumary statistics table to a file,
#one summary statistics table one table
foreach my $ss (@ssNameVect)
{  
   my $header = "$ss\n\t";
   $header .= join "\t",@uniqueLoci;
   $header .= "\n";
   print $header;
   
   foreach my $taxa (@uniqueTaxa)
   {  
      print "$taxa\t";
      foreach my $gene (@uniqueLoci)
      {
         my $taxaGene = $taxa."+".$gene;
         if (exists $indexHash{$taxaGene})
         {  
	    my $ssKey = $ss.".".$indexHash{$taxaGene};
            print $ssHash{$ssKey}."\t";
         }else
         {  
            print "NA\t";
         }
      }
      print "\n";
   }
   print "\n\n";
}

exit(0);

# The argument of this function is the filename of msbayes
# configuration file (SAMPLE_TBL is all it needs).
# It reads in SAMPLE_TBL section, and create two dimensional array (matrix)
# and return this array.
# So $result[0][0] is the taxon-pair name of the 1st pair, $result[0][1]
# is the gene name, $result[0][2] is the thetaScaler.
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
