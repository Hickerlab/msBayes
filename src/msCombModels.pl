#!/usr/bin/perl -w

use strict;

my $usage = "\nUsage: $0 [-h] [-l 'model0,model1,...'] file1 [file2...]\n\n".
    "  -h: help\n".
    "  -l: model labels (names), comma separated string enclosed in single quotes\n\n".
    "This scripts combine several prior files created by simulations.\n".
    "For example, it can be used when you run simulations under several\n".
    "differnt types of models and you want to do the model selection.\n".
    "It assumes that each file corresponds to a different model and it adds\n".
    "an additional 'PRI.model' column as the first column (tab delimited).\n".
    "\n".
    "If -l option isn't given, the model names automatically are assigned\n". 
    "to be sequential numbers (1, 2, 3, ...).\n".
    "-l option takes a text string separated by comma. For example,\n".
    "msCombModels.pl -l 'modelA,modelB,modelC' file1 file2 file3 > output\n".
    "will label the first column with modelA, modelB and modelC instead of\n".
    "1, 2 and 3.  The number of elements in -l option should be equal to\n".
    "the number of files given.";

use Getopt::Std;

our ($opt_h, $opt_l);
getopts('hl:') || die "$usage\n";

die "$usage\n" if (@ARGV < 1 || defined ($opt_h));

my @modelLabel = 1..scalar(@ARGV);
if (defined($opt_l)) {
    $opt_l =~ s/^\s+//;
    $opt_l =~ s/\s+$//;
    @modelLabel = split /\s*,\s*/, $opt_l;
    if (@modelLabel != @ARGV) {
	warn "$usage\n";
	die "\nERROR: number of model label (-l) = " . scalar(@modelLabel) .
	    " doesn't match the number of files: " . scalar(@ARGV) . "\n";
    }
}

my $fileNum = 1;
while(my $file = shift @ARGV) {

    open IN, "<$file" || die "Can't open $file\n";

    my $model=$modelLabel[$fileNum - 1];
    warn "INFO: Model name for a file: $file = \"$model\"\n";

    while(<IN>) {
	chomp;
	if ($. == 1) {
	    if ($fileNum == 1) {
		print "PRI.model\t$_\n";
	    } else {
		next;
	    }
	} else {
	    print "$model\t$_\n";
	}
    }
    
    close IN;
    $fileNum++;
}
