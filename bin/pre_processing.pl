#!/usr/bin/env perl
use FindBin;                # Import FindBin to determine script's directory
use lib "$FindBin::Bin/bin";
use Parallel::ForkManager;
use HMMmodE;
my $max_processors = 24;    ### Give the number of processors to run the job
my $fork           = new Parallel::ForkManager($max_processors);
my $file           = shift;
open( LIST, $file );
my @file = <LIST>;
close LIST;

foreach (@file) {
    $fork->start and next;    # do the fork
    chomp $_;

    #system "formatdb -p T -i $_ -o T";
    system "makeblastdb -in $_ -dbtype prot -parse_seqids";

    #system "blastall -i $_ -p blastp -d $_ -o $_.blast";
    system "blastp -query $_ -db $_ -out $_.blast";
    $fork->finish;            # do the exit in the child process
}
$fork->wait_all_children;
foreach (@file) {
    chomp $_;
    system "mcxdeblast --score=b --sort=a --bcut=5 $_.blast";
    system "mcxassemble -b $_.blast -q -r max --map"
        ;    # transform raw cooccurrence data to mcl matrix format
             # -b is the base name of files to be processed and output
     # -q Warning options, we turn these on if we expect the raw data to be free of repeats

    system
        "mcl $_.blast.sym -scheme 7 -I 1.2 -te 24 -o $_.clus --append-log=yes -V all"
        ;    # performing mcl clustering of blast.sym files
             # -scheme use a preset resource scheme,high schemes result/
     # /in more n more expensive computations that may possibly be more accurate,defualt is 4.
    system "rm $_.blast $_*.p* $_*.err $_*.map $_*.raw $_*.hdr";
    my $clusfile   = $_ . ".clus";
    my $tabfile    = $_ . ".blast.tab";
    my $clusterObj = HMMmodE->new(
        -file     => $_,
        -clusType => 'MCL',
        -clusFile => $clusfile,
        -tabFile  => $tabfile
    );
    my $ref_array = $clusterObj->nsCluster();

    #   print "ns clusters @$ref_array \n";
    my @ns = @$ref_array;
    $ns_hash{$_} = \@ns;
    foreach my $clu_no (@ns) {
        my $nameref = $clusterObj->{'sample_outfile'}{$clu_no};

        #                print ".... $nameref\n";
    }
}
END { warn time - $^T, " seconds elapsed\n" }
