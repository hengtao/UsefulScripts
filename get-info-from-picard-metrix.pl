#!/usr/bin/perl -w
use warnings;
use File::Glob;
use Getopt::Long;
use Encode;
use File::Path;
use File::Basename;

my ($inputdir, $outfile, $help);
GetOptions(
    "inputdir=s"      => \$inputdir,
    "outfile=s"      => \$outfile,
    "help|h!"        => \$help
);
if ($help)
{
    usage();
    exit(1);
}

&info_from_picard_metrix($inputdir, $outfile);

sub info_from_picard_metrix{
    my ($inputdir, $outfile) = @_;
    open(OUT, ">$outfile") or die "cannot write into this file:$!\n";
    print OUT "sampleID\tCov1X\tCov5X\tCov10X\n";
    my @metrixs = glob("$inputdir/*metrix");
    foreach my $metrix (@metrixs){
    my ($base,$path,$type) = fileparse($metrix);
    open IN, $metrix or die "cannot open this file:$!\n";
    while(<IN>){
        chomp;
        if(/GENOME_TERRITORY/){
            my $line = <IN>;
            chomp $line;
            my @array = split(/\t/, $line);
            my $Cov1X = &keep_decimal($array[12]);
            my $Cov5X = &keep_decimal($array[13]);
            my $Cov10X = &keep_decimal($array[14]);
            print OUT "$base\t$Cov1X\t$Cov5X\t$Cov10X\n";
        }
    }
    close IN;
  }
  close OUT;
}

sub keep_decimal{
    sprintf( "%.5f", shift ) * 100 . '%';
}
sub usage{
    print << "EOF";

This script was used to evaluate the quality of the raw data\n
Usage:\n
perl runclean-01.pl     --inputdir      <inputdir contains metrix files, suffix of file is "metrix">
                        --outfile      <formatted output file>
EOF
}


