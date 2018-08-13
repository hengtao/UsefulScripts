#!/usr/bin/perl -w
use warnings;
use glob;
use Getopt::Long;
use Encode;
use JSON;
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

&info_from_picard_metrix();

sub info_from_picard_metrix{
  my ($inputdir, $outfile) = @_;
  open OUT ">$outfile" or die "cannot write into this file:$!\n";
  print OUT "sampleID\tCov1X\tCov5X\tCov10X\n"
  my @metrixs = glob "$inputdir/*metrix";
  foreach my $metrix(@metrixs){
  	my ($base,$path,$type) = fileparse($metrix);
	open IN $metrix or die "cannot open this file:$!\n";
	while(<IN>){
		chomp;
		if(/GENOME_TERRITORY/){
			<IN>;
			chomp;
			split(/\t/);
			print OUT "$base\t$_[12]\t$_[13]\t$_[14]\n";
		}
	}
	close IN;
  }
  close OUT;
}

sub usage{
	print << "EOF";

This script was used to evaluate the quality of the raw data\n
Usage:\n
perl runclean-01.pl     --inputdir      <inputdir contains metrix files, suffix of file is "metrix"> 
                        --outfile      <formatted output file>
EOF
}
