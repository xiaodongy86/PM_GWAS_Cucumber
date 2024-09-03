#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#######################################################################################

my ($infile,$outfile);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($infile and $outfile);
##----------------------------------------------------------------------------------------------

# Main Body
##----------------------------------------------------------------------------------------------

if ($infile=~/\.gz/) {
	open (IN,"zcat $infile|") or die $!;
}else{
	open (IN,"$infile") or die $!;
}
if ($outfile=~/\.gz$/) {
        open (OUT,"|gzip>$outfile") or die $!;
}else{
        open (OUT,">","$outfile") or die $!;
}

while (<IN>) {
	chomp;
	if (/^#/) {
		print OUT $_,"\n";
	}else{
		my @info=split /\t/,$_;
		if($info[6]=~/PASS/){
			print OUT $_,"\n";
		}
	}
}
close (OUT) ;
close (IN) ;


##----------------------------------------------------------------------------------------------
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
##----------------------------------------------------------------------------------------------
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  Options:
	-i	<infile>	infile,
	-o	<outfile>	outfile,

E.g:
	perl $0 -i infile -o outfile 

USAGE
	print $usage;
	exit;
}
