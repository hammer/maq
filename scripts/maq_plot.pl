#!/usr/bin/perl -w

# Author: lh3

# This script will call "gnuplot" almost everywhere.

use strict;
use warnings;
use Getopt::Std;

&usage if (@ARGV < 1);

my $version = "0.1.1";
my $command = shift(@ARGV);
my %func = (abpair=>\&abpair, depth=>\&depth_plot, mixed=>\&mixed);
die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}}();
exit(0);

sub depth_plot
{
	my %opts;
	getopts('c:D:ax:', \%opts);
	die qq(Usage: maq_plot.pl depth [-c chr] [-D 100] [-a] [-x 2] <prefix> <cns2win-output>\n) if (@ARGV < 2); 
	my $chr = $opts{c};
	my $max_depth = $opts{D}? $opts{D} : 100;
	my $xtics = $opts{x}? $opts{x} : 2;
	my $back = 5;
	my $prefix = shift(@ARGV);
	my ($fh1, $fh2, $fh3, $fh, $fhin);
	open($fhin, $ARGV[0]) || die;
	my ($mean, $var, $n) = (0, 0, 0);
	my $col = (defined $opts{a})? 5 : 4;
	while (<$fhin>) {
		my @t = split;
		$chr = $t[0] if (!$chr);
		if ($chr eq $t[0] && $t[7] >= 1e-6) { # then, not N
			++$n;
			$mean += $t[$col];
		}
	}
	$mean /= $n;
	open($fh1, ">$prefix.depth") || die;
	open($fh2, ">$prefix.outlier") || die;
	open($fh3, ">$prefix.gc") || die;
	seek($fhin, 0, 0);
	$n = 0;
	while (<$fhin>) {
		my @t = split;
		if ($chr eq $t[0]) {
			if ($t[$col] < 2 * $mean && $t[$col] > 0.5 * $mean) {
				if ($t[7] >= 1e-6) {
					++$n;
					$var += ($t[$col] - $mean) * ($t[$col] - $mean);
				}
				$fh = ($t[$col] < 1.5 * $mean)? $fh1 : $fh2;
			} else { $fh = $fh2; }
			my $depth = ($t[$col] >= $max_depth-$back)? $max_depth - $back * (1.0 - rand()) : $t[$col];
			print $fh "$t[1]\t$depth\t$t[7]\n";
			if (rand() < 0.1) {
				print {$fh3} "$depth\t$t[7]\t$t[1]\n";
			}
		}
	}
	$var = sqrt($var / $n);
	close($fh1); close($fh2); close($fh3); close($fhin);
	open($fh, "| gnuplot") || die;
	my $label = sprintf("Depth (%.2f +/- %.2f)", $mean, $var);
	print $fh qq(
set t po eps co so 24
set grid;
set ylab "Depth";
set xlab "Percent GC";
set out "$prefix.gc.eps";
plot "$prefix.gc" u 2:1 w d t "";
set t po eps co so lw 2 18
set xlab "Coordinate (Mbp)";
set ylab "Depth";
set y2ran [0:100];
set y2tic 10;
set xtics $xtics;
set yran [-20:$max_depth];
set size 5, 1;
set grid;
set out "$prefix.depth.eps";
set style line 3 lt 3 lw 3;
plot "$prefix.depth" u 1:2 t "$label" w d 1, "$prefix.outlier" u 1:2 t "Outlier" w p 4, "$prefix.gc" u 3:2 smo be axis x1y2 ls 3;
exit;
);
	close($fh);
}

sub abpair
{
	my %opts;
	getopts('m:M:c:', \%opts);
	die qq(Usage: maq_plot.pl abpair [-m 250] [-M 1000000] [-c chr] <prefix> <abpair-output>\n) if (@ARGV < 2);
	my $min_dist = (defined $opts{m})? $opts{m} : 250;
	my $max_dist = (defined $opts{M})? $opts{M} : 1000000;
	my $prefix = $ARGV[0];
	my $chr = $opts{c};
	my ($fh1, $fh2, $fh, $fhin);
	open($fhin, "sort +3 -nr $ARGV[1] |") || die;
	open($fh1, ">$prefix.badori") || die;
	open($fh2, ">$prefix.baddist") || die;
	my (%hash1, %hash2);
	while (<$fhin>) {
		my @t = split;
		shift(@t);
		$chr = $t[0] if (!$chr);
		my $dist = $t[4] - $t[3] + 1;
		if ($t[0] eq $chr && $dist >= $min_dist && $dist <= $max_dist) {
			next if ($hash1{$t[3]} || $hash2{$t[4]});
			$hash1{$t[3]} = $hash2{$t[4]} = 1;
			$fh = ($t[1] == 2)? $fh2 : $fh1;
			my $q = $t[2] + rand() - 0.5;
			$q = 0 if ($q < 0);
			print $fh $t[3]/2000000.0 + $t[4]/2000000.0, "\t$q\t$dist\n";
		}
	}
	close($fhin); close($fh1); close($fh2);
	# plot
	open($fh, "| gnuplot") || die("Fail to execute gnuplot");
	print $fh qq(
set t po eps co so;
set xlab "$chr coordinate (Mbp)";
set ylab "MapQ of abnormal pairs";
set out "$prefix.abpair.eps";
set size 5, 1;
set xtics 2;
set grid;
plot "$prefix.badori" t "bad orientation" w p 2, "$prefix.baddist" t "bad distance" w p 4;
exit;
);
	close($fh);
}

sub mixed
{
	my %opts = (m=>300, M=>1000000, D=>100);
	getopts('c:D:m:M:', \%opts);
	die qq(
Usage:   maq_plot.pl mixed [options] <prefix> <cns2win-output> <abpair-output>

Options: -c STR        target chromosome [the first one]
         -D INT        maximum depth [100]
         -m INT        minumum pair distance [300]
         -M INT        maximum pair distance [1000000]
\n) if (@ARGV < 3); 
	my $chr_cmd = $opts{c}? "-c $opts{c}" : '';
	system(qq{$0 depth $chr_cmd -D $opts{D} $ARGV[0] $ARGV[1]}) && die;
	system(qq{$0 abpair $chr_cmd -m $opts{m} -M $opts{M} $ARGV[0] $ARGV[2]}) && die;
	my $prefix = $ARGV[0];
	my $fh;
	open($fh, "| gnuplot") || die;
	print $fh qq(
set t po eps co so;
set xlab "Coordinate (Mbp)";
set ylab "MapQ of abnormal pairs";
set out "$prefix.eps";
set size 5, 1;
set y2ran [-20:$opts{D}];
set y2lab "Depth";
set y2tics 10;
set xtics 2;
set grid;
plot "$prefix.badori" t "bad orientation" w p 2, "$prefix.outlier" t "Outlier" axis x1y2 3, "$prefix.baddist" t "bad distance" w p 4;
exit;
);
	close($fh);
}

sub usage
{
	die qq(
Usage:   maq_plot.pl <command> <arguments>

Command: abpair     plot abnormal pairs along the chromosome
         depth      depth plot
         mixed      mixed plot
\n);
}
