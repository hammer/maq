#!/usr/bin/perl -w

# Author: lh3
# Version: 0.1.3

use strict;
use warnings;
use Sys::Hostname qw/hostname/;
use Cwd qw/getcwd/;
use Getopt::Std;

if (@ARGV < 1) {
  warn("Arch-OS: " . &get_cpu_sys . "\n");
  die("farm-run.pl [-R LSF_resources] [-c n_CPUs] [-q LSF_queue] <config> [<job_id>]\n");
}
my %opts;
getopts("R:c:q:", \%opts);
my $conf_file = $ARGV[0];
die("Fail to find config file '$conf_file'") unless (-f $conf_file);

use vars qw/%fr_config/;
require qq/$conf_file/;

$SIG{SEGV} = $SIG{TERM} = $SIG{KILL} = \&clean_disk; # this is a hook, but temporarily useless

&check_config(\%fr_config); # validate the config file
&fix_path(\%fr_config); # remove those "/nfs" in $PATH and $PERL5LIB

my $n_jobs = $fr_config{n_jobs};
$fr_config{-R} = (defined($opts{R}))? qq{-R"$opts{R}"} : '';
$fr_config{-q} = (defined($opts{q}))? qq{-q $opts{q}} : '';
$fr_config{-c} = (defined($opts{c}))? '%'.($opts{c}<$n_jobs?$opts{c}:$n_jobs) : '';
if (!$fr_config{-R} && $fr_config{LSF_resource}) {
	$fr_config{-R} = qq{-R"$fr_config{LSF_resource}"};
}
if (!$fr_config{-q} && $fr_config{LSF_queue}) {
	$fr_config{-q} = qq{-q $fr_config{LSF_queue}};
}

if (@ARGV == 1) {
	$conf_file =~ /([^\s\/]+)(\/)?$/;
	$_ = $1; s/\.pl$//;
	mkdir("$_.err"); mkdir("$_.out");
	my $command = "echo '$0 $conf_file".' ${LSB_JOBINDEX}\' '."| bsub -J$_"."\"[1-$n_jobs]$fr_config{-c}\" $fr_config{-R} $fr_config{-q} -o $_.out/\%I.out -e $_.err/\%I.err";
	print "$command\n"; # Just print it out. We can pipe to "| sh" to execute it.
} else {
	print STDERR "Hostname: ".hostname.", Arch-OS: ".&get_cpu_sys."\n";
	&farm_action(\%fr_config, $ARGV[1]);
}

exit 0;

#############################################

# check whether the config file is complete
sub check_config
{
	my $conf = shift;
	my @check = ('run_list', 'n_jobs', 'action', 'binary_path', 'script_path');
	for (@check) {
		die(qq([check_config] missing "$_")) unless (defined($conf->{$_}));
	}
	die("\$conf->{binary_path} should start with '/'") if ($conf->{binary_path} !~ /^\//);
	die("\$conf->{script_path} should start with '/'") if ($conf->{script_path} !~ /^\//);
}
# action happens here
sub farm_action
{
	my $conf = shift;
	my $job_id = shift;
	my $pwd = getcwd;
	my @list = &prepare_run_list($conf, $job_id);
	foreach my $dir (@list) {
		# $dir = "$pwd/$dir" if ($dir !~ /^\//); # this is not necessary.
		if (-d $dir) {
			$dir =~ /([^\s\/]+)(\/)?$/;
			my $rel_dir = $1; # relative directory name
			print STDERR "Processing $rel_dir...\n";
			chdir($dir);
			&{$conf->{action}}($rel_dir);
			chdir($pwd);
		} else {
			warn("[farm_run_action] '$dir' is not found.");
		}
	}
}
# the directories in this list will be processed on the computing node
sub prepare_run_list
{
	my $conf = shift;
	my $job_id = shift;
	my ($fh, $count, @list);
	$count = 0;
	open($fh, "$conf->{run_list}") || die("[gen_fam_list] fail to open $conf->{run_list}");
	while (<$fh>) {
		chomp;
		push(@list, $_) if ($count % $conf->{n_jobs} + 1 == $job_id);
		++$count;
	}
	close($fh);
	return @list;
}
# clean the mess if terminating signals are received
sub clean_disk
{ # no need to do this now.
}
# Add script_home and binary_home to the path.
sub fix_path
{
	my $conf = shift;
	my (@t, @s);
	my $dir = &get_cpu_sys;
	@t = ($ENV{PATH})? split(":", $ENV{PATH}) : ();
	@s = ("$conf->{binary_path}/$dir", "$conf->{script_path}", "$conf->{binary_path}");
	$ENV{PATH} = join(":", @s, @t);
	@t = ($ENV{PERL5LIB})? split(":", $ENV{PERL5LIB}) : ();
	@s = ("$conf->{script_path}");
	$ENV{PERL5LIB} = join(":", @s, @t);
}
# an improve version of mkdir(). It can make "dir1/dir2/dir3" even if dir1 does not exist.
sub gendir
{
	my $path = shift;
	return unless($path);
	$_ = $path;
	my $tmp = (/^\//)? '/' : '';
	s/([^\/\n]+\/)/$tmp.=$1,mkdir($tmp),$1/eg;
	return $path;
}
# get some system information
sub get_cpu_sys
{
	my $dir = `uname -p 2>/dev/null`;
	$dir = `uname -m 2>/dev/null` if (!$dir || $dir =~ "unknown");
	$dir .= '-'.`uname -s`;
	$dir = lc($dir);
	$dir =~ s/\s//g;
	return $dir;
}
