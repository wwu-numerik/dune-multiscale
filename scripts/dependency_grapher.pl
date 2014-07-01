#!/usr/bin/perl
#
# Please comment any substantive changes with your email address and the date,
# and add a short unique suffix to the version string.

our $VERSION = '2.0';


=head1 NAME

cppdeps - Perform static analysis of C++ code

=head1 SYNOPSIS

cppdeps [options] [file ...]

=head1 DESCRIPTION

Prints a table listing the interdependencies of all C++ source and header
files found in the search directory. Dependency is determined by #include
directives found in each file. You can specify a dependency without using
#include, by using '#pragma depends "filename"' instead. Other preprocessor
macros (e.g. #ifdef) are not honoured. Perl's @INC array is used when trying
to locate named files.

If files are listed on the command line, then those files will be inspected,
otherwise a directory (usually the current directory) will be searched
recursively for files with names ending in standard C++ extensions.

The first stage lists each component (i.e. source/header pair) in alphanumeric
order with its direct dependencies indented below. The second stage lists the
"levels" and the components which are members of each one. Direct circular
dependency is indicated by an asterisk ('*'). Indirect circular dependency is
not supported and behaviour in this case is undefined. Finally, it prints a
coupling metric computed from the amount of interdependency found. Smaller is
better for this metric; anything < 1.5 is good, anything > 2.0 is getting bad.

The members of level 1 are those components which can be tested independently
of any other component. The members of level 2 are those components which can
be tested independently of any component except one or more in level 1. The
members of level N are those components which can be tested independently of
any component in levels > N, but may depend on components in levels < N.

=head1 OPTIONS

=over 8

=item B<-d> | B<--dir> I<directory>

Specifies the directory through which to search (recursively). In the absence
of this option, the default is the current working directory.

=item B<-h> | B<--help>

Prints this manual page (synopsis and options only) and exits.

=item B<-H> | B<--morehelp>

Prints this manual page (including detailed description) and exits.

=item B<-I> | B<--include> I<directory>

Adds a directory name to Perl's @INC array, which is used when searching for
source files included from other source files.

=item B<-i> | B<--include-extract> I<perlcode>

Specifies a snippet of Perl code which will extract the name of an include
file from the $_ variable and return it. The $_ variable will contain a line
of code taken from the input file, which may or may not include a '#include'
directive. If it does not, the code snippet should return zero.

Be sure to quote the code to avoid expansion by your shell.

=item B<-l> | B<--lang> I<language>

Specifies the programming language of the files being processed. The effects
of this specification can be conditionally overridden by the
B<include-extract>, B<--regex> and B<--wildcard> options.

The supported language specifiers are B<c++> and B<java>. The default is
B<c++>.

=item B<-p> | B<--png> I<pngfile>

In addition to printing a textual analysis to stdout, uses GraphViz to create
a PNG image file with the specified name, showing the interdependencies
diagramatically.

=item B<-x> | B<--imageopt> I<option>

Specifies an option that will influence the generation of the PNG diagram:

=over 8

=item B<level>

Causes components to be arranged according to their levels.

=item B<dir>

Causes components to be grouped according to their parent directory.

=back

=item B<-r> | B<--regex> I<pattern>

Specifies a Perl regular expression against which to match filenames while
searching. The default is '\.(h|H|hh|HH|hpp|hxx|c|C|cc|CC|cpp|cxx)$'. This
option is ignored if filenames are supplied on the command line.

Be sure to quote the pattern to avoid expansion by your shell.

=item B<-v> | B<--version>

Prints the version number and exits.

=item B<-w> | B<--wildcard> I<pattern>

Specifies a simple pattern string against which to match filenames while
searching. Only the asterisk ('*') is recognised as a special character.
This option is ignored in the presence of the B<--regex> option, or if
filenames are supplied on the command line.

Be sure to quote the pattern to avoid expansion by your shell.

=item B<-W> | b<--warnings>

Causes verbose warnings to be printed to stderr. These are limited to
reporting any processed or included file which could not be found or could
not be read from.

=back

=head1 TODO

- Add an option to list indirect as well as direct dependencies.

- Support indirect circular dependencies.

=head1 BIBLIOGRAPHY

I<Large Scale C++ Software Design>, John Lakos, 1996,
Addison-Wesley Publishing.

=head1 AUTHOR

The Oktalist E<lt>I<mat@oktalist.com>E<gt>

=head1 COPYING

This program is free software. You are free to modify and/or redistribute
it under the same terms as Perl itself.

=cut


use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

GetOptions(\our %Opts, qw(
	dir|d=s
	help|h
	morehelp|H
	include|I=s@
	include-extract|i=s
	lang|l=s
	png|p=s
	imageopt|x=s@
	regex|r=s
	version|v
	warnings|W
	wildcard|w=s
)) or pod2usage(-verbose => 1, -exitval => 1);

$Opts{morehelp} and pod2usage(-verbose => 2, -exitval => 0);
$Opts{help} and pod2usage(-verbose => 1, -exitval => 0);
$Opts{version} and (print("$0 version $VERSION\n"), exit 0);

$Opts{dir} ||= '.';
unshift @INC, $Opts{dir}, @{ $Opts{include} || [] };

$Opts{lang} ||= 'c++';

defined $Opts{wildcard} and $Opts{regex} ||= wildcard2regex($Opts{wildcard});

$Opts{regex} ||= {
	'c++' => '\.(h|H|hh|HH|hpp|hxx|c|C|cc|CC|cpp|cxx)$',
	'java' => '\.java$'
}->{ $Opts{lang} };

our $MatchSub = {
	'c++' => sub {
		my $line = shift;
		if ($line =~ m'^\s*#\s*(include|pragma\s+depends)\s+["<](.*?)[">]')
			{ return $2 }
		else
			{ return 0 }
	}, 'java' => sub {
		my $line = shift;
		if ($line =~ m'^\s*import\s+([\w\.]+)')
			{ my $f = $1; $f =~ tr{\.}{/}; return "$f.java" }
		else
			{ return 0 }
	}
}->{ $Opts{lang} };

defined $Opts{'include-extract'}
	and $MatchSub = sub { $_ = shift; eval $Opts{'include-extract'};
			defined $@ and (print STDERR "error in include-extract: $@", exit 1) };

if (@ARGV)
{
	CppDeps::processFile($_, $Opts{dir}, qr"$Opts{regex}") foreach @ARGV;
}
else
{
	recurseDir($Opts{dir}, \&CppDeps::processFile, qr"$Opts{regex}");
}

CppDeps::printDeps();

my @imageopts = split ',', join(',', @{ $Opts{imageopt} || [] });
$Opts{png} and CppDeps::writePng($Opts{png}, @imageopts);


# sub recurseDir($dirname, $callbackfunc)
#
# Calls the coderef $callbackfunc->($filename) for each file $filename in the
# directory $dirname, and calls recurseDir($subdirname, $callbackfunc) for each
# subdirectory in the directory $dirname.
#
sub recurseDir
{
	my $dirname = shift;
	my $callbackfunc = shift;

	if (opendir my $dir, $dirname)
	{
		while (defined(my $direlemname = readdir $dir))
		{
			next if $direlemname =~ m'^\.';

			if (-d "$dirname/$direlemname")
			{
				recurseDir("$dirname/$direlemname", $callbackfunc, @_);
			}
			else
			{
				$callbackfunc->($direlemname, $dirname, @_);
			}
		}
	}
	else
	{
		warn "Can't open directory '$dirname'";
	}
}

# sub wildcard2regex
#
# Converts a simple asterisk wildcard into a regular expression string.
#
sub wildcard2regex
{
	my $wildcard = quotemeta shift;

	$wildcard =~ s{\\\*}{.*};
	$wildcard = "^$wildcard\$";

	return $wildcard
}


package CppDeps;

use strict;
use warnings;

our %FileDeps;
# The component dependency bit-matrix. A true value in $FileDeps{A}->{B}
# indicates that component A depends on component B.

our %FileLevels;
# Hash associating components with levels. The numeric value $FileLevels{A}
# denotes the level of component A.

our @Levels;
# A list of lists containing the components contained within each level.
# $Levels[N] is an arrayref containing the component names in level N.

our %FilesChecked;
# Hash for remembering which files have been parsed so we don't parse any file
# twice. A true value for $FilesChecked{A} indicates file A has been parsed.

our %CircularDepWarned;
# A hash for keeping track of which circular dependencies have been warned
# about, so we don't warn about the same one twice.

# sub CppDeps::processFile($filename [ , $regex [ , $currentDir ] ])
#
# Searches @INC for a file named $filename and calls findDeps($fullPathToFile)
# on the first one found. Does nothing if $filename does not end in a standard
# C or C++ file extension. $regex is an optional argument; a file will be
# ignored if its name does not match. $currentDir is an optional argument; it
# will be searched before @INC. Returns the full path of the file it found.
# If no file was found, returns the empty string.
#
sub processFile
{
	my $filename = shift;
	my $currentDir = shift;
	my $regex = shift;

	if (defined $regex) { return '' unless $filename =~ $regex; }

	foreach my $incdir ($currentDir, @INC)
	{
		next unless defined $incdir;

		my $fqfilename = "$incdir/$filename";
		if (-f $fqfilename and -r $fqfilename)
		{
			defined findDeps($fqfilename)
				and return $fqfilename;
		}
	}

	if ($::Opts{warnings})
	{
		my $dirs = defined($currentDir) ? "'$currentDir' or \@INC" : '@INC';
		print STDERR "warning: no readable file '$filename' in $dirs\n";
	}
	return '';
}

# sub CppDeps::findDeps($filename)
#
# Reads the file $filename, and for each "#include" or "#pragma depends"
# directive found, sets the dependency flag within the %FileDeps matrix and
# calls processFile($includedFilename, $currentDir). Does nothing if $filename
# has already been read.
#
sub findDeps
{
	my $filename = shift;

	$filename = normalizePath($filename);
	my $basename = extractBasename($filename);
	my $pathname = extractPathname($basename);

	if (exists $FilesChecked{$filename})
	{
		return 0;
	}
	$FilesChecked{$filename} = 1;

	if (not exists $FileDeps{$basename})
	{
		$FileDeps{$basename} = {};
	}

	open my $fh, "< $filename" or return undef;

	while (defined(my $line = <$fh>))
	{
		chomp $line;
		my $incfilename = $::MatchSub->($line);

		if ($incfilename)
		{
			$incfilename = processFile($incfilename, $pathname);

			if ($incfilename)
			{
				$incfilename = normalizePath($incfilename);
				my $incbasename = extractBasename($incfilename);

				if ($basename ne $incbasename)
				{
					$FileDeps{$basename}{$incbasename} = 1;
				}
			}
		}
	}
	return 1;
}

# sub CppDeps::printDeps()
#
# Pretty-prints the %FileDeps matrix, then populates the %FileLevels hash and
# the @Levels array and then pretty-prints that too. Finally, computes and
# prints a coupling metric based on the amount of interdependence found.
#
sub printDeps
{
	%FileDeps or (print("No files found. Use the -h option for help.\n"), return);

	foreach my $key (sort keys %FileDeps)
	{
		my @deps = sort keys %{ $FileDeps{$key} };

		print "$key:\n";
		print "\t$_\n" foreach @deps;
		print "\n";

		if (0 == @deps)
		{
			$FileLevels{$key} = 1;
		}
	}

	for (0..keys %FileDeps)
	{
		foreach my $key (keys %FileDeps)
		{
			my $basename = extractBasename($key);

			my $level = exists($FileLevels{$key})
				    ? $FileLevels{$key}
				    : 0;

			foreach my $dep (keys %{ $FileDeps{$key} })
			{
				my $depbasename = extractBasename($dep);

				my $deplevel = exists($FileLevels{$dep})
						 ? $FileLevels{$dep}
						 : -1;

				if ($level < $deplevel + 1
				    and $basename ne $depbasename)
				{
					if (exists $FileDeps{$dep}->{$key})
					{
						$FileLevels{$key} = $deplevel;

						warnCircularDep($key, $dep);
					}
					else
					{
						$FileLevels{$key} = $deplevel + 1;
					}
				}
			}
		}
	}

	foreach my $key (keys %FileLevels)
	{
		push @{ $Levels[ $FileLevels{$key} ] }, $key;
	}

	foreach my $level (0..@Levels)
	{
		next unless defined $Levels[$level];

		# follow indirect dependencies
		foreach my $key (@{ $Levels[$level] })
		{
			foreach my $dep (keys %{ $FileDeps{$key} })
			{
				$FileDeps{$key} = {
					%{ $FileDeps{$key} },
					%{ $FileDeps{$dep} }
				};
			}
		}

		my @deps = sort @{ $Levels[$level] };

		print "LEVEL $level:\n";

		foreach my $dep (@deps)
		{
			if (exists $CircularDepWarned{$dep})
			{
				print "\t$dep *\n";
			}
			else
			{
				print "\t$dep\n";
			}
		}

		print "\n" if @deps;
	}

	if (0 == @Levels)
	{
		print "# Dependency graph is not levelizable, as it has no leaves.\n\n";
	}

	my $filecount = 0;
	my $depcount = 0;
	foreach my $key (keys %FileDeps)
	{
		$filecount++;
		$depcount += keys %{ $FileDeps{$key} };
	}
	my $nccd = $depcount / ($filecount * log $filecount);
	printf "COUPLING_METRIC = %.4f\n", $nccd;
}

# sub CppDeps::writePng($filename)
#
# Uses GraphViz to draw and save a PNG diagram of the dependency graph.
#
sub writePng
{
	my $filename = shift;
	my %opts;
	$opts{$_} = 1 foreach @_;

	require GraphViz;

	my $graph = new GraphViz (
		directed => 1,
		bgcolor => 'white',
		name => 'cppdeps',
		concentrate => 1,
		rankdir => 1,
		node => {
			shape => 'box'
		},
		edge => {
			style => 'bold',
			arrowsize => 2
		}
	);

	foreach my $level (0..(@Levels - 1))
	{
		foreach my $key (@{ $Levels[$level] })
		{
			my %nodeopts;
			$opts{level} and $nodeopts{rank} = $level;
			$opts{dir} and $nodeopts{cluster} = extractPathname($key);
			$opts{dir} and $nodeopts{label} = extractFilename($key);

			$graph->add_node($key, %nodeopts);
		}
	}

	foreach my $key (keys %FileDeps)
	{
		foreach my $dep (keys %{ $FileDeps{$key} })
		{
			my $draw = 1;
			foreach my $dep2 (keys %{ $FileDeps{$key} })
			{
		# remove direct dependency wherever an indirect dependency exists
				if (exists $FileDeps{$dep2}{$dep}) { $draw = 0; last }
			}

			$draw and $graph->add_edge($key => $dep);
		}
	}

	$graph->as_png($filename);
}

# sub CppDeps::warnCircularDep($filename1, $filename2)
#
# Sets a bit in the %CircularDepWarned matrix, identifying that a circular
# dependency exists between $filename1 and $filename2.
#
sub warnCircularDep
{
	my $file1 = shift;
	my $file2 = shift;

	$CircularDepWarned{$file1}->{$file2} = 1;
	$CircularDepWarned{$file2}->{$file1} = 1;
}

# sub CppDeps::normalizePath($path)
#
# Returns $path with Windows filesystem "\" replaced with Unix style "/" and any
# number of leading "./" substrings stripped from the front.
#
sub normalizePath
{
	my $path = shift;

	$path =~ s{\\}{/}g;
	$path =~ s{\G\./}{}g;

	return $path;
}

# sub CppDeps::extractBasename($filename)
#
# Returns $filename with any file extension removed from the end.
#
sub extractBasename
{
	my $file = shift;

	$file =~ s{\.[^\./]+$}{};

	return $file;
}

# sub CppDeps::extractPathname($filename)
#
# Returns $filename with everything beyond the final "/" removed, leaving only
# the directory in which $filename is located.
#
sub extractPathname
{
	my $file = shift;

	$file =~ s{/[^/*]+$}{};

	return $file;
}

# sub CppDeps::extractFilename($filename)
#
# Returns $filename with everything before the final "/" removed, leaving only
# the name of the file relative to the directory in which it is located.
#
sub extractFilename
{
	my $file = shift;

	$file =~ s{^.*/}{};

	return $file;
}

