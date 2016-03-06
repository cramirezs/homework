#!/usr/bin/perl -w
=pod

=head1 NAME

Fasta files generator

=head1 AUTHORS

Martinez Reyes Jose Damian and Ciro Ramirez Suastegui.

=head1 VERSION

	2.0

=head1 Created

<February 17, 2016 -23:34- cramirez/josemr>

=head1 USAGE

program_name [-h] file.pdb file.fa|file.fasta

=head1 DESCRIPTION

This program take a .pdb file and
generate a fasta file along with
the help of a .fa or .fasta file
in order to obtain the header.

=head1 OPTIONS OR PARAMETERS

=over 4

=item B<-h>

Print help

=back

=head1 REQUIREMENTS

.fasta and .pdb files familiarization.

=head1 INPUT FORMAT

.pdb and .fa|.fasta

=head1 OUTPUT FORMAT

.fasta

=head1 LANGUAGE

English

=head1 SEE ALSO

http://eead-csic-compbio.github.io/bioinformatica_estructural/node32.html

=cut

use strict;
use warnings;
use File::Basename;  # permit obtain file names
use autodie; # die if problem reading or writing a file
use Getopt::Long; # Let us take data from the line
my %opts = (); # Hash declaration
# Add keys with their values
GetOptions (\%opts, 'h'); # be sure the field exist
if(($opts{'h'})){
   &PrintHelp();
}

my $seq1 = $ARGV[0] || die "# usage: $0 <PDB file>\n";
my $fa = $ARGV[1] || die "# usage: $1 <fasta file\n>";
my $header;
my $seq;
my $name;

my %hash=('ALA'=>'a', 'ARG'=>'r','ASN'=>'n', # We firs construct a
		  'ASP'=>'d', 'CYS'=>'c', 'GLN'=>'q', # hash containing the
		  'GLU'=>'e', 'GLY'=>'g', 'HIS'=>'h', # character equivalences
		  'ILE'=>'i', 'LEU'=>'l', 'LYS'=>'k',
		  'MET'=>'m', 'PHE'=>'f', 'PRO'=>'p',
		  'SER'=>'s', 'THR'=>'t', 'TRP'=>'w',
		  'TYR'=>'y','VAL'=>'v');

open(SEQ, $seq1) || die "# cannot open input $seq1 : $!\n";
my $rep = 0;
while(<SEQ>)
{
	if(/^ATOM.{13,13}(\w{3,3}).*/) # stract the aa's abbreviations
	{
		if ($rep ne $hash{$1}) # check if the abbrevations ain't repeted
		{
			$seq = $seq.$hash{$1};	# construct the sequence of aa's
									# checking the equivalent
			$rep = $hash{$1};	# saving the last aa in order to avoid
								# to repeat them
        }
	}
}
close(SEQ);
open(FA, $fa) || die "# cannot open input $fa : $!\n";
while (<FA>)
{
	if (/^>(.*)/) # stracting the header of the .fa file
	{
		$header = $1;
	}
$name = basename($fa)."sta"; # taking the domain name from
							  # the .fa file
}
close(FA);
# Writing the fasta file
open(FASTA, ">$name") || die "# cannot create the *.fasta file : $!\n";
print FASTA ">$header\n$seq\n";
close(FASTA);
print("Your fasta file is $name\n");


##################################################################
#### Despliega la ayuda en linea con opcion ah ###################
##################################################################
sub PrintHelp {
   system "pod2text -c $0 "; # Convert POD data to formatted ASCII text 
   exit();
}
