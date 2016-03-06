#!/usr/bin/perl -w

# Start pod
=pod

=head1 NAME

Nearest Neighbor dG calculator

=head1 AUTHOR

Ciro Ramirez Suastegui and Jose Damian Martinez Reyes.
From: Bruno Contreras-Moreira.

=head1 VERSION

	2.0

=head1 Created

<February 16, 2016 -15:12- cramirez>

=head1 USAGE

program_name [-h] file

=head1 DESCRIPTION

This program calculate the stability of DNA duplex
based on the Gibbs energy in order to identify
candidate regions to promoter.

=head1 OPTIONS OR PARAMETERS

=over 4

=item B<-h>

Print help

=back

=head1 REQUIREMENTS

Basic biochemistry. See: Kanhere & Bansal, (2005) algorithm;
SantaLucia, (1998) and Breslauer et al., (1986).

=head1 INPUT FORMAT

.bin

=head1 OUTPUT FORMAT

Plane text

=head1 LANGUAGE

English

=head1 SEE ALSO

http://eead-csic-compbio.github.io/bioinformatica_estructural/node25.html

=cut

# prog1.1 
# Bruno Contreras-Moreira
# Nearest Neighbor dG calculator

use strict;
use Getopt::Long; # Let us take data from the line
my %opts = (); # Hash declaration
# Add keys with their values
GetOptions (\%opts, 'h'); # be sure the field exist
if(($opts{'h'})){
   &PrintHelp();
}

# global variables
my $T           = 37; # temperature(C)
my $windowL     = 15;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/kmol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H',   0, 'S',-1.4 } );
my %seqs; # Sequences
my @names; #Sequence name
my $g_deltaG;

my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";

print "# parameters: Temperature=$T gC Window=$windowL\n\n";

open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
print("Reading sequences...\n");
print("Calculating dGs...\n");
while(<SEQ>)
{
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		my ($name,$seq) = ($1,$2);
		#printf("sequence %s (%d nts)\n",$name,length($seq));
		
		# Call the subroutine and save the parameters needed
		$g_deltaG = duplex_deltaG($seq, $T);
		$seqs{$name}{'seq'} = $seq; # Save string in a global hash %seqs
		$seqs{$name}{'dG'} = $g_deltaG;
		printf("%s \n", $seqs{$name});
	}
}
close(SEQ);

#######################
##### Subroutines #####
#######################

# calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
# parameters: 1) DNA sequence string; 2) Celsius temperature
# returns; 1) free energy scalar
# uses global hash %NNparams
sub duplex_deltaG
{
   	my ($seq,$tCelsius) = @_;
	
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);
	my @sequence = split(//,uc($seq));
	my $tK = 273.15 + $tCelsius;
	my $n;
	
	# add dG for overlapping dinucleotides
	for($n=0;$n<$#sequence;$n++)
	{
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.
				complement($sequence[$n].$sequence[$n+1]);
			
			if(!defined($NNparams{$DNAstep})) # Ask if there's the string
			{
				$DNAstep = reverse($DNAstep); # If no, it reverses
			}
			
			$dG = ((1000*$NNparams{$DNAstep}{'H'})-
					($tK*$NNparams{$DNAstep}{'S'}))
					/ 1000 ; # Calculating the dG value
			
			$total_dG += $dG;
	}
	
	# add correction for helix initiation
	$nt = $sequence[0]; # first pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))
					/ 1000;
	
	$nt = $sequence[$#sequence]; # last pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) }
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))
					/ 1000;
					
	# please complete for symmetry correction
	for(my $n=0;$n<length($seq)-15;$n++)
	{
		if (substr($seq, $n, 7) eq complement(reverse(substr($seq, $n+8, 7))))
		{
			$total_dG += ((1000*$NNparams{'sym'}{'H'})-
				($tK*$NNparams{'sym'}{'S'}))
				/ 1000;
		}
	}

	return $total_dG;
}
print("\nRESULTS\n--------------------\nName\tInit.\tdG\n");
foreach my $i (keys(%seqs))
{
	print($i,"\t", reverse(substr($seqs{$i}{'seq'}, 0, 5)),"\t", $seqs{$i}{'dG'}, "\n");
}
print("--------------------\n");

############# Obtain the complement string #############
sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }

##################################################################
#### Despliega la ayuda en linea con opcion -h ###################
##################################################################
sub PrintHelp {
   system "pod2text -c $0 "; # Convert POD data to formatted ASCII text 
   exit();
}
