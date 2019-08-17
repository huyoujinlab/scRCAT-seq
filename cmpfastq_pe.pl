#!/usr/bin/perl

##############################################################################
#
#                                 cmpfastq
#
# Concept: Stephen Newhouse (stephen.newhouse@kcl.ac.uk)
# Author:  David To (david.to@kcl.ac.uk)
# Written: 30th June 2010
# Edited: Simon Topp 1st Aug 2011 to make more memory efficient - note now only
#         works on files produced from the same lane and machine.
#
# DESCRIPTION:
# This script is designed to compare two fastq files.
# It produces 4 output files
# - 2 files that is a list of common sequences of both files
# - 2 files that are unique to each file
#
#
# Copyright 2010 NIHR Biomedical Research Centre for Mental Health
#                South London and Maudsley NHS Foundation Trust &
#                Institute of Psychiatry
#                Kings College London
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see:
# <http://www.gnu.org/licenses/>.
##############################################################################

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

## CONSTANTS ##
my $TRUE              = 1;
my $FALSE             = 0;
my $DEBUG             = $FALSE;
my $EXITSTATUS        = 0;

# Default umask
umask 027;

# Define variables

# Get the options passed at the command prompt
GetOptions();

##############################
#   M A I N   P R O G R A M  #
##############################

# Check to see we received two files in the arguments
if(scalar(@ARGV) != 2)
{
        print STDERR "Incorrect number of arguments\n";
        Usage();
        exit(1);
}

my $fail = $FALSE;

# Check to see if the files exist
foreach my $file (@ARGV)
{
        if(!-e $file)
        {
                print STDERR "File $file didn't exist\n";
                $fail = $TRUE;
        }
}

# If any of the files didn't exist, let's kill it
if($fail)
{
        exit(1);
}

# Read the file names in from the command line
my $file1 = shift(@ARGV);
my $file2 = shift(@ARGV);

my $date = `date`;
chomp($date);
print STDERR "BEGIN cmpfastq3 on $file1 $file2 at $date\n";

# Index the first file.
my %fastqIndex1 = %{IndexFastq($file1)};

# Compare the two files
CompareFastq($file1, $file2, \%fastqIndex1);

$date = `date`;
chomp($date);
print STDERR "END cmpfastq3 on $file1 $file2 at $date\n";

exit($EXITSTATUS);

# Subroutines
sub Usage
{
    my $base = basename($0);
    print "Usage: $base [dh] file1 file2\n";
    print "\td:\tDebug mode on (default off)\n";
    print "\th:\tPrint this usage\n";
}

sub GetOptions
{
    # Get the options passed at the command prompt
    my %options=();
    getopts("dh", \%options);

    if(defined($options{'d'}))
    {
        $DEBUG = $TRUE;
    }

    if(defined($options{'h'}))
    {
        Usage();
        exit($EXITSTATUS);
    }
}

sub IndexFastq
{
    my $file = shift;
    my %fastqIndex;

    open(IN, $file) or die("Could not open $file\n");
    my $pos = tell(IN);
    my $seqCounter = 1;
    my $lineCounter = 1;
    while(my $line = <IN>)
    {
        chomp($line);

        # Each block is going to be of 4 lines
        # Let's get the seq ID from the sequence name
        # Assuming both reads come from same machine and flowcell lane

        if($line =~ m/^@(.*?):(\d+):(\d+):(\d+):(\d+)\#.*/)  #CASAVA pre 1.7
        {
            #@HWUSI-EAS100R:6:73:941:1973#0/1
            #HWUSI-EAS100R  the unique instrument name
            #6     flowcell lane
            #73    tile number within the flowcell lane
            #941   'x'-coordinate of the cluster within the tile
            #1973  'y'-coordinate of the cluster within the tile
            #0     index number for a multiplexed sample (0 for no indexing)
            #/1   the member of a pair, /1 or /2 (paired-end or mate-pair reads only)
            #print STDOUT "Sequences file1 read:$seqCounter\r" if ($seqCounter%1000 ==0);
            $seqCounter++;
            $fastqIndex{$3}{$4}{$5} = $pos;
            # Skip the next 3 lines
            for(my $i=0; $i<3; $i++)
            {
                <IN>;
                $lineCounter++;
            }
        }

        elsif($line =~ m/^@(.*?):(\d+):(.*?):(\d+):(\d+):(\d+):(\d+)\s\d+:[YN]/) # CASAVA 1.8+
        {
            #@HWI-ST1025:64:D09CVACXX:8:1101:1133:1921 2:N:0:TGACCA
            #HWI-ST1025   Instrument
            #64           Run
            #D09CVACXX    Flowcell
            #8            Lane
            #1101         Tile
            #1133         X coord
            #1921         Y coord
            #2            Pair number [1/2]
            #N            Fails QC [Y/N]
            #0            No control bits [even numbers]
            #TGACCA       Index tag

            #print STDOUT "Sequences file1 read:$seqCounter\r" if ($seqCounter%1000 ==0);
            $seqCounter++;
            $fastqIndex{$5}{$6}{$7} = $pos;
            # Skip the next 3 lines
            for(my $i=0; $i<3; $i++)
            {
                <IN>;
                $lineCounter++;
            }
        }
        elsif($line =~ m/^\#/)
        {
            print STDERR "File: $file\[$lineCounter]: Skipping comment line: $line\n" if($DEBUG);
        }
        elsif($line =~ m/^$/)
        {
            print STDERR "File: $file\[$lineCounter]: Skipping empty line: $line\n" if($DEBUG);
        }
        else
        {
            print STDERR "File: $file\[$lineCounter]: Could not match the sequence ID from the name: $line\n" if($DEBUG);
        }
        $pos = tell(IN);
        $lineCounter++;
    }
    close(IN);

    return \%fastqIndex;
}
print STDOUT "\n";
sub CompareFastq
{
    my $file1          = shift;
    my $file2          = shift;
    my $fastqIndex1Ref = shift;
    my %fastqIndex1    = %{$fastqIndex1Ref};
    my %found1;
    my $seqCounter = 1;
    # We don't want to have to open/close file handles each time, so let's open them here
    open(F1COUT, ">$file1-common.out") or die("Could not write to file: $file1-common.out\n");
    open(F2COUT, ">$file2-common.out") or die("Could not write to file: $file2-common.out\n");
    open(F1UOUT, ">$file1-unique.out") or die("Could not write to file: $file1-unique.out\n");
    open(F2UOUT, ">$file2-unique.out") or die("Could not write to file: $file2-unique.out\n");

    open(F1IN, $file1) or die("Could not open $file1\n");
    open(F2IN, $file2) or die("Could not open $file2\n");
    while(my $line = <F2IN>)
    {
        chomp($line);

        # Skip empty lines or comments
        if($line =~ m/^$/g or $line =~ m/^\s*\#/)
        {
            next;
        }

        # Each block is going to be of 4 lines
        # Let's get the seq ID from the sequence name

        if($line =~ m/^@(.*?):(\d+):(\d+):(\d+):(\d+)\#.*/)  #CASAVA  V1.7 and previous

        {
            #print STDOUT "Sequences file2 read:$seqCounter\r" if ($seqCounter%1000 ==0);
            $seqCounter++;
            if(defined($fastqIndex1{$3}{$4}{$5}))
            {
                $found1{$3}{$4}{$5} = $TRUE;

                # Print out from file1
                seek(F1IN, $fastqIndex1{$3}{$4}{$5}, 0);
                for(my $i=0;$i<4;$i++)
                {
                    my $tmpLine = <F1IN>;
                    print F1COUT $tmpLine;
                }

                # Print out from file 2
                print F2COUT $line . "\n";
                for(my $i=0; $i<3; $i++)
                {
                    my $tmpLine = <F2IN>;
                    print F2COUT $tmpLine;
                }
            }
            else
            {
                # Print out from file 2
                print F2UOUT $line . "\n";
                for(my $i=0; $i<3; $i++)
                {
                    my $tmpLine = <F2IN>;
                    print F2UOUT $tmpLine;
                }
            }
        }

        elsif($line =~ m/^@(.*?):(\d+):(.*?):(\d+):(\d+):(\d+):(\d+)\s\d+:[YN]/) # CASAVA 1.8+
        {
            #print STDOUT "Sequences file2 read:$seqCounter\r" if ($seqCounter%1000 ==0);
            $seqCounter++;
            if(defined($fastqIndex1{$5}{$6}{$7}))
            {
                $found1{$5}{$6}{$7} = $TRUE;

                # Print out from file1
                seek(F1IN, $fastqIndex1{$5}{$6}{$7}, 0);
                for(my $i=0;$i<4;$i++)
                {
                    my $tmpLine = <F1IN>;
                    print F1COUT $tmpLine;
                }

                # Print out from file 2
                print F2COUT $line . "\n";
                for(my $i=0; $i<3; $i++)
                {
                    my $tmpLine = <F2IN>;
                    print F2COUT $tmpLine;
                }
            }
            else
            {
                # Print out from file 2
                print F2UOUT $line . "\n";
                for(my $i=0; $i<3; $i++)
                {
                    my $tmpLine = <F2IN>;
                    print F2UOUT $tmpLine;
                }
            }
        }
        else
        {
            print STDERR "Could not match the sequence ID from the name: $line\n";;
            next;
        }
    }
    close(F1COUT);
    close(F2COUT);
    close(F2UOUT);
    close(F2IN);
    print STDOUT "\n";
    # Now let's worry about the sequences that weren't common in file 1

        # File 1
    my $lastCounter = 0;
    foreach my $tile (keys %fastqIndex1)
    {
        foreach my $xcoord (keys %{$fastqIndex1{$tile}})
        {
            foreach my $ycoord (keys %{$fastqIndex1{$tile}{$xcoord}})
            {
                print STDOUT "Sequences last bit compared:$lastCounter\r" if ($lastCounter%1000 == 0);
                $lastCounter++;
                if(!defined($found1{$tile}{$xcoord}{$ycoord}))
                {
                    seek(F1IN, $fastqIndex1{$tile}{$xcoord}{$ycoord}, 0);
                    for(my $i=0;$i<4;$i++)
                    {
                        my $tmpLine = <F1IN>;
                        print F1UOUT $tmpLine;
                    }
                }
            }
        }
    }
    close(F1UOUT);
    close(F1IN)
}
