#!/usr/bin/perl

use Getopt::Long;
my $input_file = "input.txt";
my $output_file = "output.txt";
GetOptions ( "input=s" => \$input_file,
             "output=s" =>\$output_file);
open(SOURCE,"$input_file");
open(TEMP,">$output_file");
while(<SOURCE>) {
  chomp($_);
  @string = split(/-|s+/,$_);
    foreach(@string) {
       if($_ =~ /[a-z\W]/) {
       }elsif ($_ =~ /[A-Z]/) {
 	 if($_ =~ /[A-Z0-9]{2}/) {
 	    print TEMP "$_\n";
       	}
       }
   }
}
