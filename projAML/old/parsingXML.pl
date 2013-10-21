#! usr/bin/perl -w

use strict;
use XML::Simple;
use Data::Dumper;
use FileHandle;

my($outfile) ="clinic.txt";
my $xml =new XML::Simple;

my $data = $xml->XMLin("GSE1159_family.xml");

my($out) = new FileHandle ">$outfile";
print $out Dumper($data);
close $out;

# my $count = 1;
# foreach($data)
# {
# 	if($count <= 10){print;}
# 	$count++;
# }