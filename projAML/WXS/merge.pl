#!/usr/bin/perl -w
#perl merge.pl -f1 file1 -f2 file2 -f3 file3 -f4 file4
%hash_a = ();
%hash_b = ();
%hash_c = ();
%hash_union = ();
open(FILE1,"PATABB-NoA.freq");
while(<FILE1>) {
  @current_line = split(/\s+/, $_);
  $key = $current_line[0] . "," . $current_line[1];
  @{$hash_a{$key}} = [];
  push (@{$hash_a{$key}},$current_line[2]);
  push (@{$hash_a{$key}},$current_line[3]);
  if(!(exists $hash_union{$key})) {
    ${$hash_union{$key}} = 0; 
  }
}
close(FILE1);

open(FILE2,"PATABB-TuA.freq");
while(<FILE2>) {
  @current_line = split(/\s+/, $_);
  $key = $current_line[0] . "," . $current_line[1];
  @{$hash_b{$key}} = [];
  push (@{$hash_b{$key}},$current_line[2]);
  push (@{$hash_b{$key}},$current_line[3]);
  if(!(exists $hash_union{$key})) {
    ${$hash_union{$key}} = 0; 
  }
}
close(FILE2);

open(FILE3,"PATABB-ReA.freq");
while(<FILE3>) {
  @current_line = split(/\s+/, $_);
  $key = $current_line[0] . "," . $current_line[1];
  @{$hash_c{$key}} = [];
  push (@{$hash_c{$key}},$current_line[2]);
  push (@{$hash_c{$key}},$current_line[3]);
  if(!(exists $hash_union{$key})) {
    ${$hash_union{$key}} = 0; 
  }
}
close(FILE3);

open(FILE4,">PATABB.freq");
foreach (sort keys %hash_union) {
  $key = $_;
  @column_12 = split(/,/,$key);
  $column_1 = $column_12[0];
  $column_2 = $column_12[1];

  if(!(exists $hash_a{$key})) {
    @{$hash_a{$key}} = [];
    push (@{$hash_a{$key}},0);
    push (@{$hash_a{$key}},0);
  }
  if(!(exists $hash_b{$key})) {
    @{$hash_b{$key}} = [];
    push (@{$hash_b{$key}},0);
    push (@{$hash_b{$key}},0);
  }
    if(!(exists $hash_c{$key})) {
    @{$hash_c{$key}} = [];
    push (@{$hash_c{$key}},0);
    push (@{$hash_c{$key}},0);
  }
  print FILE4 $column_1 . " " . $column_2 . "\t" . ${$hash_a{$key}}[1] . "\t" .  ${$hash_a{$key}}[2] . "\t" . ${$hash_b{$key}}[1] . "\t" . ${$hash_b{$key}}[2] . "\t" . ${$hash_c{$key}}[1] . "\t" . ${$hash_c{$key}}[2] . "\n";
}
close(FILE4);