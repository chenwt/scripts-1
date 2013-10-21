#!/usr/bin/perl -w

open(INF, "omim_result.txt");
open(OUTF, ">gene_AML_OMIM.txt");
while(<INF>) {
	if(m/^\W\d+/) {
	    m/\w+$/;
	    print OUTF "$&\n";
	}
}

close(INF);
