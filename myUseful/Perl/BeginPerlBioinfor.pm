
################################################################################
####################### sub match_positions ##########################
######## this subroutine is to find the position of a given seq regexpr
#input: the regular expression of a dna seq; return: an array of positions

sub match_positions {
        my($regexp, $dna) = @_;

        use strict;
        use warnings;
        use BeginPerlBioinfor;

        my @positions = ();

        while ($$dna =~ /$$regexp/ig){
                push( @positions, (pos($$dna) - length($&) + 1) );
        }
        return @positions;
}

################---------------------------------------------------------------

######################################################################################
######################### sub parseREBASE ################################
# this subroutine is to parse a rebase file
# input: rebase filename, return: a hash

sub parseREBASE {
        my($rebasefile) = @_;

        use strict;
        use warnings;
        use BeginPerlBioinfor;

        my @rebasefile = ( );
        my $name;
        my $site;
        my $regexp;
        my %rebase_hash = ( );
        my $count_header;
        @rebasefile = Get_file_data($rebasefile);

        foreach (@rebasefile){
                ++$count_header;
                ($count_header < 10) and next;

                /^\s*$/ and next;

                my @fields = split (" ", $_);

                $name = shift @fields;
                $site = pop @fields;
                $regexp = IUB_to_regexp($site);
                $rebase_hash{ $name} = "$site $regexp";
        }
        return %rebase_hash;
}
#########-----------------------------------------------------------------------------

#####################################################################################
##################### sub IUB_to_regexp ##################################
## this subroutine is to change the icu code seq into regular expression
# input:$iub seq (ATGCNM....), output:$regexp

sub IUB_to_regexp {

        my($iub_seq) = @_;

        my $reg = '';
        my $regexp = '';
        my %iub_hash = (
        A => 'A', C => 'C',
        G => 'G', T => 'T',
        R => '[GA]', Y => '[CT]',
        M => '[AC]', K => '[GT]',
        S => '[GC]', W => '[AT]',
        B => '[CGT]', D => '[AGT]',
        H => '[ACT]', V => '[ACG]',
        N => '[ACGT]',
        );

        $iub_seq =~ s/\^//g;
	$iub_seq = uc $iub_seq;
        for (my $i = 0; $i < length($iub_seq); ++$i){
                $reg = substr($iub_seq, $i, 1);
                $regexp .= $iub_hash{$reg};
        }

        return $regexp;
}

###########--------------------------------------------------------------------------




##############################################################################
 
############## sub Rev_com_dna ################################
# this subroutine is to generate the reverse complement of a dna seq
# input: $dna, return: $revcom_dna

sub Rev_com_dna {
        my($dna) = @_;
        my $revcom_dna = reverse $dna;
        $revcom_dna =~ tr/ATGCatgc/TACGtacg/;

#### the following is another method....thought much more complex...
#       my $base
#       my $revcom_dna = "";
#       my $length = length($dna);
#       for( my $i =0; $i < $length; ++$i){
#               $base = substr($dna, $i,1);
#               if(    $base =~ /A/ig){
#                       $revcom_dna .= 'T';
#               }elsif($base =~ /[ TU]/ig){
#                       $revcom_dna .= 'A';
#               }elsif($base =~ /G/ig){
#                       $revcom_dna .= 'C';
#               }elsif($base =~ /C/ig){
#                       $revcom_dna .= 'G';
#               }else{
#                       print "Cannot recognize base: $base!\n";
#               }
#       }
        
        return $revcom_dna;

}
#########----------------------------------------------------------------------


############################################################
###################### sub countATGC ################
sub countATGC {
        my ($infile) = @_;
        chomp $infile;
        open (INFILE, $infile);
        my (@dna) = <INFILE>;
        close INFILE;
        my ($dna) = join ("", @dna);
        $dna =~ s/\s//g;
        my $count_A = ($dna =~ tr/Aa//);
        my $count_T = ($dna =~ tr/Tt//);
        my $count_G = ($dna =~ tr/Gg//);
        my $count_C = ($dna =~ tr/Cc//);
        return ($count_A, $count_T, $count_G, $count_C);
}

########----------------------------------------------------

####################################################################################
######################## sub Condon2aa ################################
# this subroutine is to translate one condon to it's amino acid,input a 3 nucleotides, return a amino acid

sub Condon2aa {
        my($condon) = @_;

        $condon = uc $condon;

        my(%genetic_code) = (
        'TCA' => 'S',   'TCC' => 'S',
        'TCG' => 'S',   'TCT' => 'S',
        'TTC' => 'F',   'TTT' => 'F',
        'TTA' => 'L',   'TTG' => 'L',
        'TAC' => 'Y',   'TAT' => 'Y',
        'TAA' => '_',   'TAG' => '_',
        'TGC' => 'C',   'TGT' => 'C',
        'TGA' => '_',   'TGG' => 'W',
        'CTA' => 'L',   'CTC' => 'L',
        'CTG' => 'L',   'CTT' => 'L',
        'CCA' => 'P',   'CCC' => 'P',
        'CCG' => 'P',   'CCT' => 'P',
        'CAC' => 'H',   'CAT' => 'H',
        'CAA' => 'Q',   'CAG' => 'Q',
        'CGA' => 'R',   'CGC' => 'R',
        'CGG' => 'R',   'CGT' => 'R',
        'ATA' => 'I',   'ATC' => 'I',
        'ATT' => 'I',   'ATG' => 'M',
        'ACA' => 'T',   'ACC' => 'T',
        'ACG' => 'T',   'ACT' => 'T',
        'AAC' => 'N',   'AAT' => 'N',
        'AAA' => 'K',   'AAG' => 'K',
        'AGC' => 'S',   'AGT' => 'S',
        'AGA' => 'R',   'AGG' => 'R',
        'GTA' => 'V',   'GTC' => 'V',
        'GTG' => 'V',   'GTT' => 'V',
        'GCA' => 'A',   'GCC' => 'A',
        'GCG' => 'A',   'GCT' => 'A',
        'GAC' => 'D',   'GAT' => 'D',
        'GAA' => 'E',   'GAG' => 'E',
        'GGA' => 'G',   'GGC' => 'G',
        'GGG' => 'G',   'GGT' => 'G',
        );
        if( exists $genetic_code{$condon} ){
                return $genetic_code{ $condon};
        }else {
                print STDERR "Bad condon\' $condon \'!\n";
                exit;
        }
}


#############-----------------------------------------------------
#################################################################################
################sub dna2peptide #############################
# this subroutine is to translate dna seq into aa chain,
# input a value of $dna, return a $protein

sub dna2peptide {
        use strict;
        use warnings;

        my($dna) = @_;
        my($protein) = '';
        for (my $i = 0; $i < (length($dna) - 2); $i += 3) {
                $protein .= Condon2aa(substr($dna, $i, 3));
        }
        return $protein;
}
#########-----------------------------------------------------------################


############################################################################
####################### sub fasta2seq ############################
# this subroutine is to change fasta data format into dna sequence, wiht a given fastafile, and return a unformatted sequence
# input a array of @fasta, return a value of $seq

sub fasta2seq {
        my(@fasta) = @_;

        use strict;
        use warnings;
        my $seq = '';

#       print "in sub:\n",@fasta,"\n";

        foreach my $line (@fasta){
#               print $line, "\n";
                if(    $line =~ /^\s*$/ ){
                        next;
                }elsif( $line =~ /^\s*#/){
                        next;
                }elsif($line =~ /^>/ ){
                        next;
                }else{
                        $seq .= $line;
                }
        }
        $seq =~ s/\s//g;
        return $seq;
}
#######--------------------------------------------##########################


###############################################################################
######################### sub Get_file_data ########################
# this subroutine is to get fasta file
# input a $filename, return an array of @filedata

sub Get_file_data {
        my($filename) = @_;

        use strict;
        use warnings;

        my @filedata = ( );
        unless( open( GET_FILE_DATA, $filename)) {
                print STDERR "Cannot open file \'$filename\'\n";
                exit;
        }
        @filedata = <GET_FILE_DATA>;
        close GET_FILE_DATA;

        return @filedata;
}
#########-------------------------------------------#################################

#######################################################################################
#######################sub Print_format_seq ############################
###this subroutine is to format dna seq as a printable formation, return an printable seq set with a seq and line length as input
## passing by reference, input $seq, and $line_length, print out formated seq, without any return

sub Print_format_seq {
        my($seq, $line_length) = @_;
        my $i;
        for( $i = 0; $i < length($$seq); $i += $$line_length){
                print substr($$seq, $i, $$line_length), "\n";
        }
}
##########------------------------------------------------------------------------
1;
