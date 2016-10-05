#!/usr/bin/env perl

# --------------------------------------------------------------------
# %% Setting up %%
#

use strict;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use Try::Tiny::Retry;

use constant PROGRAM_NAME            => 'pMLST-1.4.pl';
use constant PROGRAM_NAME_LONG       => 'Calculate pMLST profile for a sequence or genome';
use constant VERSION                 => '1.4';

#Global variables
my $BLAST;
my $BLASTALL;
my $FORMATDB;
my $MLST_DB;
my ($Help, $Organism, $InFile, $dir);
my $IFormat = "fasta";
my $OFormat = "ST";
my %ARGV    = ('-p' => 'blastn', '-a' => '5' , '-F' => 'F');

#Getting global variables from the command line
&commandline_parsing();

my $REG_EXP   = "([_]?[a-zA-Z0-9]+[_]?)([-_])([a-zA-Z0-9]+)";              # Regular expression for matching allele names. $1 should match gene name, $2 the connector (- or _), $3 the allele number.

if (defined $Help || not defined $Organism || not defined $InFile) {
   print_help();
   exit;
}

#If there are not given a path to the database or BLAST the program assume that the files are located in the curet directury
if (not defined $BLAST) {
   $BLASTALL = "blastall";
   $FORMATDB = "formatdb";
}
if (not defined $MLST_DB) {
   $MLST_DB = "database";
}
if (not defined $dir) {
  $dir = ".";
}

# Making tmp directory for BLST output
my $tmp_dir = "$dir/tmp";
mkdir("$dir/tmp");

# --------------------------------------------------------------------
# %% Main Program %%
#

# Run BLAST and find best matching Alleles
my ($Seqs_mlst, $Seqs_input, @Blast_lines);

retry{
   $Seqs_mlst   = read_seqs(-file => $MLST_DB.'/'.$Organism.'.fsa', format => 'fasta');  
   $Seqs_input  = $InFile ne "" ? read_seqs(-file => $InFile, -format => $IFormat) :
                                    read_seqs(-fh => \*STDIN,   -format => $IFormat);

   @Blast_lines = get_blast_run($tmp_dir, $Organism, -d => $Seqs_input, -i => $Seqs_mlst, %ARGV);
}
catch{ die $_ };

## AVAILABLE ORGANISMS ##
# hash mapping the mlst profiles to organism names
my %mlstProfiles = ('incf' => 'IncF', 'inchi1' => 'IncHI1', 'inchi2' => 'IncHI2', 'inci1' => 'IncI1', 'incn' => 'IncN',
				'fia' => 'FIA', 'fib' => 'FIB', 'fii' => 'FII', 'fic' => 'FIC', 'fiik' => 'FIIK', 'fiis' => 'FIIS', 'fiiy' => 'FIIY',
				'hcm1043' => 'hcm1043', 'hcm1064' => 'hcm1064', 'hcm1099' => 'hcm1099', 'hcm1116' => 'hcm1116', 'hcm1178ac' => 'hcm1178ac', 'hcm1259' => 'hcm1259',
				'smr0018' => 'smr0018', 'smr0199' => 'smr0199',
				'repi1' => 'repI1', 'arda' => 'ardA', 'trba' => 'trbA', 'sogs' => 'sogS', 'pill' => 'pilL',
				'repn' => 'repN', 'traj' => 'traJ', 'kora' => 'korA');

#Declaring variables
my @RESULTS_AND_SETTINGS_ARRAY; #will contain the typing results some setting and the hash with the results for each gene
my @GENE_RESULTS_ARRAY; #will contain the typing results some setting and the hash with the results for each gene
my %GENE_ALIGN_HIT_HASH; #will contain the sequence alignment lines
my %GENE_ALIGN_HOMO_HASH; #will contain the sequence alignment homolog string
my %GENE_ALIGN_QUERY_HASH; #will contain the sequence alignment allele string
my %GENE_RESULTS_HASH; #will contain the results to be printed for each gene
my %MLST;
my %PERC_IDENT;
my %QUERY_LENGTH;
my %HSP_LENGTH;
my %GAPS;
my %Q_STRING;
my %HIT_STRING;
my %HOMO_STRING;
my %QUERY_START;
my %CONTIG_NAME;
my %HIT_STRAND;
my %HIT_START;
my %HIT_END;
my %HIT_LENGTH;

for my $blast_line (@Blast_lines) {  # Notice that a properly formatted mlst blastdb will have the gene name as the description and the allele as id
  chomp $blast_line;
  my @blast_elem = split ("\t",$blast_line);
  my $qid = $blast_elem[0];
  my $query_length = $blast_elem[1];
  my $hsp_length = $blast_elem[2];
  my $gaps = $blast_elem[3];
  my $ident = $blast_elem[4];
  my $e = $blast_elem[5];
  my $bits = $blast_elem[6];
  my $calc_score = $query_length - $hsp_length + $gaps + 1; #Notice that I add 1 to the calc_score since I later need it to evaluate to true, which 0 wouldn't.
  my $q_string = $blast_elem[7];
  my $hit_string = $blast_elem[8];
  my $homo_string = $blast_elem[9];
  my @seq_inds = $blast_elem[10];
  my $hit_strand = $blast_elem[11];
  my $hit_start = $blast_elem[12];
  my $hit_end = $blast_elem[13];
  my $contig_name = $blast_elem[14];
  my $query_strand = $blast_elem[15];
  my $query_start = $blast_elem[16];
  my $hit_length = $blast_elem[17];
  
  #print "$qid, $query_length, $hsp_length, $ident, $e, $bits, $calc_score, $q_string, $hit_string, $homo_string, $hit_strand, $hit_start, $hit_end, $contig_name, $query_strand, $query_start\n"; 
  
  $qid =~ tr/[a-z]/[A-Z]/;
 
  unless (exists $MLST{$qid} and $MLST{$qid} <= $calc_score) {  ##The lines are sorted per query ID so that the HSP with the lowerst e value is first. For each query ID, only the HSP with the lowest calc_score is saved (adk1.., adk1.., adk1.. -> adk1..). 
    $MLST{$qid}  = $calc_score; #%MLST and %PERC_IDENT are later used to find the best query ID at the next level (adk1.., adk2.., adk3.. -> adk2..). 
    $PERC_IDENT{$qid}  = $ident;
    $QUERY_LENGTH {$qid}  = $query_length; #%QUERY_LENGTH, %HSP_LENGTH, and %GAPS are outputted, if the match is not perfect
    $HSP_LENGTH{$qid}  = $hsp_length;
    $GAPS{$qid} = $gaps;
    $Q_STRING{$qid} = $q_string;
    $HIT_STRING{$qid} = $hit_string;
    $HOMO_STRING{$qid} = $homo_string;
    $QUERY_START{$qid} = $query_start;
    $CONTIG_NAME{$qid} = $contig_name ;
    $HIT_STRAND{$qid} = $hit_strand ;
    $HIT_START{$qid} = $hit_start;
    $HIT_END{$qid} = $hit_end;
    $HIT_LENGTH{$qid} = $hit_length;
  }
}

#Declaring variables
my %Profile;
my $Connector;
#Here, the best matching query ID is picked (adk1.., adk2.., adk3.. -> adk2..).
for (sort keys %MLST) {
  /$REG_EXP/;
  $Connector = $2;
  #print "$1, $2, $3\n";
  unless ($Profile{$1}->{score}) { #$1 is for instance adk. Firs time around, we enter this unless..
    $Profile{$1}->{score}  = $MLST{$_};
    $Profile{$1}->{allele} = $1.$Connector.$3;
    #print "Foerste gang: ID: $1 - $_: $MLST{$_}  $PERC_IDENT{$_}\n";
  } if ($MLST{$_} < $Profile{$1}->{score}) { #If there is an allele (another adk) with a better (lower) $calc_score that one is saved
    $Profile{$1}->{score}  = $MLST{$_};
    $Profile{$1}->{allele} = $1.$Connector.$3;
    #print "KUN CALC_SCORE: ID: $1 - $_: $MLST{$_}  $PERC_IDENT{$_}\n";
  } if ($MLST{$_} == $Profile{$1}->{score} and $PERC_IDENT{$_} > $PERC_IDENT{$Profile{$1}->{allele}}) { #If there is an allele (another adk) with an equal $calc_score, but a better (higher) percent ID that one is saved 
    #print "OGSAA PERCENT ID ID: $1 - $_: $MLST{$_}  $PERC_IDENT{$_}\n";
    $Profile{$1}->{score}  = $MLST{$_};
    $Profile{$1}->{allele} = $_;
  }
  
  #The following made the program not function in the right way - should it be implanted in som other way?
  #Because we only want conservative hits, an allele with %ID < 85 and HSP covers < 66% of allele length, is per difinition allele_0 ##New in version 1.3
  #my $length_tester = $HSP_LENGTH{$_} / $QUERY_LENGTH {$_};
  #if (($PERC_IDENT{$_} < 85 || $length_tester < 0.66)){
  #    $Profile{$1}->{allele} = $1.$Connector.0;
  #}
}

# Read profile file
open (TBL, $MLST_DB."/".$Organism.".txt.clean");
my %Table;
my $TableHdr = <TBL>; #The header in the profile table
$TableHdr =~ tr/[a-z]/[A-Z]/;
my @TableHdr = split /\t|\n/, $TableHdr;
my $SequenceType;
while (<TBL>) { #For each of the next lines except the header
  my @Line = split /\t|\n/, $_;
  my @ST = ();
  for (1..$#TableHdr) {
    push @ST, $TableHdr[$_].$Connector.$Line[$_];
  }
  $Table{join("\t", @ST)} = "ST-".$Line[0];
}

# Match the Alleles to an ST and output
my @ST;
for (1..$#TableHdr) {
  #print "\n\nST - $TableHdr[$_]\n\n"; #Debug
  push @ST, $Profile{$TableHdr[$_]}->{allele};
}

#Here, comments are outputtet concerning each of the picked alleles. In addition in this version of the program, the mlst allele sequence and match in the genome is outputtet
foreach my $allele (@ST){ #$allele is for instance adk-6
  #replace '-' and '_' if in the gene name
  # format the gene name by removing the "_" and/or "-" symbols
  #print "allele - $allele<br><br>"; #Debug

  my $gene = $allele;
  if ($gene =~ m/-/)
  {
   my @values = split('-', $gene);
   $gene = $values[0];
  }
  if ($gene =~ m/_/)
  {
   my @values = split('_', $gene);
   $gene = $values[0];
  }
  $GENE_RESULTS_HASH{$gene} = ();
  $GENE_ALIGN_HIT_HASH{$gene} = ();
  $GENE_ALIGN_HOMO_HASH{$gene} = ();
  $GENE_ALIGN_QUERY_HASH{$gene} = ();
  #Declaring variables:
  my $sub_contig;
  my $sub_rev_contig;
  my $final_contig;
  my $spaces_hit;
  my $spaces_match_string;
  my $seq_lower_case;
  my $final_seq_lower_case;
  my $variant = "None";
  my @gaps_in_mlst_allele;
  my @gaps_in_hit;
  my $no_gaps_mlst = 0;
  my $no_gaps_hit = 0;

  my @genArray = split("-", $allele);
  #If two elements were not generated by this (a name and a number), the connector was not "-", and I try "_"
  if (scalar @genArray < 2){
    @genArray = split("_", $allele);
  }
  my $gen = $genArray[0];
  #If there is a perfect match to an MLST allele, it is easy to get the right sequences
  if ($MLST{$allele} == 1 && $PERC_IDENT{$allele} == 100){
     #@#printf "$allele is a perfect match - allele length: $QUERY_LENGTH{$allele}\tHSP length: $HSP_LENGTH{$allele}\tGaps: $GAPS{$allele}\tPercent ID: %6.2f\n", $PERC_IDENT{$allele};
      push(@{$GENE_RESULTS_HASH{$gene}}, "perfect");
      push(@{$GENE_RESULTS_HASH{$gene}}, 100);
      push(@{$GENE_RESULTS_HASH{$gene}}, $QUERY_LENGTH{$allele});
      push(@{$GENE_RESULTS_HASH{$gene}}, $HSP_LENGTH{$allele});
      push(@{$GENE_RESULTS_HASH{$gene}}, $GAPS{$allele});
      push(@{$GENE_RESULTS_HASH{$gene}}, $allele);
     $final_seq_lower_case = $Q_STRING{$allele};
     $final_contig = $HIT_STRING{$allele};
  } 
  #If there is no perfect match to a given pMLST allele, things become a bit more complicated... 
  else {
    #@#printf "*WARNING*: No perfect match for $gen. $allele is the closest match - allele length: $QUERY_LENGTH{$allele}\tHSP length: $HSP_LENGTH{$allele}\tGaps: $GAPS{$allele}\tPercent ID: %6.2f\n" ,$PERC_IDENT{$allele};
    push(@{$GENE_RESULTS_HASH{$gene}}, "warning");
    if ($PERC_IDENT{$allele}!=100){
	  push(@{$GENE_RESULTS_HASH{$gene}},sprintf("%.2f", $PERC_IDENT{$allele}));
      }
    else{
	  push(@{$GENE_RESULTS_HASH{$gene}}, 100);
    }
    push(@{$GENE_RESULTS_HASH{$gene}}, $QUERY_LENGTH{$allele});
    push(@{$GENE_RESULTS_HASH{$gene}}, $HSP_LENGTH{$allele});
    push(@{$GENE_RESULTS_HASH{$gene}}, $GAPS{$allele});
    push(@{$GENE_RESULTS_HASH{$gene}}, $allele);
	
    #Identifying gaps in the pMLST allele string and hit string
    @gaps_in_mlst_allele = Getting_gaps($Q_STRING{$allele});
    @gaps_in_hit = Getting_gaps($HIT_STRING{$allele});
    
    $no_gaps_mlst = scalar @gaps_in_mlst_allele;
    $no_gaps_hit = scalar @gaps_in_hit;


    #Getting the complete mlst allele (even thought it may not all be part of the HSP)   
    my @array_for_getting_mlst_seq = ($allele);
    my $Seqs_ref = grep_ids(-seqs => $Seqs_mlst, -ids => \@array_for_getting_mlst_seq);  
    for (@{ $Seqs_ref }) {
       $seq_lower_case = lc($_->seq);
    }

    #Getting the right contig
    my @array_for_getting_genome_seq = ($CONTIG_NAME{$allele});
    my $Seqs_genome_ref = grep_ids(-seqs => $Seqs_input, -ids => \@array_for_getting_genome_seq);  
    for (@{ $Seqs_genome_ref }) {
       my $contig = lc($_->seq);
       my $length_contig = length($contig);
       
       #Getting the right sub_contig depends on which strand the hit is on
       #If the hit is on the +1 
      if ($HIT_STRAND{$allele} == 1){
         
		if (($QUERY_START{$allele} == 1) && (($QUERY_LENGTH{$allele} + $no_gaps_mlst) == $HSP_LENGTH{$allele})){  
          $variant = 1;
		  $sub_contig = substr($contig, ($HIT_START{$allele} - 1 ),  ($HSP_LENGTH{$allele} - $no_gaps_hit) );
		}
		elsif (($QUERY_START{$allele} == 1) && (($QUERY_LENGTH{$allele} + $no_gaps_mlst) > $HSP_LENGTH{$allele})) { 
		  if (($HIT_START{$allele} + ($QUERY_LENGTH{$allele} + $no_gaps_mlst)) > ($length_contig + $no_gaps_hit)){ 
			$variant = 2;
			#$major_variants_detector = 1;
			$sub_contig = substr($contig, ($HIT_START{$allele} - 1 ),  (($HIT_LENGTH{$allele} + $no_gaps_hit) - $HIT_START{$allele} + 1));
		  }
		  elsif (($HIT_START{$allele} + ($QUERY_LENGTH{$allele} + $no_gaps_mlst)) <= ($length_contig + $no_gaps_hit)) {
			$variant = 3;
			$sub_contig = substr($contig, ($HIT_START{$allele} - 1 ),  ($QUERY_LENGTH{$allele} + $no_gaps_mlst - $no_gaps_hit));
		  }
		}
		elsif ($QUERY_START{$allele} > 1){
		  if (($HIT_START{$allele} - $QUERY_START{$allele}) < 0){
			$variant = 4;
			#$major_variants_detector = 1;
			$sub_contig = substr($contig, 0,  ( $HIT_START{$allele} + (($QUERY_LENGTH{$allele} + $no_gaps_mlst) - $QUERY_START{$allele})));
			#If, as here, the HSP only starts some nucleotides within the mlst allele, a number of spaces must be written before the matching string from the genome. Likewise, the match-string (the "|| ||| ||") should be preceeded by spaces
			$spaces_hit = $QUERY_START{$allele} - $HIT_START{$allele} + 1;
			$spaces_match_string = $QUERY_START{$allele};
		  }
		  else {    
			if ((($HIT_START{$allele} - $QUERY_START{$allele}) + ($QUERY_LENGTH{$allele} + $no_gaps_mlst)) < ($HIT_LENGTH{$allele} + $no_gaps_hit)){
			  $variant = "5a";
			  $sub_contig = substr($contig, ($HIT_START{$allele} - $QUERY_START{$allele}), ($QUERY_LENGTH{$allele} + $no_gaps_mlst));
			  $spaces_match_string = $QUERY_START{$allele};
			}
			else {
			  $variant = "5b";
			  #$major_variants_detector = 1;
			  $sub_contig = substr($contig, ($HIT_START{$allele} - $QUERY_START{$allele}),  ((($HIT_LENGTH{$allele} + $no_gaps_hit) - $HIT_START{$allele}) + $QUERY_START{$allele} -1 ));
			  $spaces_match_string = $QUERY_START{$allele};
			}
		  }
		}
		#else {
		#  print "New option not taken into account!\n";
		#}
	  } 
      #If the hit is on the -1 strand
	  elsif ($HIT_STRAND{$allele} == -1){
		if (($QUERY_START{$allele} == 1) && (($QUERY_LENGTH{$allele} + $no_gaps_mlst) == $HSP_LENGTH{$allele})){
		  $variant = 6;
		  $sub_contig = substr($contig, ($HIT_START{$allele} - 1 ), ($HSP_LENGTH{$allele} - $no_gaps_hit));
        }
		elsif (($QUERY_START{$allele} == 1) && (($QUERY_LENGTH{$allele} + $no_gaps_mlst ) > $HSP_LENGTH{$allele})) {
		  if(($HIT_START{$allele}-(($QUERY_LENGTH{$allele} + $no_gaps_mlst) - ($HSP_LENGTH{$allele} - $no_gaps_hit))) >= 0){
		  $variant = 7;
		  $sub_contig = substr($contig, ($HIT_START{$allele}-(($QUERY_LENGTH{$allele} + $no_gaps_mlst) - ($HSP_LENGTH{$allele} - $no_gaps_hit))) - 1,  (($QUERY_LENGTH{$allele} + $no_gaps_mlst )));
		  }
		  elsif (($HIT_START{$allele}-(($QUERY_LENGTH{$allele} + $no_gaps_mlst ) - ($HSP_LENGTH{$allele} - $no_gaps_hit))) < 0){
			$variant = 8;
			#$major_variants_detector = 1;
			$sub_contig = substr($contig, 0 , ($HIT_START{$allele} + ($HSP_LENGTH{$allele}- $no_gaps_hit) -1));
		  }
		}
		elsif ($QUERY_START{$allele} > 1){
		  if (($HIT_START{$allele} + $QUERY_LENGTH{$allele}) > $length_contig){ 
			#If, as here, the HSP only starts some nucleotides within the mlst allele, a number of spaces must be written before the matching string from the genome
			$variant = 10;
			#$major_variants_detector = 1;
			$spaces_hit = $QUERY_START{$allele} - ($length_contig - $HIT_END{$allele}) ;
			$spaces_match_string = $QUERY_START{$allele};
			$sub_contig = substr($contig, ($HIT_START{$allele} -1 - ($QUERY_LENGTH{$allele} - $HSP_LENGTH{$allele} - $QUERY_START{$allele}+1)), ((($QUERY_LENGTH{$allele} - $QUERY_START{$allele})+($length_contig - $HIT_END{$allele}))+1));
		  }
		  elsif (($HIT_START{$allele} + $QUERY_START{$allele}) <= $length_contig) {
			if (($QUERY_START{$allele} - 1 + $HSP_LENGTH{$allele}) == ($QUERY_LENGTH{$allele} + $no_gaps_mlst)){
			  $variant = "9a";
			  $sub_contig = substr($contig, ($HIT_START{$allele}-1),  ($QUERY_LENGTH{$allele} + $no_gaps_mlst - $no_gaps_hit));
			  $spaces_match_string = $QUERY_START{$allele};
			}
			elsif (($QUERY_START{$allele} + $HSP_LENGTH{$allele}) < $QUERY_LENGTH{$allele}){
			  $spaces_match_string = $QUERY_START{$allele};
			  if ($HIT_START{$allele} - ($QUERY_LENGTH{$allele} - $HSP_LENGTH{$allele} - $QUERY_START{$allele}) < 0) {
				$variant = "9b";
				#$major_variants_detector = 1;
				$sub_contig = substr($contig, 0 , $HSP_LENGTH{$allele} + $QUERY_START{$allele} + $HIT_START{$allele});
			  }
			  elsif ($HIT_START{$allele} - ($QUERY_LENGTH{$allele} - $HSP_LENGTH{$allele} - $QUERY_START{$allele}) >= 0) {
				$variant = "9c";
				$sub_contig = substr($contig, ($HIT_START{$allele} -2 - ($QUERY_LENGTH{$allele} - $HSP_LENGTH{$allele} - $QUERY_START{$allele} )),  $QUERY_LENGTH{$allele});
			  }
			}
		  }
		  #$final_contig = $sub_contig;  
		}
		else {    
		  print "New option not taken into account!\n";       
        }
		$sub_rev_contig = reverse($sub_contig); #NOTE I SHOULD THINK ABOUT WHEN TO ADD THE "-" IN THE CREATION OF THE final_contig. BEFORE OR AFTER sub_contig is REVERSED?
		$sub_rev_contig =~ tr/acgt/tgca/;
		$sub_contig = $sub_rev_contig;
	  }
    }
    #Adding gaps to the sub_contig (if there are any) leading to the creation of final_contig
    if ($no_gaps_hit > 0){
      my $hsp_length1 = (length $sub_contig) + $no_gaps_hit;
      my $sign1;
      my $flag1 = 0;
      my $start1 = 0;
      for (my $i = 0 ; $i < $hsp_length1 ; ++$i){
		$flag1 = 0;
        foreach my $pos (@gaps_in_hit){
		  if ($variant != 4){
			if ($i == ($pos + $QUERY_START{$allele} -1)){
              $sign1 = "-";
              $flag1 = 1;
              $start1 = $start1 - 1;
            }
		  }
		  else {
			if ($i == ($pos + $HIT_START{$allele} -1)){
			  $sign1 = "-";
              $flag1 = 1;
              $start1 = $start1 - 1;
            }
		  }
        }
        unless ($flag1 == 1) {
          $sign1 = substr($sub_contig, $start1 , 1);
        }
        $final_contig .= $sign1;
        ++$start1;
      }
    }
    else {
      $final_contig = $sub_contig;
    }
    #Adding gaps to the mlst allele sequence, $seq_lower_case (if there are any) leading to the creation of final_seq_lower_case
    if ($no_gaps_mlst > 0){
      my $hsp_length2 = (length $seq_lower_case) + $no_gaps_mlst;
      my $sign2;
      my $flag2 = 0;
      my $start2 = 0;
      for (my $i = 0 ; $i < $hsp_length2 ; ++$i){
        $flag2 = 0;
        foreach my $pos (@gaps_in_mlst_allele){
          if ($i == ($pos + $QUERY_START{$allele} -1)){
            $sign2 = "-";
            $flag2 = 1;
            $start2 = $start2 - 1;
          }
        }
        unless ($flag2 == 1) {
          $sign2 = substr($seq_lower_case, $start2 , 1);
        }
        $final_seq_lower_case .= $sign2;
        ++$start2;
      }
    }   
    else {
      $final_seq_lower_case = $seq_lower_case;
    }
  } 
  #Nicely printing the sequences
  #print "(Variant: $variant)\n";
  for (my $j = 0 ; $j < $QUERY_LENGTH{$allele} ; $j+=60){
    #Printing the pMLST allele
	#print "pMLST allele seq: ";
	my $mlst_substr = substr($final_seq_lower_case, $j , 60);
    #print $mlst_substr . "\n";
	push(@{$GENE_ALIGN_QUERY_HASH{$gene}}, $mlst_substr);
	 
	#Printing spaces before the match string with the "|| |||| || ||" and the string itself
	#print "                 ";
	my $string_spaces_match_string = "";  #For saving the spaces as a string of spaces instead of a number of spaces
	for (my $i = 1 ; $i < $spaces_match_string ; ++$i){
	  $string_spaces_match_string .= " ";  
	}
	my $homo_string =  $string_spaces_match_string . $HOMO_STRING{$allele};
    my $homo_string_substr = substr($homo_string, $j , 60);
    #print $homo_string_substr . "\n";
	push(@{$GENE_ALIGN_HOMO_HASH{$gene}}, $homo_string_substr);
	 
	#Printing the match in the genome
	my $string_spaces_hit = "";  #For saving the spaces as a string of spaces instead of a number of spaces
	for (my $i = 1 ; $i < $spaces_hit ; ++$i){
	  $string_spaces_hit .= " ";  
	}
	#print "Hit in genome:   ";
	my $string_final_contig = $string_spaces_hit . $final_contig;
    my $string_final_contig_substr = substr($string_final_contig, $j , 60);
    #print $string_final_contig_substr . "\n\n";
	push(@{$GENE_ALIGN_HIT_HASH{$gene}}, $string_final_contig_substr);
	#print "GAHA - $gene  #  $string_final_contig_substr<br><br>";
  }
}


#If there are no known ST, this is outputtet
my $SeqType;
#print join("\t", @ST); #Debug
if ($Organism eq "IncF" or $Organism eq "incf") {
  $SeqType = join("\t", @ST);
}
else {
  $SeqType = $Table{join("\t", @ST)};
}
#print "Seqtype: $SeqType\n"; #Debug
unless ($SeqType =~m/^ST-.+/){
    print "";
  #print "Unknown ST";
}

#print "ST OFormat ::\t".$OFormat."\n";
#print $Table{join("\t", @ST)}, "\t", join("\t", @ST)."\n";

#Here, it is examined if the ST is associated with a clonal complex. If so, the clonal complex is outputtet
if ($SeqType =~m/^ST-[\d]+/){
  my @split_SeqType = split('-',$SeqType);
  my $OnlyNo = $split_SeqType[1];
  open (FH, $MLST_DB.'/'.$Organism.'.txt.clpx');
  while (defined ( my $line = <FH>)){
    my @split_line = split (' ' , $line);
    if ($split_line[0] == $OnlyNo){
      print "Clonal complex: " . $split_line[1] . "\n";
    }
  }
}

#print Dumper(%GENE_RESULTS_HASH);

#let's check the hash content
if ($SeqType eq ""){
    push(@RESULTS_AND_SETTINGS_ARRAY, "Unknown ST");
}
else{
    push(@RESULTS_AND_SETTINGS_ARRAY, $SeqType);
}
push(@RESULTS_AND_SETTINGS_ARRAY, $Organism);
push(@RESULTS_AND_SETTINGS_ARRAY, $mlstProfiles{$Organism});
push(@RESULTS_AND_SETTINGS_ARRAY, $InFile);


#############################################################


#### TXT PRINTING ####
my $seqtype = print_txt_results(\@RESULTS_AND_SETTINGS_ARRAY, \%GENE_RESULTS_HASH, \%GENE_ALIGN_QUERY_HASH, \%GENE_ALIGN_HOMO_HASH, \%GENE_ALIGN_HIT_HASH);

exit;


# --------------------------------------------------------------------
# %% Land of the Subroutines %%
#

###################################
sub commandline_parsing {
    while (scalar @ARGV) {
        if ($ARGV[0] =~ m/^-d$/) {
            $MLST_DB = $ARGV[1];
            shift @ARGV;
            shift @ARGV;
        }
        elsif ($ARGV[0] =~ m/^-b$/) {
            $BLAST = $ARGV[1];
            $BLASTALL = "$BLAST/bin/blastall";
            $FORMATDB = "$BLAST/bin/formatdb";
            shift @ARGV;
            shift @ARGV;
        }
        elsif ($ARGV[0] =~ m/^-s$/) {
            $Organism = $ARGV[1];
            shift @ARGV;
            shift @ARGV;
        }
        elsif ($ARGV[0] =~ m/^-i$/) {
            $InFile = $ARGV[1];
            shift @ARGV;
            shift @ARGV;
        }
        elsif ($ARGV[0] =~ m/^-o$/) {
            $dir = $ARGV[1];
            mkdir $dir;
            shift @ARGV;
            shift @ARGV;
        }
        elsif ($ARGV[0] =~ m/^-h$/) {
            $Help = 1;
            shift @ARGV;
        }
        else {
         &print_help();
         exit
        }
    }
}

# Run blast and parse output
# Arguments should be a hash with arguments to blast in option => value format
# Returns text lines of blast output

sub get_blast_run {
   my ($tmp_dir, $org, %args) = @_;
   #my $fh = $tmp_dir;
   my $file = "blast_$org.fsa";
   #my ($fh, $file) = tempfile( DIR => '/tmp', UNLINK => 1);
   output_sequence(-file => ">$tmp_dir/$file", seqs => delete $args{-d}, -format => 'fasta');
   die "Error! Could not build blast database" if (system("$FORMATDB -p F -i $tmp_dir/$file"));
   system("rm -r formatdb.log");
   system("rm -r  $tmp_dir/$file.n*");
   my $query_file = "$file.blastpipe";

   #open QUERY, ">> $query_file" || die("Error! Could not perform blast run");
   output_sequence(-file => ">$tmp_dir/$query_file", seqs => $args{-i}, -format => 'fasta');
   #close QUERY;
   
   delete $args{-i};

   my $cmd = join(" ", %args);
   my $file2 = "$tmp_dir/$file.blast_output";
   system("$BLASTALL -d $tmp_dir/$file -i $tmp_dir/$query_file -o $file2 $cmd");

   my $report = new Bio::SearchIO( -file   => $file2,
                                   -format => "blast"
                                 );
   # Go through BLAST reports one by one
   my @blast;
   while(my $result = $report->next_result) {
      # Go through each matching sequence
      while(my $hit = $result->next_hit)    {
         # Go through each each HSP for this sequence
         while (my$hsp = $hit->next_hsp)  {
            push(@blast, $result->query_accession ."\t".
                        $result->query_length ."\t".
                        $hsp->hsp_length ."\t".
                        $hsp->gaps ."\t".
                        $hsp->percent_identity ."\t".
                        $hsp->evalue ."\t".
                        $hsp->bits ."\t".
                        $hsp->query_string ."\t".
                        $hsp->hit_string ."\t".
                        $hsp->homology_string ."\t".
                        $hsp->seq_inds ."\t".
                        $hsp->strand('hit') ."\t".
                        $hsp->start('hit') ."\t".
                        $hsp->end('hit') ."\t".
                        $hit->name ."\t".
                        $hsp->strand('query') ."\t".
                        $hsp->start('query') ."\t".
                        $hit->length ."\n");
         }
      }
   }
	system("rm -r $tmp_dir/$file");
   system("rm -r $tmp_dir/$query_file");
   return @blast;
}

###################################
# Finds sequences with specific ids in an array of Bio::Seq objects
# Args:
#   -seqs => A reference to an array with Bio::Seq objects
#   -ids  => A reference to an array of ids
#   -v    => Like for grep, if defined, return the ids which didn't match
# Returns
#   An array (or array reference) with the Bio::Seq objects having the requested ids

sub grep_ids {
  my %argv = @_;
  my %ids;
  for my $id (@{ $argv{-ids} }) {
    $ids{$id} = 1;
  }
  my @out;
  for my $seq (@{ $argv{-seqs} }) {
    if (exists $ids{$seq->id()}) {
      push @out, $seq unless (defined $argv{-v});
    } elsif (defined $argv{-v}) {
      push @out, $seq;
    }
  }
  return wantarray ? @out : \@out;
}


###################################
# Reads one sequence file in a format supported by BioPerl
# Arguments should be an array:
#   Filenames to be loaded, ideally the @ARGV array
#   -fh         => The filehandle to read from, defaults to ARGV
#   -format     => The file format, defaults to "fasta"
#   <...>       => Additional options to Bio::SeqIO
# Returns:
#   A reference to an array of Bio::Seq objects

sub read_seqs {
  my %args      = @_;
  $args{-fh}    = \*ARGV unless (exists $args{-fh} or exists $args{-file});
  my (@seqs, %ids);
  $args{-format} = "fasta" unless (defined $args{-format});
  my $seq_in = Bio::SeqIO->new(%args);
  while (my $seq = $seq_in->next_seq) {
    push @seqs, $seq;
  }
  return (wantarray ? @seqs : \@seqs);
}



#####################################
#This sub takes a nucleotide string as input and identifies gaps ("-").
#An array with positions of the gaps is returned
sub Getting_gaps {
  my $input_string = $_[0];
  my @split_input_string = split('',$input_string);
  my @gap_positions;
  for (my $i = 0 ; $i < (scalar @split_input_string) ; ++$i){
    if ($split_input_string[$i] eq "-"){
      push(@gap_positions,$i);
    }
  }
  
  return (@gap_positions);
}



###################################
# Output in sequence formats supported by bioperl
# Arguments should be a hash:
#   seqs    => A reference to an array of sequences. Values are ids or Bio::Seq objects
#   tempdir => If set, writes to a temp file in specified directory instead. Uses File::Temp
#   Any additional arguments will be forwarded to the Bio::SeqIO->new() call.
#      If these don't include either -fh or -file, STDOUT will be used (but see tempdir)
# Returns:
#   The filehandle and filename

sub output_sequence {
   my %args = @_;
   my $seqs_ref = delete $args{seqs};
   my $i = 1;
   #$args{-fh} = \*STDOUT unless (exists $args{-fh} or exists $args{-file});
   #if (exists $args{tempdir}) {
   #   my $tempdir = delete $args{tempdir};
   #   ($args{-fh}, $args{-file}) = tempfile(DIR => $tempdir, SUFFIX => ".".$args{-format})
   #}
   my $file = delete $args{-file} if (exists $args{-fh} && exists $args{-file}); # Stupid BioPerl cannot handle that both might be set...
   print %args;
   my $seq_out = Bio::SeqIO->new(%args);
   $args{-file} = $file if (defined $file);
   for my $seq (@{ $seqs_ref }) {
      $seq_out->write_seq($seq);
   }
   return ($args{-fh}, $args{-file});
}



# --------------------------------------------------------------------
# %% Help Page/Documentation %%
#

sub print_help {
  my $ProgName     = PROGRAM_NAME;
  my $ProgNameLong = PROGRAM_NAME_LONG;
  my $Version      = VERSION;
  my $CMD = join(" ", %ARGV);
  print <<EOH;

NAME
	$ProgName - $ProgNameLong

SYNOPSIS
	$ProgName [Options]
	
DESCRIPTION
	Calculates the pMLST profile based on a BLAST alignment of the input
	sequence file and the specified allele set. If possible the ST will be
	given, or if unknown, that field will be left empty

        Notice that although the options mimic that the input sequences are
	aligned against the alleles, it is in fact the other way around. First,
	the input is converted to a blast database, against which is aligned the
	alleles from the species specified with '-s'.

	Notice also that the default options for BLAST are changed to suit the
	pMLST alignment. 

OPTIONS

	-h HELP
                    Prints a message with options and information to the screen
    -d DATABASE
                    The path to where you have located the database folder
    -b BLAST
                    The path to the location of blast-2.2.26 if it is not added
                    to the user's path (see the install guide in 'README.md')
    -i INFILE
                    Your input file which needs to be preassembled partial
                    or complete genomes in fasta format
    -o OUTFOLDER
                    The folder you want to have your output files places.
                    If not specified the program will create a folder named
                    'Output' in which the result files will be stored
    -s SCHEME
                    The pMLST scheme you want to use. Details about the
                    available schemes can be found in the config file

Example of use with the 'database' folder located in the current directory and Blast added to the user's path
    
    perl pMLST-1.4.pl -i INFILE.fasta -o OUTFOLDER -s incf 

Example of use with the 'database' and 'blast-2.2.26' folders loacted in other directories

    perl pMLST-1.4.pl -d path/to/database -b path/to/blast-2.2.26 -i INFILE.fasta -o OUTFOLDER -s incf 
    
VERSION
    Current: $Version

AUTHORS
    Carsten Friis, carsten\@cbs.dtu.dk, Mette Voldby Larsen, metteb\@cbs.dtu.dk

EOH
}

#ADD SPACES BEFORE AND AFTER GIVEN STRING TO FIT MINIMUM AMOUNT OF SPACES
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}
sub AlignLeft { #string, minimumStrLenght 
	my($str, $num) = @_;
	my $returnString = $str;
	# Adding spaces after string
	my $spacesneeded = $num-length($str);
	while ($spacesneeded > 0) {
		$returnString .= " ";
		$spacesneeded--;
	}
	return $returnString;
}
sub AlignCenter { #string, minimumStrLenght 
	my($str, $num) = @_;
	my $returnString = "";
	my $spacesneeded = &roundup(($num-length($str))/2);
	# Adding spaces before string
	while ($spacesneeded > 0) {
		$returnString .= " ";
		$spacesneeded--;
	}
	# Adding spaces after string
	$returnString .= $str;
	$spacesneeded = int(($num-length($str))/2);
	while ($spacesneeded > 0) {
		$returnString .= " ";
		$spacesneeded--;
	}
	return $returnString;
}
sub AlignRight { #string, minimumStrLenght 
	my($str, $num) = @_;
	my $returnString = "";
	my $spacesneeded = $num-length($str);
	# Adding spaces before string
	while ($spacesneeded > 0) {
		$returnString .= " ";
		$spacesneeded--;
	}
	$returnString .= $str;
	return $returnString;
}
sub formatName {
	my($nam) = @_;
	if(index($nam,'-') > -1){
      my @split = split('-',lc($nam));
		if(exists $mlstProfiles{$split[0]}){
		  return join('_',$mlstProfiles{$split[0]},uc($split[1]));
		}else{
		  return lc($nam);
		}
	}elsif(index($nam,'_') > -1){
      my @split = split('_',lc($nam));
		if(exists $mlstProfiles{$split[0]}){
		  return join('_',$mlstProfiles{$split[0]},uc($split[1]));
		}else{
		  return lc($nam);
		}
	}else{
		return exists $mlstProfiles{lc($nam)} ? $mlstProfiles{lc($nam)} : lc($nam);
	}
}


# --------------------------------------------------------------------
# print_txt_results creates 3 files: Hit_in_genome_seq.fsa, pMLST_allele_seq.fsa and results.txt
#
sub print_txt_results{
# %% results.txt is a text result table and list of alleles  %%
# Generates a tab separated text table containing the scripts results and alignment
# eg.
# pMLST Results                                                                   
#                                                                                 
# Sequence Type: Unknown ST                                                       
# pMLST Profile: incf                                                             
#                                                                                 
# --------------------------------------------------------------------------------
#    GENE       % IDENTITY    HSP Length    Allele Length   GAPS     BEST MATCH   
# --------------------------------------------------------------------------------
#    FIA           95.24          21             329          0           FIA_6   
#    FIB           95.45          22             373          0           FIB_1   
#    FIC           94.74          19             200          0           FIC_1   
#    FII           91.30          23             156          0          FII_25   
#    FIIK         100.00          13              82          0          FIIK_4   
#    FIIS          94.12          17             212          0          FIIS_4   
#    FIIY          95.00          20             164          0          FIIY_1   
#    HCM1178AC     95.00          20             164          0     HCM1178AC_1   
# --------------------------------------------------------------------------------
#                                                                                 
#                                                                                 
# Please note that one or more loci do not match perfectly to any previously      
# registered pMLST allele. We recommend verifying the results by traditional      
# methods for pMLST!                                                              
#                                                                                 
# fia: WARNING, ID: 95.24%, HSP/Length: 21/329, Gaps: 0, Best match: fia_6        
#                                                                                 
# pMLST allele seq:   ttgccctcagaagaagctgttttttcacaaagttatccctgcttattgactcttttttat
#                                    ||||||||||||||| ||||!                        
# Hit in genome:      ggctcaaaagtgcggtctgttttttcacaaatttatcaagtactttatcaaaatgcccat

	my ($resultsAndSettingsArray, $geneResultsHash, $geneAlignQueryHash, $geneAlignHomoHash, $geneAlignHitHash) = @_;
	#print "print_txt_results::START\n"; #DEBUG
	my $txtresults = "";
   my $tabr = "";
	my $allelealign = "";
	my $hits = "";
   my $scheme = @{$resultsAndSettingsArray}[2];
	
	# INITIALIZING WARNING VARIABLE
	my $txtwarning = 0;
	my $stwarning = 0;
	
	#if one or more alleles do not match perfectly, a -like* will be adde to the sequence results
	foreach my $key (sort { if(lc($a) eq 'fii'){return -1;}else{return lc($a) cmp lc($b);}} (keys %{$geneResultsHash})) {
		my $array = ${$geneResultsHash}{$key};
		if ($scheme ne "IncF" and $scheme ne "incf") {
		  if( @$array[0] ne "perfect" ){ $stwarning = 1; }
		}
	}
	
	#Finding the right ST for IncF
	my @FAB = ('F','A','B');
	my @IAB = ('-','-','-');
	my $st;
	
	if ($scheme eq "IncF" or $scheme eq "incf") {
	  my @ST = split("\t", @{$resultsAndSettingsArray}[0]);
	  my $i = 0;
	  
	  foreach my $key (sort { if(lc($a) eq 'fii'){return -1;}else{return lc($a) cmp lc($b);}} (keys %{$geneResultsHash})) {
		my $array = ${$geneResultsHash}{$key};
		my $identity = @$array[1];
		my $hspLen = @$array[3];
		my $allLen = @$array[2];
		
		if ($ST[$i] =~ m/FII\_([A-Za-z]*)(\d+)/) {
		  if ($identity == 100 and $allLen > 0 and $hspLen / $allLen == 1) {
			if ($1 ne '') {
               $FAB[0] = $1;
            }
			$IAB[0] = $2;
		  }
		  elsif ($identity > 85 and $allLen > 0 and $hspLen / $allLen > 0.66) {
			if ($1 ne '') {
               $FAB[0] = $1;
            }
            $IAB[0] = "$2*";
		  }
		}
		if ($ST[$i] =~ m/FIA\_([A-Za-z]*)(\d+)/) {
		  if ($identity == 100 and $allLen > 0 and $hspLen / $allLen == 1) {
			if ($1 ne '') {
               $FAB[1] = $1;
            }
			$IAB[1] = $2;
		  }
		  elsif ($identity > 85 and $allLen > 0 and $hspLen / $allLen > 0.66) {
			if ($1 ne '') {
               $FAB[1] = $1;
            }
            $IAB[1] = "$2*";
		  }
		}
		if ($ST[$i] =~ m/FIB\_([A-Za-z]*)(\d+)/) {
		  if ($identity == 100 and $allLen > 0 and $hspLen / $allLen == 1) {
			if ($1 ne '') {
               $FAB[2] = $1;
            }
			$IAB[2] = $2;
		  }
		  elsif ($identity > 85 and $allLen > 0 and $hspLen / $allLen > 0.66) {
			if ($1 ne '') {
               $FAB[2] = $1;
            }
            $IAB[2] = "$2*";
		  }
		}
		
		$i++;
	  }
	  
	  #Writing the ST in the right way
	  $st ='[';
	  $i = 0;
	  while ($i <= 2) {
		$st .= join('', $FAB[$i],$IAB[$i]);
		if ($i!= 2){
		  $st .= ":";
		}	
		$i++;
	  }
	  $st .= "]";
	  
	  if ($st =~ m/^\[F-:A-:B-]$/) {
		$st = "Unknown ST";
	  }
	  
	}
	else {
	  $st =  @{$resultsAndSettingsArray}[0] ;
	}
	
	## PRINTING HEADER / SETTINGS
	$txtresults .= "pMLST Results\n\n";
   $txtresults .= "pMLST Profile: ".@{$resultsAndSettingsArray}[2]."\n";
	
	if ($stwarning == 1 and $st ne "Unknown ST") {
	   $txtresults .= "Sequence Type: Unknown ST\n";
	   $txtresults .= "Closest match: $st\n";
		
		$tabr .= "Sequence Type: Unknown ST\n";
	   $tabr .= "Closest match: $st\n";
	}
	else {
		$txtresults .= "Sequence Type: ".$st."\n";
		$tabr .= "Sequence Type: ".$st."\n";

	} 
	#$txtresults .= "pMLST Profile: ".@{$resultsAndSettingsArray}[2]."\n\n";

	# PRINTING RESULT TABLE

	#$txtresults .= "--------------------------------------------------------------------------------\n";
   $tabr .= "Gene\t% Identity\tHSP Length\tAllele Length\tGaps\tBest match\n";
	$txtresults .= "********************************************************************************\n";
	$txtresults .= "   GENE       % IDENTITY    HSP Length    Allele Length   GAPS     BEST MATCH   \n";
	$txtresults .= "********************************************************************************\n";
	foreach my $key (sort { if(lc($a) eq 'fii'){return -1;}else{return lc($a) cmp lc($b);}} (keys %{$geneResultsHash})) {
	  my $array = ${$geneResultsHash}{$key};
	  my $locus = &formatName($key);
	  my $identity =@$array[1];
	  my $hspLen =@$array[3];
	  my $allLen =@$array[2];
	  my $gaps =@$array[4];
	  my $matchAll;
	  if (@$array[5] eq "") {
	    $matchAll = "$key". "_0";  
	  } else {
	    $matchAll = &formatName(@$array[5]);
	  }		

	  # ADDING SPACES TO LOCUS' WHICH ARE SHORTER THAN 6 CHARS
	  my $spacesneeded = 6-length($key);
	  while ($spacesneeded > 0) {
	    $locus .= " ";
	    $spacesneeded--;
	  }
	  
	  # Removing the result if the proportion of allele length that is
	  # found in hit is less than 66% or identity is less than 85%
	  #if ($Conservative_match && (($identity <= 85) || ($hspLen / $allLen <= 0.66)) ){
	  if (($identity <= 85) || ($hspLen / $allLen <= 0.66) ){
	    ($identity, $hspLen, $allLen, $gaps) = ('-  ', '-', '-', '-');  # $matchAll and the last '-'was removed, as we have defined matchAll earlier
	  }else{
		$identity = sprintf("%.2f", $identity);
	  }
		
     $tabr .= "$locus\t$identity\t$hspLen\t$allLen\t$gaps\t$matchAll\n";
	  $txtresults .=  "   ".&AlignLeft($locus,12)." ".&AlignRight($identity,6)."       ".&AlignCenter($hspLen,7)."         ".&AlignCenter($allLen,7)."       ".&AlignCenter($gaps,3)." ".&AlignRight($matchAll,14)."\n";
		
	  if( @$array[0] ne "perfect" ){ $txtwarning = 1; } # WARNING IS ADDED IF NOT ALL MATCHES ARE PERFECT
	}#end foreach
	
	$txtresults .=  "================================================================================\n\n";
	
	# PRINTING WARNING (if any!)
	if( $txtwarning == 1 ){
	   $txtresults .= "* Please note that one or more loci do not match perfectly to any previously\n".
		               "registered pMLST allele. We recommend verifying the results by traditional\n".
						  	"methods for pMLST!\n\n\nExtended Output:\n".
							"--------------------------------------------------------------------------------\n\n";
	}
	
	# PRINTING THE EXTENDED OUTPUT
	foreach my $key (sort { if(lc($a) eq 'fii'){return -1;}else{return lc($a) cmp lc($b);}} (keys %{$geneResultsHash})) {

    my $array = ${$geneResultsHash}{$key};
    my $identity =@$array[1];
    my $hspLen =@$array[3];
    my $allLen =@$array[2];
    my $gaps =@$array[4];
    my $matchAll = &formatName(@$array[5]);
	 my $qStart = $QUERY_START{uc($matchAll)};
    my $qend = $qStart + $hspLen - 1;

    # print header line for each allele
    my $outStr = &formatName($key);
    if (@$array[0] eq "perfect" ){
      $outStr .= ": PERFECT MATCH, ";
    }elsif (@$array[0] eq "warning" ){
      $outStr .= ": WARNING, ";
    }

    #if (! $Conservative_match || (($identity > 85) && ($hspLen / $allLen > 0.66)) ){
	 if (($identity > 85) && ($allLen > 0 and $hspLen / $allLen > 0.66) ){
      $txtresults .= $outStr."ID: ".$identity."%, HSP/Length: ".$hspLen."/".$allLen.", Gaps: ".$gaps.", Best match: ".$matchAll."\n\n";
      $allelealign .= ">".$matchAll."\n";
      $hits .= ">".$outStr."ID: ".$identity."%, HSP/Length: ".$hspLen."/".$allLen.", Pos: ".$qStart."..".$qend.", Gaps: ".$gaps.", Best match: ".$matchAll."\n";
      #now print the alleles
      my $queryArray = ${$geneAlignQueryHash}{$key};
      my $homoArray = ${$geneAlignHomoHash}{$key};
      my $hitArray = ${$geneAlignHitHash}{$key};

      for (my $i=0; $i < scalar(@$hitArray); $i++){
        my $tmpQuerySingleLine = @$queryArray[$i];
        my $tmpHomoSingleLine = @$homoArray[$i];
        my $tmpHitSingleLine = @$hitArray[$i];

        $txtresults .= "pMLST allele seq:   ".$tmpQuerySingleLine."\n";
        $txtresults .= "                    ".$tmpHomoSingleLine."\n";
        $txtresults .= "Hit in genome:      ".$tmpHitSingleLine."\n\n";
        $allelealign .= $tmpQuerySingleLine."\n";
        $hits .= $tmpHitSingleLine."\n";
		}#end for
		$txtresults .= "\n--------------------------------------------------------------------------------\n\n";
    }
  }#end foreach
	$tabr .= "pMLST Profile: ".@{$resultsAndSettingsArray}[2]."\n";

	#WRITING results.txt
	open (TXTRESULTS, '>'."$dir/results.txt") || die("Error! Could not write to results.txt");
	print TXTRESULTS $txtresults;
	close (TXTRESULTS);

	#WRITING Hit_in_genome_seq.fsa
	open (HIT, '>'."$dir/Hit_in_genome_seq.fsa") || die("Error! Could not write to Hit_in_genome_seq.fsa");
	print HIT $hits;
	close (HIT);

	#WRITING pMLST_allele_seq.fsa
	open (ALLELE, '>'."$dir/pMLST_allele_seq.fsa") || die("Error! Could not write to pMLST_allele_seq.fsa");
	print ALLELE $allelealign;
	close (ALLELE);
   
   #WRITING standard_output.tab
   open (TABRESULTS, '>'."$dir/".'results_tab.txt') || die("Error! Could not write to results_tab.txt");
   print TABRESULTS $tabr;
   close (TABRESULTS);
	
	return $st;
}#end sub(print_txt_results)