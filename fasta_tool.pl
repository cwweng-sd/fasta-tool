# This script used to parse the FASTA file for various respects of processing.
#
# Date:2016/11/7
# Upgrade:2017/7/12
# Author:C.W.Weng (Jeff)
 
#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

my $opts=parse_params();
do_pileup_to_fasta($opts);

exit;

#----------------

sub error{
    my(@msg)=@_;
    if(scalar@msg){croak(@msg);}
    die
        "Usage: fasta_tool.pl [OPTIONS] < in.fa\n",
        "\n",
	"Options:\n",
        "   -h,  -?, --help                               This help message.\n",
        "   -p,  --prerequisite-process                   If the input of FASTA file doesn't be represented the format of each sequence on a single line,\n",
        "                                                 the flag must be used to convert the regular FASTA file.\n",
        "                                                 Each sequence will be re-labelled as the following format:seq1, seq2, seq3 and so on.\n",
        "\n",
        "   -c,  --count-bases                            Calculate the ratio of bases.\n",
        "\n",
        "   -f,  --fetch-seq                              Fetch the applicable sequences that match the forward/reverse primers in the 5'/3' end,\n",
        "                                                 and then trim primers.\n",
        "\n",
        "   -g,  --get-seq                                Get target sequences from the label list.\n",
        "                                                 The label of target sequences need be listed per line in a file,\n",
        "                                                 furthermore, the file must be named seq_label.txt.\n",
        "\n",
        "   -l,  --cal-len                                Print the length of sequence.\n",
        "\n",
        "   -m,  --format-output-fasta                    Convert the input of FASTA file to standard FASTA format,\n",
        "                                                 that is represented as a series of lines, each of which be 70 nucleotides.\n",
        "\n",
        "   -s,  --fetch-subseq                           Fetch specified region of sequence from the input FASTA file.\n",
        "\n",
        "   -sn, --search-nu-seq-occu [sequence pattern]  Retrieve the position of last occurrence of a match of nucleotide sequence.\n",
        "\n",
        "   -sa, --sort-asc-by-length                     Sort the length of sequences by ascending order.\n",
	"\n",
	"   -sd, --sort-desc-by-length                    Sort the length of sequences by descending order.\n",
	"\n",
	"If the flag of -f is used, please provide the primer file including forward and reserve sequences.\n",
        "Please refer the example file, i.e. primer.txt.\n",
        "\n",
        "If the flag of -s is used, please provide the specified region of sequences.\n",
	"Please refer the example file, i.e. chr_position.txt.\n",
        "\n",
	"If the flag of -g is used, please provide the label of target sequences.\n",
        "Please refer the example file, i.e. seq_label.txt.\n",
	"\n";
}

sub parse_params{
    my %opts=();

    $opts{fh_in}=*STDIN;

    while(my $arg=shift(@ARGV)){
        if($arg eq '-p'  || $arg eq '--prerequisite-process') { $opts{convert_fasta}=1; next; }
        if($arg eq '-f'  || $arg eq '--fetch-seq') { $opts{fetch_seq}=1; next; }
        if($arg eq '-g'  || $arg eq '--get-seq') { $opts{get_seq}=1; next; }
        if($arg eq '-l'  || $arg eq '--cal-len') { $opts{cal_len}=1; next; }
        if($arg eq '-m'  || $arg eq '--format-output-fasta') { $opts{format_fasta}=1; next; }
        if($arg eq '-c'  || $arg eq '--count-bases') { $opts{count_base}=1; next; }
        if($arg eq '-s'  || $arg eq '--fetch-subseq') { $opts{fetch_subseq}=1; next; }
        if($arg eq '-sn' || $arg eq '--search-nu-seq-occu') { $opts{search_nu_seq_occu}=1; next; }
	if($arg eq '-sa' || $arg eq '--sort-asc-by-length') { $opts{sort_asc_by_length}=1; next; }
        if($arg eq '-sd' || $arg eq '--sort-desc-by-length') { $opts{sort_desc_by_length}=1; next; }
	if($arg eq '-?'  || $arg eq '-h' || $arg eq '--help') { error(); }

        if($arg=~m/\D+/){
            if(check_fit_iupac_nu($arg)){ $opts{sequence_pattern}=$arg; next; }
            else{ print "Please input the correct IUPAC format of nucletoide sequence.\n"; next; }
        }

        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }

    return \%opts;
}

sub do_pileup_to_fasta{
    my ($opts)=@_;

    my $convert_fasta=$$opts{convert_fasta} ? 1 : 0;
    my $fetch_seq=$$opts{fetch_seq} ? 1 : 0;
    my $get_seq=$$opts{get_seq} ? 1 : 0;
    my $cal_len=$$opts{cal_len} ? 1 : 0; 
    my $format_fasta=$$opts{format_fasta} ? 1 : 0;
    my $count_base=$$opts{count_base} ? 1 : 0;
    my $fetch_subseq=$$opts{fetch_subseq} ? 1 : 0;
    my $search_nu_seq_occu=$$opts{search_nu_seq_occu} ? 1 : 0;
    my $sort_asc_by_length=$$opts{sort_asc_by_length} ? 1 : 0;
    my $sort_desc_by_length=$$opts{sort_desc_by_length} ? 1 : 0;

    if($convert_fasta){ convert_fasta_format(); }
    if($fetch_seq){ fetch_seq_pool(); }
    if($get_seq){ get_target_seq(); }
    if($cal_len){ cal_seq_len(); }
    if($format_fasta){ format_output_fasta(); }
    if($count_base){ cal_base_ratio(); }
    if($fetch_subseq){ fetch_sub_seq(); }
    if($search_nu_seq_occu){ search_nu_seq_occu(); }
    if($sort_asc_by_length){ sort_asc_by_length(); }
    if($sort_desc_by_length){ sort_desc_by_length(); }
}

sub convert_fasta_format{
    open(OUTPUT,">converted_fasta.fa")||die "Unable to output the converted FASTA file!\n";
    my $fh_in=$$opts{fh_in};
    my $seq="";
    my $seq_label;
    my $seq_count=0;

    while(my $line=<$fh_in>){
        chomp($line);
        my $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq ">"){
            if($seq ne ""){
                print OUTPUT ">$seq_label\n$seq\n";
            }
            $seq_count+=1;
            $seq_label="seq".$seq_count;
            $seq="";
        }else{
            $seq=$seq.$line;
        }
        if(eof()){
            print OUTPUT ">$seq_label\n$seq\n";
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub fetch_seq_pool{
    open(OUTPUT1,">extracted_seq.fa")||die "Unable to output the applicable sequences!\n";
    open(OUTPUT2,">extracted_trim_seq.fa")||die "Unable to output the applicable trimming sequences!\n";

    # process the primer file
    my @primers=parse_primers();
    my $forward_primer=$primers[0];
    my $revcomp_forward_primer=$primers[1];
    my $reverse_primer=$primers[2];
    my $revcomp_reverse_primer=$primers[3];
    my @forward_primer_ary=();
    my @revcomp_forward_primer_ary=();
    my @reverse_primer_ary=();
    my @revcomp_reverse_primer_ary=();
    if(length($forward_primer)==length($reverse_primer)){
        @forward_primer_ary=trim_primer($forward_primer);
        @revcomp_forward_primer_ary=trim_primer_revcomp($revcomp_forward_primer);
        @reverse_primer_ary=trim_primer($reverse_primer);
        @revcomp_reverse_primer_ary=trim_primer_revcomp($revcomp_reverse_primer);
        if($#forward_primer_ary==$#revcomp_forward_primer_ary && $#forward_primer_ary==$#reverse_primer_ary &&
           $#forward_primer_ary==$#revcomp_reverse_primer_ary){
            for(my $i=0;$i<=$#forward_primer_ary;$i++){
                print length($forward_primer_ary[$i])." mer\n";
                print "Forward primer: $forward_primer_ary[$i]\n";
                print "Revcomp of forward primer: $revcomp_forward_primer_ary[$i]\n";
                print "Reverse primer: $reverse_primer_ary[$i]\n";
                print "Revcomp of reverse primer: $revcomp_reverse_primer_ary[$i]\n";
            }
        }
    }
    print "\n";
    print "-----------------------------------------";
    print "\n";

    # fetch the applicable sequences
    my $fh_in=$$opts{fh_in};
    my $fasta_pattern;
    my $seq_title;
    my $seq;
    my $seq_len;
    my $trim_seq;
    my $trim_seq_len;
    my $trim_pair_seq;
    my $prime5seq_front_match_fprimer;
    my $prime5seq_end_match_fprimer;
    my $prime3seq_front_match_rprimer_revcomp;
    my $prime3seq_end_match_rprimer_revcomp;
    my $prime5seq_front_match_rprimer;
    my $prime5seq_end_match_rprimer;
    my $prime3seq_front_match_fprimer_revcomp;
    my $prime3seq_end_match_fprimer_revcomp;
    while(my $line=<$fh_in>){
        chomp($line);
        $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq '>'){
            $seq_title=substr $line,1;
        }else{
            $seq=$line;
            $seq_len=length($seq);
            if($#forward_primer_ary==$#revcomp_forward_primer_ary && $#forward_primer_ary==$#reverse_primer_ary &&
               $#forward_primer_ary==$#revcomp_reverse_primer_ary){
                my @prime_5_seq_front=();
                my @prime_5_seq_end=();
                my @prime_3_seq_front=();
                my @prime_3_seq_end=();
                for(my $i=0;$i<=$#forward_primer_ary;$i++){
                    # for @prime_5_seq_front, @prime_3_seq_end
                    my $sub_seq_len_posi_front_start=$i;
                    my $sub_seq_len_posi_front_end=length($forward_primer_ary[$i]);
                    my $sub_seq_len_nega_end_start=-length($forward_primer_ary[0]);
                    my $sub_seq_len_nega_end_end=length($forward_primer_ary[$i]);
                    my $prime_5_seq_front=substr $seq,$sub_seq_len_posi_front_start,$sub_seq_len_posi_front_end;
                    my $prime_3_seq_end=substr $seq,$sub_seq_len_nega_end_start,$sub_seq_len_nega_end_end;
                    push @prime_5_seq_front,$prime_5_seq_front;
                    push @prime_3_seq_end,$prime_3_seq_end;
                    # for @prime_5_seq_end, @prime_3_seq_front
                    my $sub_seq_len_posi_end_start=0;
                    my $sub_seq_len_posi_end_end=length($forward_primer_ary[0])-$i;
                    my $sub_seq_len_nega_front_start=-length($forward_primer_ary[$i]);
                    my $prime_5_seq_end=substr $seq,$sub_seq_len_posi_end_start,$sub_seq_len_posi_end_end;
                    my $prime_3_seq_front=substr $seq,$sub_seq_len_nega_front_start;
                    push @prime_5_seq_end,$prime_5_seq_end;
                    push @prime_3_seq_front,$prime_3_seq_front;
                }
                print "$seq_title\n";
                if($#forward_primer_ary==$#prime_5_seq_front){
                    $prime5seq_front_match_fprimer="false";
                    print "prime5seq-front to forward: ";
                    for(my $i=0;$i<=$#forward_primer_ary;$i++){
                        if($prime_5_seq_front[$i] eq $forward_primer_ary[$i]){
                            $prime5seq_front_match_fprimer="true";
                            print length($prime_5_seq_front[$i])."\t";
                            if(length($seq)==$seq_len){
                               $trim_seq=substr $seq,length($prime_5_seq_front[$i]);
                               $trim_seq_len=length($trim_seq);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#forward_primer_ary==$#prime_5_seq_end){
                    $prime5seq_end_match_fprimer="false";
                    print "prime5seq-end to forward: ";
                    for(my $i=0;$i<=$#forward_primer_ary;$i++){
                        if($prime_5_seq_end[$i] eq $forward_primer_ary[$i]){
                            $prime5seq_end_match_fprimer="true";
                            print length($prime_5_seq_end[$i])."\t";
                            if(length($seq)==$seq_len){
                               $trim_seq=substr $seq,length($prime_5_seq_end[$i]);
                               $trim_seq_len=length($trim_seq);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#reverse_primer_ary==$#prime_5_seq_front){
                    $prime5seq_front_match_rprimer="false";
                    print "prime5seq-front to reverse: ";
                    for(my $i=0;$i<=$#reverse_primer_ary;$i++){
                        if($prime_5_seq_front[$i] eq $reverse_primer_ary[$i]){
                            $prime5seq_front_match_rprimer="true";
                            print length($prime_5_seq_front[$i])."\t";
                            if(length($seq)==$seq_len){
                               $trim_seq=substr $seq,length($prime_5_seq_front[$i]);
                               $trim_seq_len=length($trim_seq);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#reverse_primer_ary==$#prime_5_seq_end){
                    $prime5seq_end_match_rprimer="false";
                    print "prime5seq-end to reverse: ";
                    for(my $i=0;$i<=$#reverse_primer_ary;$i++){
                        if($prime_5_seq_end[$i] eq $reverse_primer_ary[$i]){
                            $prime5seq_end_match_rprimer="true";
                            print length($prime_5_seq_end[$i])."\t";
                            if(length($seq)==$seq_len){
                               $trim_seq=substr $seq,length($prime_5_seq_end[$i]);
                               $trim_seq_len=length($trim_seq);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#revcomp_forward_primer_ary==$#prime_3_seq_front){
                    $prime3seq_front_match_fprimer_revcomp="false";
                    print "prime3seq-front to revcomp forward: ";
                    for(my $i=0;$i<=$#revcomp_forward_primer_ary;$i++){
                        if($prime_3_seq_front[$i] eq $revcomp_forward_primer_ary[$i]){
                            $prime3seq_front_match_fprimer_revcomp="true";
                            print length($prime_3_seq_front[$i])."\t";
                            if(length($trim_seq)==$trim_seq_len){
                               $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_front[$i]);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#revcomp_forward_primer_ary==$#prime_3_seq_end){
                    $prime3seq_end_match_fprimer_revcomp="false";
                    print "prime3seq-end to revcomp forward: ";
                    for(my $i=0;$i<=$#revcomp_forward_primer_ary;$i++){
                        if($prime_3_seq_end[$i] eq $revcomp_forward_primer_ary[$i]){
                            $prime3seq_end_match_fprimer_revcomp="true";
                            print length($prime_3_seq_end[$i])."\t";
                            if(length($trim_seq)==$trim_seq_len){
                               $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_end[$i]);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#revcomp_reverse_primer_ary==$#prime_3_seq_front){
                    $prime3seq_front_match_rprimer_revcomp="false";
                    print "prime3seq-front to revcomp reverse: ";
                    for(my $i=0;$i<=$#revcomp_reverse_primer_ary;$i++){
                        if($prime_3_seq_front[$i] eq $revcomp_reverse_primer_ary[$i]){
                            $prime3seq_front_match_rprimer_revcomp="true";
                            print length($prime_3_seq_front[$i])."\t";
                            if(length($trim_seq)==$trim_seq_len){
                               $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_front[$i]);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if($#reverse_primer_ary==$#prime_3_seq_end){
                    $prime3seq_end_match_rprimer_revcomp="false";
                    print "prime3seq-end to revcomp reverse: ";
                    for(my $i=0;$i<=$#revcomp_reverse_primer_ary;$i++){
                        if($prime_3_seq_end[$i] eq $revcomp_reverse_primer_ary[$i]){
                            $prime3seq_end_match_rprimer_revcomp="true";
                            print length($prime_3_seq_end[$i])."\t";
                            if(length($trim_seq)==$trim_seq_len){
                               $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_end[$i]);
                            }
                            last;
                        }
                    }
                    print "\n";
                }
                if(($prime5seq_front_match_fprimer eq "true" && $prime3seq_end_match_rprimer_revcomp eq "true") ||
                   ($prime5seq_front_match_rprimer eq "true" && $prime3seq_end_match_fprimer_revcomp eq "true") ||
                   ($prime5seq_front_match_fprimer eq "true" && $prime3seq_front_match_rprimer_revcomp eq "true") ||
                   ($prime5seq_front_match_rprimer eq "true" && $prime3seq_front_match_fprimer_revcomp eq "true") ||
                   ($prime5seq_end_match_fprimer eq "true" && $prime3seq_front_match_rprimer_revcomp eq "true") ||
                   ($prime5seq_end_match_rprimer eq "true" && $prime3seq_front_match_fprimer_revcomp eq "true") ||
                   ($prime5seq_end_match_fprimer eq "true" && $prime3seq_end_match_rprimer_revcomp eq "true") ||
                   ($prime5seq_end_match_rprimer eq "true" && $prime3seq_end_match_fprimer_revcomp eq "true")){
                    print OUTPUT1 ">$seq_title\n$seq\n";
                    print OUTPUT2 ">$seq_title\n$trim_pair_seq\n";
                    print "output\n";
                }
                print "\n";
            }
        }
    }

    close $fh_in;
    close OUTPUT1;
    close OUTPUT2;
}

sub cal_seq_len{
    open(OUTPUT,">fasta_seq.len")||die "Unable to output the length of sequence!\n";
    my $fh_in=$$opts{fh_in};
    my $fasta_pattern;
    my $seq_title;
    my $seq;
    my $seq_len;

    while(my $line=<$fh_in>){
        chomp($line);
        $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq '>'){
            $seq_title=substr $line,1;
        }else{
            $seq=$line;
            $seq_len=length($seq);
            print OUTPUT "$seq_title\t$seq_len\n";
        }
    }
    
    close $fh_in;
    close OUTPUT;
}

sub get_target_seq{
    open(TARGET_SEQ,"seq_label.txt")||die "Please provide the list of sequence labels in the current directory.\n";
    my @target_seq=<TARGET_SEQ>;
    close TARGET_SEQ;

    open(OUTPUT,">target_seq.fa")||die "Unable to output the target sequences!\n";
    my $fh_in=$$opts{fh_in};
    my $fasta_pattern;
    my $seq_title;
    my $seq;
    
    while(my $line=<$fh_in>){
        chomp($line);
        $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq '>'){
            $seq_title=substr $line,1;
        }else{
            $seq=$line;
            foreach my $target(@target_seq){
                chomp($target);
                if($seq_title eq $target){
                    print OUTPUT ">$seq_title\n$seq\n";
                }
            }
        }
    }
    
    close $fh_in;
    close OUTPUT;
}

sub cal_base_ratio{
    open(OUTPUT,">base_ratio.txt")||die "Unable to output the ratio of bases!\n";
    my $fh_in=$$opts{fh_in};
    my $fasta_pattern;
    my $seq;
    my $seq_title;
    my $count_A;
    my $count_T;
    my $count_C;
    my $count_G;
    my $count_N;

    # calculate the ratio of bases
    while(my $line=<$fh_in>){
        chomp($line);
        $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq '>'){
            $seq_title=substr $line,1;
        }else{
            $seq=$line;
            $count_A=0;
            $count_T=0;
            $count_C=0;
            $count_G=0;
            $count_N=0;
            for(my $i=0;$i<length($seq);$i++){
                my $base=substr $seq,$i,1;
		if($base eq "A" || $base eq "a"){ $count_A=$count_A+1; }
                if($base eq "T" || $base eq "t"){ $count_T=$count_T+1; }
                if($base eq "C" || $base eq "c"){ $count_C=$count_C+1; }
                if($base eq "G" || $base eq "g"){ $count_G=$count_G+1; }
                if($base eq "N" || $base eq "n"){ $count_N=$count_N+1; }
            }
            # output the ratio of bases
            print OUTPUT "$seq_title\tA:$count_A\tT:$count_T\tC:$count_C\tG:$count_G\tN:$count_N\n";
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub fetch_sub_seq{
    open(SUBSEQ_REGION,"chr_position.txt")||die "Please provide the specified region of sequences in the current directory.\n";
    my @subseq_region=<SUBSEQ_REGION>;
    close SUBSEQ_REGION;

    open(OUTPUT,">subseq.fa")||die "Unable to output the specified region of sequences!\n";
    my $fh_in=$$opts{fh_in};
    my $fasta_pattern;
    my $seq_title;
    my $seq;

    while(my $line=<$fh_in>){
        chomp($line);
        $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq '>'){
            $seq_title=substr $line,1;
        }else{
            $seq=$line;
            foreach my $subseq_reg(@subseq_region){
                chomp($subseq_reg);
                my @region_info=split '\t',$subseq_reg;
                if($region_info[0] eq "sequence"){
                    next;
                }
                my $subseq_len=$region_info[2]-$region_info[1];
                my $subseq_start=$region_info[1]-1;
                if($subseq_len>0){
                    if($seq_title eq $region_info[0]){
                        my $subseq=substr $seq,$subseq_start,$subseq_len;
                        print OUTPUT ">$region_info[0]:$region_info[1]-$region_info[2]\n$subseq\n";
                    }
                }elsif($subseq_len<=0){
                    print "Current parsed sequence is $seq_title. ";
                    print "Focus on $region_info[0]:$region_info[1]-$region_info[2], ";
                    print "the start position of specified region larger than the end position! ";
                    print "Please modify the 'chr_position.txt' file.\n";
                    last;
                }
            }
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub format_output_fasta{
    open(OUTPUT,">formatted_fasta.fa")||die "Unable to output the formatted FASTA file!\n";
    my $fh_in=$$opts{fh_in};
    my $fasta_pattern;
    my $seq;

    while(my $line=<$fh_in>){
        chomp($line);
    	$fasta_pattern=substr $line,0,1;
   	if($fasta_pattern eq '>'){
    	    print OUTPUT "$line\n";
        }else{
            $seq=$line;
            my $output_len=70;
            while(my $chunk=substr($seq,0,$output_len,"")){
                print OUTPUT "$chunk\n";
            }
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub search_nu_seq_occu{
    open(OUTPUT,">seq_occu.txt")||die "Unable to output the position of last occurrence of a match of nucleotide sequence!\n";
    my $fh_in=$$opts{fh_in};
    my $nu_seq_pattern=$$opts{sequence_pattern};
    my $fasta_pattern;
    my $seq_title;
    my $seq;

    print OUTPUT "\t$nu_seq_pattern\n";
    while(my $line=<$fh_in>){
        $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq '>'){
            chomp($line);
            $seq_title=substr $line,1;
            print OUTPUT "$seq_title\t";
        }else{
            $seq=$line;
            chomp($seq);
            if($seq=~m{($nu_seq_pattern+)}g){
                my $pos=pos($seq)-length $1;
                print OUTPUT "$pos\n";
            }else{
                print OUTPUT "na\n";
            }
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub parse_primers{
    my @total_primers=();
    my $forward_primer;
    my $revcomp_forward_primer;
    my $reverse_primer;
    my $revcomp_reverse_primer;
    open(PRIMER_FILE,"primer.txt")||die "Please provide the primer file in the current directory.\nThe example primer file can be a reference.\n";
    my @primers=<PRIMER_FILE>;
    close PRIMER_FILE;
    foreach my $primer(@primers){
        chomp($primer);
        my $primer_direction=substr $primer,0,1;
        if($primer_direction eq 'F'){
            $primer=~s/F://g;
            $forward_primer=$primer;
            $revcomp_forward_primer=revcomp($forward_primer);
            next;
        }
        if($primer_direction eq 'R'){
            $primer=~s/R://g;
            $reverse_primer=$primer;
            $revcomp_reverse_primer=revcomp($reverse_primer);
            next;
        }
    }
    push @total_primers,$forward_primer;
    push @total_primers,$revcomp_forward_primer;
    push @total_primers,$reverse_primer;
    push @total_primers,$revcomp_reverse_primer;
    # retern the primers: forward primer, revcomp of forward primer, reverse primer, revcomp of reverse primer
    return @total_primers;
}

sub check_fit_iupac_nu{
    my $input_seq=$_[0];
    my %base_type;
    my %iupac_nu=iupac_nu_info();
    my $count_nonfit=0;
    my $check_status="true";
    for(my $i=0;$i<length($input_seq);$i++){
        my $base=substr $input_seq,$i,1;
        $base_type{$base}=$base;
    }
    for my $base(keys %base_type){
        if(defined($iupac_nu{$base})){ next;
        }else{ $count_nonfit+=1; }
    }
    if($count_nonfit>0){ $check_status="false"; }
    return $check_status;
}

sub trim_primer{
    my $primer=$_[0];
    my $mer_len=abs(3-length($primer));
    my @primer_ary=();
    for(my $i=0;$i<=$mer_len;$i++){
        my $primer_seq=substr $primer,$i;
        push @primer_ary,$primer_seq;
    }
    return @primer_ary;
}

sub trim_primer_revcomp{
    my $primer=$_[0];
    my $mer_len=abs(3-length($primer));
    my @revcomp_primer_ary=();
    for(my $i=0;$i<=$mer_len;$i++){
        my $revcomp_primer_seq=substr $primer,0,length($primer)-$i;
        push @revcomp_primer_ary,$revcomp_primer_seq;
    }
    return @revcomp_primer_ary;
}

sub revcomp{
   my $primer=$_[0];
   $primer=~tr/ACGTacgt/TGCAtgca/;
   $primer=reverse($primer);
   return $primer; 
}

sub sort_asc_by_length{
    open(OUTPUT,">sort_asc_fasta.fa")||die "Unable to output the sorted FASTA file by ascending order!\n";
    my $fh_in=$$opts{fh_in};
    my $seq;
    my $seq_label;
    my $seq_len;
    my %seqlabel2seq;
    my %seqlabel2len;

    while(my $line=<$fh_in>){
        chomp($line);
        my $fasta_pattern=substr $line,0,1;
        if($fasta_pattern eq ">"){
            $seq_label=substr $line,1;
	}else{
            $seq=$line;
	    $seq_len=length($seq);
            $seqlabel2seq{$seq_label}=$seq;
	    $seqlabel2len{$seq_label}=$seq_len;
	}
    }
    foreach my $seqlabel(sort {$seqlabel2len{$b} <=> $seqlabel2len{$a}} keys %seqlabel2len){
        print OUTPUT ">$seqlabel\n";
	print OUTPUT "$seqlabel2seq{$seqlabel}\n";
    }

    close $fh_in;
    close OUTPUT;
}

sub sort_desc_by_length{
    open(OUTPUT,">sort_desc_fasta.fa")||die "Unable to output the sorted FASTA file by descending order!\n";
    my $fh_in=$$opts{fh_in};
    my $seq;
    my $seq_label;
    my $seq_len;
    my %seqlabel2seq;
    my %seqlabel2len;

    while(my $line=<$fh_in>){
        chomp($line);
	my $fasta_pattern=substr $line,0,1;
	if($fasta_pattern eq ">"){
	    $seq_label=substr $line,1;
	}else{
	    $seq=$line;
	    $seq_len=length($seq);
	    $seqlabel2seq{$seq_label}=$seq;
	    $seqlabel2len{$seq_label}=$seq_len;
        }
    }
    foreach my $seqlabel(sort {$seqlabel2len{$a} <=> $seqlabel2len{$b}} keys %seqlabel2len){
        print OUTPUT ">$seqlabel\n";
	print OUTPUT "$seqlabel2seq{$seqlabel}\n";
    }

    close $fh_in;
    close OUTPUT;
}

sub iupac_nu_info{
    my %iupac_nu=(
        "A" => "Adenine",     "C" => "Cytosine",    "G" => "Guanine",     "T" => "Thymine",
        "U" => "Uracil",      "R" => "A or G",      "Y" => "C or T",      "S" => "G or C",
        "W" => "A or T",      "K" => "G or T",      "M" => "A or C",      "B" => "C or G or T",
        "D" => "A or G or T", "H" => "A or C or T", "V" => "A or C or G", "N" => "any base",
        "." => "gap",         "-" => "gap");
    return %iupac_nu;
}

sub iupac_amino_acid_info{
    my %iupac_aa_code=(
        "A" => "Alanine",       "C" => "Cysteine",  "D" => "Aspartic Acid", "E" => "Glutamic Acid",
        "F" => "Phenylalanine", "G" => "Glycine",   "H" => "Histidine",     "I" => "Isoleucine",
        "K" => "Lysine",        "L" => "Leucine",   "M" => "Methionine",    "N" => "Asparagine",
        "P" => "Proline",       "Q" => "Glutamine", "R" => "Arginine",      "S" => "Serine",
        "T" => "Threonine",     "V" => "Valine",    "W" => "Tryptophan",    "Y" => "Tyrosine");
    return %iupac_aa_code; 
}

sub iupac_amino_acid_code_info{
    my %iupac_aa_three_letter_code=(
        "A" => "Ala", "C" => "Cys", "D" => "Asp", "E" => "Glu",
        "F" => "Phe", "G" => "Gly", "H" => "His", "I" => "Ile",
        "K" => "Lys", "L" => "Leu", "M" => "Met", "N" => "Asn",
        "P" => "Pro", "Q" => "Gln", "R" => "Arg", "S" => "Ser",
        "T" => "Thr", "V" => "Val", "W" => "Trp", "Y" => "Tyr");
    return %iupac_aa_three_letter_code;
}
