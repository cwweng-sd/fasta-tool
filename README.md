# fasta-tool
<br>
This perl script can be used to parse the FASTA file for various respects of processing.<br>
<br>
Date: 2016/11/7<br>
Upgrade: 2017/7/12<br>
Author: Jeff<br>
<br>
Usage: fasta_tool.pl [OPTIONS] < in.fa <br>
<br>
Options:<br>
-h, -?, --help<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This help message.<br>
-p, --prerequisite-process<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If the input of FASTA file doesn't be represented the format of each sequence on a single line, the flag must be used to convert the regular FASTA file. Each sequence will be re-labelled as the following format:seq1, seq2, seq3 and so on.<br>
-c, --count-bases<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Calculate the ratio of bases.<br>
-f, --fetch-seq<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fetch the applicable sequences that match the forward/reverse primers in the 5'/3' end, and then trim primers.<br>
-g, --get-seq<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Get target sequences from the label list. The label of target sequences need be listed per line in a file, furthermore, the file must be named seq_label.txt.<br>
-l, --cal-len<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Print the length of sequence.<br>
-m,  --format-output-fasta<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Convert the input of FASTA file to standard FASTA format, that is represented as a series of lines, each of which be 70 nucleotides.<br>
-s, --fetch-subseq<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fetch specified region of sequence from the input FASTA file.<br>
-sn, --search-nu-seq-occu [sequence pattern]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Retrieve the position of last occurrence of a match of nucleotide sequence.<br>
-sa, --sort-asc-by-length<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sort the length of sequences by ascending order.<br>
-sd, --sort-desc-by-length<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sort the length of sequences by descending order.<br>
<br>
If the flag of -f is used, please provide the primer file including forward and reserve sequences.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please refer the example file, i.e. primer.txt.<br>
<br>
If the flag of -s is used, please provide the specified region of sequences.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please refer the example file, i.e. chr_position.txt.<br>
<br>
If the flag of -g is used, please provide the label of target sequences.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please refer the example file, i.e. seq_label.txt.<br>
