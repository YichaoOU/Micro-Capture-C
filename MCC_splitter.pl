#!/usr/bin/perl -w
use strict;
#use Data::Dumper;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use IO::Handle;
use Data::Dumper;

# Specifications
my $f_counter = 0;
my %hash;
my @line_labels = qw(name seq spare qscore); # FASTQ line labels
my %psldata;
my %counters;
my %gene_hash;
my $min_read_length = 25;
my $help=0;
my $man=0;
my $reference = "_MNAse_multi_";
my $parameter_reference = ""; #Generally keep the reference blank unless you are going to run the same data twice with different parameters.
my $current_directory = cwd(); 
my $repcapture = 1;
my $dump=0; #flag to output unaligned reads to a file
my $use_limit=0; #whether to limit the script to analysing the first n lines
my $end_flag=0;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); #defines the start time of the script (printed out into the report file later on)


=head1 NAME

MCC_splitter.pl

=head1 SYNOPSIS

 This script uses a psl and FASTQ files as input
 It reads the gene name out of the psl file (which is determined by the original fasta name in the reference 'genome' of oligo sequences. The 'genome' fa file name needs to be in this format >chr11:32196518-32196638_Hba-2C
 The script splits the fastq file into multiple fastq files for each gene depending the number of input files (either Read1 or Read2 or a single merged file can be used)
 
=head1 EXAMPLE

 # perl /t1-data/user/hugheslab/jdavies/00Scripts/MCC_splitter.pl -f /t1-data/user/jdavies/a10_MNAse/MNAse_extended.fastq -p /t1-data/user/jdavies/a10_MNAse/MNAse_extended_300.psl -r MNAse_ext_300
 MNAse_psl_multi.pl -f input_R1.fq -f input_R2.fq -p psl_file_R1 -p psl_file_R2 -r reference
	perl /research/rgs01/home/clusterHome/yli11/Programs/Micro-Capture-C/MCC_splitter.pl -f test_ext.fastq -p test_ext.psl -r test -limit 0 -all

=head1 OPTIONS

 -f		Input FASTQ filenames
 -p		Input PSL files (BLAT output)
 -dir 	Output directory
 -r		Reference - a name for the experiment
 -all	Flag that tells the script to output both the capture and non capture reads.

=head1 AUTHOR

 Written by James Davies 2019

=cut



# The GetOptions from the command line
&GetOptions
(
	"f=s"=>\ my @full_fastq_files,	# -f1		Input filename 
	"p=s"=>\ my @psl_files,			# -p1 psl filename 1	
	"dir=s"=>\ $current_directory, 	# overwrites output directory over the current directory
	"r=s"=>\ $reference, 			# reference
	"pr=s"=>\ $parameter_reference, # parameter reference allows you to rerun the code on files with same names - used when opening the FASTQ files
	"all"=>\ $repcapture,			# -nosplit - report the whole read - do not remove the target sequence from the read
	"limit=i"=>\ $use_limit,			# -limit		Limit the analysis to the first n reads of the file
	"dump"=>\ $dump,
);

pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;
pod2usage(2) unless ($full_fastq_files[0]); 


# Creates a folder for the output files to go into - this will be a subdirectory of the file that the script is in

my $output_path= "$current_directory/$reference";
if (-d $output_path){}
else {mkdir $output_path};

my @FHIN; #Stores input file handles
my @FHDUMP; #Stores dump file handles

#This part of the script derives the names for the output fq files from the input file names and stores them in the array @fastq_files
#Splits out the filename from the path
#print Dumper (\@fastq_files);

my @fastq_files;

for (my$i=0; $i<=$#full_fastq_files; $i++)
{
	my $input_filename=$full_fastq_files[$i];
	my $input_path = "";
	if ($input_filename =~ /(.*)\/(\V++)/) {$input_path = $1; $input_filename = $2};
	print "Input fastq files include: $input_filename\n";
	unless ($input_filename =~ /(.*).fastq/) {die"filename does not match .fastq"};
	my $file_name=$1;
	push @fastq_files, $file_name;
}


my $report_file_name = "$output_path/$fastq_files[0]"."\_report_file.txt";
print "Report file: $report_file_name\n";
#Opens a report file based on the name of the first input fastq file
open(FHOUT, ">", $report_file_name) or die "Cannot open report file";
#open(TESTFH, ">$output_path/$fastq_files[0]\_debug_file.txt") or die "Cannot open report file \n"; #For debugging the script
print "Script $0 run at: $mday/$mon/2$year $hour:$min:$sec\n";
print FHOUT "Script $0 run at: $mday/$mon/2$year $hour:$min:$sec\n#### If you don't see a time stamp for the script finishing it has crashed somewhere! ####\n\n";
#################################################################################################################################################################################################################################
# Opens the PSL files outputted by blat for reads 1 and 2 and stores the data in a hash
#
#print Dumper (\@psl_files);

if (scalar(@psl_files)==0){die "no psl files specified"} #checks that there is a psl file specified
if (scalar(@psl_files) != scalar(@full_fastq_files)){die "The number of fastq and psl files doesn't match"} #checks that the number of fastq and psl files match

my @PSLIN;
my $current_read_name;
my $last_read_name = "first";
print FHOUT "Opening PSL files....\n";
#open(PSLTMP, ">$output_path/psl.tmp") or die "Cannot open psl file psl.tmp$!\n";
# Opens the FASTQ input files from the command line
for (my$i=0; $i<=$#psl_files; $i++)
{
	open $PSLIN[$i], "$psl_files[$i]" or die "Cannot open input fastq file $psl_files[$i]\n";
	print FHOUT "Opened $psl_files[$i]\n";
} #opens the FASTQ input files

FHOUT->flush(); #Flushes the file handle so that you know what is going on.
# Reads in the first n lines of the psl file - this needs to be enough that it will include a blat to all of the targets because the output fastq files are opened using the entries in the psldata hash. 
psl_read(20000);


#################################################################################################################################################################################################################################
# Opens a FASTQ files to output the data into for.
# This will open the same number of FASTQ files as are used for the input for each target in the psl file (i.e. 2 if R1 and R2 files are used or only one if combined reads are used)
# The file handles are stored in the hash 
print FHOUT "Opening FASTQ files\n";
my @keys = sort keys(%gene_hash);
#print Dumper (\@keys);
for (my$i=0; $i<=$#keys; $i++)
{
	for (my$j=0; $j<=$#fastq_files; $j++)
	{
		open $gene_hash{$keys[$i]}{"FH"}{$fastq_files[$j]}, ">$output_path/$keys[$i]\_$parameter_reference$fastq_files[$j].fastq" or die "can't open output fastq file $output_path/".$keys[$i]."_$fastq_files[$j]" ; #opens R1 output files
	}
}

#################################################################################################################################################################################################################################

# Opens the FASTQ input files from the command line
for (my$i=0; $i<=$#full_fastq_files; $i++)
{
	print FHOUT "Opening $full_fastq_files[$i]....\n";
	open $FHIN[$i], "$full_fastq_files[$i]" or die "Cannot open input fastq file $full_fastq_files[$i]\n";
	#Opens a dump file based on the name of the first input fastq file
	if ($dump==1){open($FHDUMP[$i], ">$output_path/$fastq_files[$i]\_dump.fastq") or die "Cannot open dump file for $fastq_files[$i]\n";}
} #opens the FASTQ input files

FHOUT->flush(); #Flushes the file handle so that you know what is going on.

# Beginning of loop through the fastq files
while ($hash{$fastq_files[0]}{$line_labels[$f_counter]} = readline $FHIN[0])  #assigns each fq line to the hash in batches of 4
{

for (my$i=1; $i<=$#fastq_files; $i++){$hash{$fastq_files[$i]}{$line_labels[$f_counter]}= readline $FHIN[$i]} # puts the data from the other fastq files into the hash 

$f_counter++; $counters{"01 FQ Lines read:"}++;

if ($f_counter==4) # when the four lines of the fastq file for each read have been read in
{
    $counters{"01 FQ reads read:"}++;
	for (my$i=0; $i<=$#fastq_files; $i++){name_trim($hash{$fastq_files[$i]}{"name"}, $fastq_files[$i], \%hash)} # trims off the last part of the read names so they should all be the same
    for (my$i=1; $i<=$#fastq_files; $i++){unless ($hash{$fastq_files[0]}{"name"} eq $hash{$fastq_files[$i]}{"name"}){die "fastq_names don't match"}} #checks that the read names are the same
    my $read_name = $hash{$fastq_files[0]}{"name"}; # specifies the read name from the first fastq file
	
	if (($counters{"01 PSL reads in memory"}<10000) and ($end_flag==0)){$end_flag=psl_read(10000)} #Loads 10000 lines from the psl file if the buffer of lines in the memory falls below 10000. The flag ==1 when the end of the psl file is reached to stop it from trying to load more lines when this happens.
	if (exists $psldata{$read_name}) #Checks whether there is a match for the read in psl file
	{
		# print $read_name;
		# print "###############";
		# print ;
		# print Dumper($psldata{$read_name});
		# print Dumper(\%psldata);
		# print "----------------";
		if($psldata{$read_name}{"best_match"}>15) #This asks if the best match is > than the target.  The best match is independent of R1 or 2
		{
		$counters{"01 FQ reads with a Blat match:"}++;
		
		my $gene = $psldata{$read_name}{"gene"};
		$counters{"04a $gene FQ reads with a PSL match:"}++;
		
			if (exists($gene_hash{$gene}{"FH"})) #Occasionally someone puts extra oligo targets into the BLAT_fasta file which appear in the PSL file at low frequency and are not detected in the first 20k lines - this picks up this possibility
			{
				for (my$i=0; $i<=$#fastq_files; $i++)  #Loops through the fastq files and outputs split reads from the file with the matches in the PSL file 
				{
					my $new_name = "_Gene>".$psldata{$read_name}{"gene"}.$psldata{$read_name}{"strand"}; #defines the new readname including the UMI and the index sequences
					
					#print FHOUT "$read_name\t$gene\n"; #unhash for debugging - can be useful for finding the line that is causing the script to crash
					
					if (exists $psldata{$read_name}{"$fastq_files[$i]"})  #checks whether the fq file is the match in the psldata hash 
				    {
						$counters{"01 Split read output $fastq_files[$i]:"}++;
						$counters{"04b $gene FQ reads with a capture:"}++;
						split_read_output_all_out2($gene_hash{$gene}{"FH"}{$fastq_files[$i]},$read_name ,$fastq_files[$i], \%hash, \%psldata, $new_name, \%counters, $gene) #default output is to output the whole read
					}
					else
					{
						$psldata{$read_name}{$fastq_files[$i]}{"basestring"}="Y," x length($hash{$fastq_files[$i]}{"seq"});  #if the read doesn't have a match it loads Y's into the psl hash for outputting - this is important when using R1 and R2 because one of the reads would ideally not have a blat match if there is a ligation junction in the middle unsequenced region.
						split_read_output_all_out2($gene_hash{$gene}{"FH"}{$fastq_files[$i]},$read_name ,$fastq_files[$i], \%hash, \%psldata, $new_name, \%counters, $gene);
						$counters{"01 Reads with a match but no output $fastq_files[$i]:"}++;
						$counters{"04b $gene FQ reads with a match but no capture:"}++;
					}
				}
			}
			else{$counters{"!! Error $gene not initially identified in first 20000 lines!!"}++}
		}
		else {$counters{"01 Reads with PSL match <15:"}++;}
		delete $psldata{$read_name};
		$counters{"01 PSL reads in memory"}--;
	}
	else
	{
		$counters{"01 Non matching FQ reads:"}++;
		# print $read_name."\t";my @keys = keys (%psldata);print $keys[0]."\n"; #to test that the psl and fastq reads names are being handled properly
		if ($dump==1)
		{
			for (my$i=0; $i<=$#fastq_files; $i++){read_output_dump($FHDUMP[$i], $read_name ,$fastq_files[$i], \%hash, \%counters)}
		}
	} #stores a record of the non matching read indices
    $f_counter=0
}

if (($use_limit !=0) and ($counters{"01 FQ Lines read:"} >($use_limit*4))){print "ending script early\n"; last}
}



if (exists($counters{"01 FQ reads with a Blat match:"}) and exists($counters{"01 FQ reads read:"})){$counters{"01 Capture efficiency"}=100*$counters{"01 FQ reads with a Blat match:"}/$counters{"01 FQ reads read:"};}
else{$counters{"01 Capture efficiency"}=0; print FHOUT "!!ERROR..... no PSL or FQ reads read!!\n\n"}
output_hash (\%counters, \*FHOUT);

# Closes the file handles 
for (my$i=0; $i<=$#keys; $i++)
{
	for (my$j=0; $j<=$#fastq_files; $j++)
	{
		close $gene_hash{$keys[$i]}{"FH"}{$fastq_files[$j]} or warn "can't close output fastq file $output_path/".$keys[$i]."_$fastq_files[$j]" ; #closes the file handles
	}
}

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
print FHOUT "Script finished at: $mday/$mon/2$year $hour:$min:$sec";
close FHOUT;
exit;


################################################################################################################################################################################################################################
# Subroutines
#
################################################################################################################################################################################################################################
# This subroutine splits the reads based on the array containing Ys and digits for the number of the blat match from the string generated by the match string subroutine

sub split_read_output_all_out2 
{
  my ($fh_ref, $read_name ,$fq_file, $hash_ref, $psldata_ref, $new_name, $counters_ref, $gene) = @_;        
		my $read_length = length ($hash{$fq_file}{"seq"});
        chomp $$hash_ref{$fq_file}{"seq"};
		chomp $$hash_ref{$fq_file}{"qscore"};
		#print $$psldata_ref{$read_name}{$fq_file}{"basearray"}."\n";
		$$counters_ref{"02 No matches ".$$psldata_ref{$read_name}{"No_matches"}}++;
	if(exists ($$psldata_ref{$read_name}{$fq_file}{"basestring"}))
	{
		my @array = split /,/, $$psldata_ref{$read_name}{$fq_file}{"basestring"};	#Splits the string of BLAT matches into the array
		#print Dumper(\@array);
        my $qStart = 0;
		#print $read_name."\n".$$psldata_ref{$read_name}{$fq_file}{"basearray"}."\n".$$hash_ref{$fq_file}{"seq"}."\n".$$hash_ref{$fq_file}{"qscore"};
        my$j=0; #used to track the left hand side of the fragment
		for (my$i=0; $i+1<(scalar(@array)); $i++)
        {
			if ($array[$i] ne $array[$i+1]) # Checks if the next letter in the array is the same as the last
            {
            my $seqLength = $i+1 -$qStart;
			#print "i:$i\t$seqLength\n";
            if ($array[$i] =~/Y/) # This deals with the parts of the read that do not blat
				{
                my $read_type = "Type:R_Coord:$j-$i";
                if ($seqLength > $min_read_length){read_output2 ($fh_ref, $read_name ,$fq_file, $hash_ref, $psldata_ref, $new_name, $counters_ref, $read_type, $qStart, $seqLength, $gene);}
				}
				elsif (($repcapture == 1)) # Deals with the parts of the read that do blat
				{
                my $read_type = "Type:C$array[$i]_Coord:$j-$i";
                if ($seqLength > $min_read_length){read_output2 ($fh_ref, $read_name ,$fq_file, $hash_ref, $psldata_ref, $new_name, $counters_ref, $read_type, $qStart, $seqLength, $gene);}
				}
             else{$$counters_ref{"00 split read output matching error 1!"}++}
			 $qStart=$i+1;
			 $j=$i
            }
			if ($i+1 ==scalar(@array)-1) # This deals with the last part of the read
			{
				my $seqLength = $i+2 -$qStart;
				if ($array[$i] =~/Y/) 
				{
                my $read_type = "Type:R_Coord:$j-$i";
                if ($seqLength > $min_read_length){read_output2 ($fh_ref, $read_name ,$fq_file, $hash_ref, $psldata_ref, $new_name, $counters_ref, $read_type, $qStart, $seqLength, $gene);}
				}
				elsif (($array[$i] =~/(\d++)/) and ($repcapture == 1))
				{
                my $read_type = "Type:C$array[$i]_Coord:$j-$i";
				my $end = $qStart+$seqLength;
                #print "read length $read_length qstart $qStart seqLength $seqLength end \n";
				if ($seqLength > $min_read_length){read_output2 ($fh_ref, $read_name ,$fq_file, $hash_ref, $psldata_ref, $new_name, $counters_ref, $read_type, $qStart, $seqLength, $gene);}
				}
				else{$$counters_ref{"00 split read output matching error 2!"}++}
			}
        }
	}
	else{$$counters_ref{"00 basestring undefined - likely problem with Match string collapse subroutine"}++}

}

# Outputs the whole read as a single sequence
sub read_output2
{
    my ($fh_ref, $read_name ,$fq_file, $hash_ref, $psldata_ref, $new_name, $counters_ref, $read_type, $qStart, $seqLength, $gene) = @_;    
			my $seq = substr($$hash_ref{$fq_file}{"seq"},$qStart, $seqLength)."\n";
			my $qscore = substr($$hash_ref{$fq_file}{"qscore"},$qStart, $seqLength)."\n";
            
            #print "\nread out: $read_name $read_type\n$seq$qscore\n";
			#$$counters_ref{"05a $gene $read_type reads outputted to FASTQ:"}++;
			print $fh_ref $$hash_ref{$fq_file}{"name"}.$new_name.$read_type.$$hash_ref{$fq_file}{"name_end"}."\n".$seq.$$hash_ref{$fq_file}{"spare"}.$qscore; #Outputs the name of the captured read
	
			#$$counters_ref{"02 Read type $read_type:"}++;
}

sub read_output_dump
{
    my ($fh_ref, $read_name ,$fq_file, $hash_ref, $counters_ref) = @_;    
			my $seq = $$hash_ref{$fq_file}{"seq"};
			my $qscore = $$hash_ref{$fq_file}{"qscore"};
            
            #print "\nread out: $read_name $read_type\n$seq$qscore\n";
			$$counters_ref{"05a dump reads outputted to FASTQ:"}++;
			print $fh_ref $$hash_ref{$fq_file}{"name"}."\n".$seq.$$hash_ref{$fq_file}{"spare"}.$qscore; #Outputs the name of the captured read
	
			#$$counters_ref{"02 Read type $read_type:"}++;
}


# Trims the end off the names, modified for SRA data
sub name_trim
{      
	my ($name, $read, $hash_ref) = @_;
	if ($name =~ /(.*)( .*)( .*)/)
	{
	$$hash_ref{$read}{"name"}=$1;
	$$hash_ref{$read}{"name_end"}=$2;
	}

	elsif ($name =~ /(.*)( .*)/)
	{
	$$hash_ref{$read}{"name"}=$1;
	$$hash_ref{$read}{"name_end"}=$2;
	}
	elsif ($name =~ /(.*)/)
	{
	$$hash_ref{$read}{"name"}=$1;
	# print $1;
	$$hash_ref{$read}{"name_end"}="";
	}
	else {die "$name read name fails to parse"}
}

# This subroutine reads n lines from the psl file. It returns the flag 0 if it hasn't reached the end of the file and the flag 1 if  the end of the file is reached
sub psl_read
{	
my ($line_target) = @_; 
my $flag = 1;
	for (my$i=0; $i<=$#psl_files; $i++)
	{
		my $counter =0;
		while (my $line = readline $PSLIN[$i])
		{
			# fastq and psl files for reads must be put in in the same order
			$counters{"01 PSL lines total:"}++;
			if ($line =~ /^(\d++\t\d++)/)
			{
				my $current_read_name=line_read2($line, \%counters, $fastq_files[$i], \%psldata, \%gene_hash); $counters{"01 PSL lines read:"}++;
			
				# This part of the code analyses the groups of lines from the same read and collapses the arrays in the hash for the match strings on the fly to stop the code from using too much RAM
				if ($last_read_name eq $current_read_name){} #continues to read in lines with the same read name
				elsif ($last_read_name eq "first"){$last_read_name=$current_read_name} #deals with the first line 
				else #actions if the readname changes
				{
					Match_string_collapse(\%psldata, $last_read_name,\%counters);  #Runs the subroutine that collapses the match strings in the psldata file to save on RAM
					$counters{"01 PSL reads in memory"}++;
					#Psl_out (\%psldata, $last_read_name, \%counters, \*PSLTMP);
					$last_read_name=$current_read_name;
				} 
			}
			else{$counters{"01 PSL file lines not read:"}++}
			$counter++;
			if ($counter > $line_target){$flag=0;last}
		}
	}
	return $flag;
}


#Reads in the data from the psl into the hash 
sub line_read2
{
  my ($line, $counters_ref, $fq_file, $psldata_ref, $gene_hash_ref) = @_;
  # next if (/^psL/ || /^match/ || /^\s+match/ || /^----/ || /^$/);
  chomp $line;
  $$counters_ref{"01 PSL lines read:"}++;
	
  my ($match, $mismatches, $repmatch, $ncount, $qNumInsert, $aBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName, $qSize, $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount, $blockSizes, $qStarts, $tStarts) = split /\t/, $line; 
		# print "$qName\n";
  my $read_name = "@".$qName; #puts the @ back onto the readname
		# print "$tName\n";
	if ($tName =~ /(.*):(.*)-(.*)_(.*)/)
    {
	my $gene = $4;
	$$gene_hash_ref{$4}{"counts"}++;
	$$counters_ref{"01 $fq_file counts loaded from PSL:"}++;
	$$counters_ref{"03a gene $4 $fq_file counts loaded from PSL:"}++;
	$$psldata_ref{$read_name}{"No_matches"}++;
	if (exists($$psldata_ref{$read_name}{"best_match"}))# if the read exists already
		 {
				# print "$read_name\t$$psldata_ref{$read_name}\n";
			$$counters_ref{"03b $fq_file multimatch:"}++;
			$$counters_ref{"03b multimatch $gene v ".$$psldata_ref{$read_name}{"gene"}}++;
			if ($match > $$psldata_ref{$read_name}{"best_match"}) # if the match is the best one from all of the reads inputted
			{
			$$psldata_ref{$read_name}{"best_match"}=$match;
			$$psldata_ref{$read_name}{"gene"}=$gene;# specifies the gene / oligo match to use for the output file
			$$psldata_ref{$read_name}{"strand"}=$strand; #be careful of putting too much weight on this - the strand can change but this should be the strand of the best match
			}
			else {} # if the match isn't the best from all the psl files
		 }
		else # if this is the first time a read is loaded with that name
		{
			# print "$read_name\t$match\n";
		$$psldata_ref{$read_name}{"best_match"}=$match;
		$$psldata_ref{$read_name}{"gene"}=$gene;
		$$psldata_ref{$read_name}{"strand"}=$strand; #be careful of putting too much weight on this - the strand can change but this should be the strand of the best match
		}
	Match_string ($psldata_ref, $read_name, $fq_file, $strand, $qSize, $qStarts, $blockSizes, $counters_ref);
	$$psldata_ref{$read_name}{"string"} .= "|$fq_file-$gene-$match";
	
	}
    else{$$counters_ref{"00! ERROR PSL gene undefined:"}++}
	# exit;
	return $read_name;
}

# This makes an array of Ys / numbers depending on the block sizes in the PSL file - the numbers are the number of the match from the $$psldata_ref{$read_name}{"No_matches"}
# The code overwrites the string - a sequence would be as follows
# YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY - nothing blats
# YYYYY1111111111111111111YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY - single blat
# YYYYYYYYYYYYYYYYY22222222222222222222YYYYYYYYYYYYYYYYYY - second blat
# YYYYYYYYYYYYYY3333333333333333333333333333333YYYYYYYYYY - third blat
# YYYYYYYYYYY444444444444444444444444YY4444444444YYYYYYYY - fourth blat

sub Match_string
{
    my ($psldata_ref, $read_name, $fq_file, $strand, $qSize, $qStarts, $blockSizes, $counters_ref) = @_;
	unless (exists $$psldata_ref{$read_name}{$fq_file}{"basearray"}[0]){$$psldata_ref{$read_name}{$fq_file}{"basearray"}=[("Y") x $qSize]}; #Makes a string of Ys of the length of the read
	my @blockSizes = split /,/, $blockSizes; # Splits the block sizes based on the commas
    my @qStarts = split /,/, $qStarts; # Splits the qStarts based on the commas
	my $total_block_size=0;
	my $ref;
	my $no_matches=$$psldata_ref{$read_name}{"No_matches"};
	$$psldata_ref{$read_name}{$fq_file}{"psl basearray"}{$no_matches}=[("Y") x $qSize]; #Makes a new base array for each new BLAT match but not for each block
    if ($strand eq "-")  #Sorts out the strange notation for the reverse strand (which treats effectively annotates the reverse complement) - it is 0 based
    {
        for (my $i=0; $i <= $#qStarts; $i++) #Loops through each of the qStarts in the array #https://genome.ucsc.edu/FAQ/FAQformat.html#format2
        {
            $qStarts[$i] = $qSize-$qStarts[$i]-$blockSizes[$i];
        }
    }
	for (my $i=0; $i < scalar(@qStarts); $i++) #Changes the "Y"s to refs at the positions in the string defined by the block sizes / qstarts
    {
		$ref = $no_matches.".".$i;
		for (my $j = $qStarts[$i]; $j<($blockSizes[$i]+$qStarts[$i]); $j++){$$psldata_ref{$read_name}{$fq_file}{"psl basearray"}{$no_matches}[$j]=$ref};

		$total_block_size += $blockSizes[$i]; #
	}
	if(exists ($$psldata_ref{$read_name}{$fq_file}{"total_block_size"}{$total_block_size}))
	{
		while (exists($$psldata_ref{$read_name}{$fq_file}{"total_block_size"}{$total_block_size})){$total_block_size=$total_block_size - 0.1} #In the instance where there are two identical sized blocks the block size is subsequently reduced slightly ... will mean that the latter reads are collapsed first
	}
	$$psldata_ref{$read_name}{$fq_file}{"total_block_size"}{$total_block_size}="$no_matches"; #puts the total block size into the hash and points back to the read number for easy mapping later
}

# This makes a string of Ys / numbers depending on the block sizes in the PSL file - the numbers are the number of the match from the $$psldata_ref{$read_name}{"No_matches"}
# The code overwrites the string in the reverse order of the total size of the blocks - The sequence would be as follows
# YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY - nothing blats
# YYYYY111YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY - single blat
# YYYYY111YYYYYYYYY22222222222222222222YYYYYYYYYYYYYYYYYY - second blat
# YYYYY111YYYYYY3333333333333333333333333333333YYYYYYYYYY - third blat
# YYYYY111YYY444444444444444444444444334444444444YYYYYYYY - fourth blat
# $$psldata_ref{$read_name}{$fq_file}{"basearray"}{$char}

# It collapses the arrays into a string - to prevent the RAM from going completely nuts.
sub Match_string_collapse
{
    my ($psldata_ref, $read_name, $counters_ref) = @_;
		foreach my $fq_file (@fastq_files)  
		{
			if (exists($$psldata_ref{$read_name}{$fq_file}))
			{
				foreach my $total_block_size (sort {$a <=> $b} keys %{$$psldata_ref{$read_name}{$fq_file}{"total_block_size"}}) #Goes through the hash in ascending order of the total size of the blocks
				{
					my $counter=0;
					
					for (my $i=0; $i < scalar(@{$$psldata_ref{$read_name}{$fq_file}{"basearray"}}); $i++) #Loops through each of the qStarts in the array
						{
							#print TESTFH $$psldata_ref{$read_name}{$fq_file}{"basearray"}[$i]; #For debugging
							#$$psldata_ref{$read_name}{$fq_file}{"total_block_size"}{$total_block_size} produces the read number corresponding to the block size
							if ($$psldata_ref{$read_name}{$fq_file}{"psl basearray"}{$$psldata_ref{$read_name}{$fq_file}{"total_block_size"}{$total_block_size}}[$i] ne "Y")
								{
									$$psldata_ref{$read_name}{$fq_file}{"basearray"}[$i]=
									$$psldata_ref{$read_name}{$fq_file}{"psl basearray"}{$$psldata_ref{$read_name}{$fq_file}{"total_block_size"}{$total_block_size}}[$i]; #Pulls out the correct base array (the $char of the base array is stored in the block size)
									$counter++
								}
						}
				}
				# Converts array to string of comma separated values and deletes the arrays after they have been used to save on RAM requirements
				$$psldata_ref{$read_name}{$fq_file}{"basestring"} = join(",", @{$$psldata_ref{$read_name}{$fq_file}{"basearray"}});
				delete $$psldata_ref{$read_name}{$fq_file}{"basearray"};
				delete $$psldata_ref{$read_name}{$fq_file}{"psl basearray"};
			}
		}
}


#This ouputs a 2column hash to a file
sub output_hash
{
    my ($hashref, $filehandleout_ref) = @_;
    foreach my $value (sort keys %$hashref)
    {
    print $filehandleout_ref "$value\t".$$hashref{$value}."\n";
    }        
}
