#!/usr/bin/perl -w

use strict;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper; #Used for debugging

=head1 NAME

 MCC_analyser.pl

=head1 SYNOPSIS
 
 This script takes in the SAM file resulting from alignment with Bowtie 2 of the output of the MCC_splitter.pl script.
 For the purpose of this script we will use the following semantics. Following sequencing READ 1 and READ 2 are combined into a single READ. This is then split into SLICES depending on how the sequence maps to the 800 bp around oligo sequence or not. The MCC_splitter.pl script parses the output from BLAT and the FASTQ file to do this.
 
 It has a number of functions:
 1. It checks that the READS contain a SLICE that maps to the oligo sequence itself rather than just the surrounding region
 2. It collapses PCR duplicates with a wobble
 3. Derives the ligation junction from the position of the SLICES in the READ and their strand
 4. It derives whether the read is up or downstream of the junction - for a ChIP nexus like approach to transcription factor identification
 
 The script outputs wig files and a custom junction file 
 
=head1 EXAMPLE

perl MCC_analyser.pl -f Gene_name.sam -s MCC1

=head1 OPTIONS

 -f		Input filename 
 -o		Oligonucleotide position filename 
 -pf	Your public folder 
 -pu	Your public url 
 -limit		Limit the analysis to the first n reads of the file
 -genome	Specify the genome (mm9 / hg18)

=head1 AUTHOR

 Written by James Davies 2020

=cut

# Specifications
# Strings
# NB. The following strings should be changed to to your public account, server url and bigwig folder.
my $public_folder = "";
my $public_url = "";
my $bigwig_folder = "";
my $email = '';

my $input_path = '';
my $oligo_filename = "";
my $sample = "MCC";
my $analysis_read;
my $last_read="first";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); #defines the start time of the script (printed out into the report file later on)
my $use_dump =0; #whether to create an output file with all of the non-aligning sequences
my $use_limit=0; #whether to limit the script to analysing the first n lines
my $genome = "mm9";
my $help=0;
my $man=0;
my $target_chr='Undefined';
my %oligo_data;

# Arrays
my @samheadder; # contains the headder of the sam file

# Hashes
my %data; # contains the parsed data from the sam file
my %coords_hash=(); # contains a list of all the coordinates of the mapped reads to exclude duplicates
my %coords_c2n=(); # contains a hash of coordinates pointing to the readname
my %coords_n2c=(); # contains a hash of readnames pointing to the coordinates
my %counters;  # contains the data for all the counters in the script, which is outputted into the report file
my %duplicates;
my %wighash_blat_rep; #reporters based on typing by BLAT
my %wighash_blat_cap;
my %wighash_oligo_rep; #reporters based on typing by mapping by bowtie outside the capture oligo coordinates
my %wighash_oligo_cap;
my %wighashall_oligo_rep; #non-deduplicated reads based on mapping by bowtie
my %wighashall_oligo_cap;
my %junction;



# The GetOptions from the command line
&GetOptions
(
	"f=s"=>\ my $input_filename_path, 		# -f		Input filename 
	"o=s"=>\ $oligo_filename,			        # -o		Oligonucleotide position filename 
	"pf=s"=>\ $public_folder,		        	# -pf		Your public folder (e.g. /hts/data0/public/username)
	"pu=s"=>\ $public_url,			        	# -pu		Your public url (e.g. sara.molbiol.ox.ac.uk/public/username)
	"bf=s"=>\ $bigwig_folder,							# -bf 	Bigwig folder location
	"s=s"=>\ $sample,				            	# -s		Sample name (and the name of the folder it goes into)
	"dump"=>\ $use_dump,			        		# -dump		Print file of unaligned reads (sam format)
	"limit=i"=>\ $use_limit,		        	# -limit		Limit the analysis to the first n reads of the file
	"genome=s"=>\ $genome,			        	# -genome	Specify the genome (mm9 / hg18)
	'h|help'=>\$help,				            	# -h or -help 	Help - prints the manual
	'man'=>\$man,					            		# -man  	prints the manual
);

pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;
pod2usage(2) unless ($input_filename_path); 

my $full_url = $public_url.$public_folder;

#Splits out the filename from the path
my $input_filename=$input_filename_path;

if ($input_filename =~ /(.*\/)(\V++)/) {$input_path = $1; $input_filename = $2};

unless ($input_filename =~ /(^([^_]++)_.*).sam/) {die"filename does not match .sam"};

# Creates files for the a report and a fastq file for unaligned sequences to go into
my $outputfilename = $1; 
my $report_filename = $1."_report.txt";
my $coord_filename = $1."_coordstring.txt";
my $target_from_filename = $2; #Parses the target out of the filename;

unless (open(REPORTFH, ">$report_filename")){die "Cannot open file $report_filename $!\n"; exit};
unless (open(CHIMOUT, ">$outputfilename\_chimeric.txt")){die "Cannot open file $outputfilename\_chimeric.txt$!\n"; exit};
unless (open(SAMOUT, ">$outputfilename\_junction.sam")){die "Cannot open file $outputfilename\_chimeric.txt$!\n"; exit};

PRINTOUT:
# Prints the statistics into the report file
print REPORTFH "Script version: $0
Input parameters
Input path/file: $input_path$input_filename
Public folder: $public_folder
Public URL: $public_url
Full URL: $full_url
Bigwig Foler: $bigwig_folder
Genome $genome
Sample name: $sample
Oligo file: $oligo_filename
#### Script started at: $mday/$mon/$year $hour:$min:$sec                              ####
#### If you don't see a time stamp for the script finishing it has crashed somewhere! ####
";


# Uploads coordinates of capture oligos and exclusion regions into the array @oligo_data
open(OLIGOFH, $oligo_filename) or die "Cannot open oligo file $oligo_filename $!\n";
#Splits out the filename from the path
my $oligo_file=$oligo_filename; my $oligo_path;
if ($oligo_file =~ /(.*)\/(\V++)/) {$oligo_path = $1; $oligo_file = $2};
unless ($oligo_file =~ /(.*)\.(.*)/) {die"Cant regex oligo filename"};
my $oligo_bed_out_filename = $1.".bed";


while ( my $line = <OLIGOFH> )
{
  #print $line;
	chomp $line;
	$line =~ s/\r//; #in case carriage returns cause a problem
  if ($line =~ /^>chr(.*):(\d++)-(\d++)_(.*)/)
	{
		$counters{"01 Number of oligos loaded"}++;
	  chomp $line;
		my $gene = $4; my $chr = $1; my $start = $2; my $end = $3;
    $oligo_data{$gene}{"chr"}=$chr; 
    $oligo_data{$gene}{"start"}=$start;
    $oligo_data{$gene}{"end"}=$end;
		my $oligo_size = $end-$start;
		$counters{"10_input_oligo_size:\t$oligo_size"}++;
		if ($gene =~ /$target_from_filename/)
		{
			print REPORTFH "From the filename - processing $gene\tchr$chr\toligo start\t$start\toligo_end\t$end\n\n";
			$target_chr=$chr;
		}
  }
	elsif ($line =~ /[ACTGactg]/){$counters{"01 DNA sequences in oligo file qq"}++;}#Should disregard the DNA sequences in the oligo file
	else {print REPORTFH "!!Oligo file line failed to parse $line\n"}
};


# Opens input sam file 
unless (open(INFH, "$input_path$input_filename")) {print "Cannot open file $input_path$input_filename $!\n"; exit};

while (my $line = <INFH>)  #loops through all the lines of the sam file
{
    chomp $line;
    $counters{"01b Total lines in sam file:"}++;

    if ($line =~ /^(@.*)/) #removes headder of sam file and stores it in the array @samheadder
    {
    push @samheadder, $line;
    print SAMOUT $line."\n";
    $counters{"01a Lines in sam file headder:"}++; next
    } 
        
    my ($name, $bitwise, $chr, $readstart, $mapq, $cigar, $cola, $colb, $colc, $sequence, $qscore, $resta, $restb, $restc) = split /\t/, $line;   # splits the line of the sam file   
  

    if (($chr =~ /chr(.*)/) and ($name =~ /(.*)_Gene>(.*)([+-])Type:(.*)_Coord:(\d++)-(\d++)$/)) # checks that the sam file matches the chr and name (non aligned reads will not do this)
    {
        
		$counters{"01c Aligning sequences:"}++;
		
	# To pull the coordinates of the slice from the original read

			$name =~ /(.*)_Gene>(.*)([+-])Type:(.*)_Coord:(\d++)-(\d++)$/;
			my $readname=$1; my $gene = $2; my $blatstrand = $3;  my $readtype = $4; my $coord_start = $5; my $coord_end = $6;
      $gene =~ s/\r//; #Again to remove carriage returns from the gene name
			$counters{"10 Gene: $gene"}++;
			
			if ($gene ne $target_from_filename){$counters{"!! potential problem - gene doesn't match target from filename $gene v $target_from_filename"}++}
			
	# Puts the gene name into the hash and checks that there aren't multiple matches
            if (exists ($data{$readname}{"gene"}))
                {
                    if ($data{$readname}{"gene"} eq $gene){}
                    else {$counters{"multiple gene match"}++; $counters{"multiple gene match $gene"}++}
                }
            else{$data{$readname}{"gene"}=$gene}
        
    # Counts the number of fragments in each read
			$data{$readname}{"number of fragments"}++;
			my $readno = $data{$readname}{"number of fragments"} -1;
		
    # Assigns the entire line of the sam file to the hash
	       	$data{$readname}{$readno}{"whole line"}= $line;
            
	# Parses the chromosome from the $chr - nb without the chr in front
            $chr =~ /chr(.*)/; $chr = $1;
            $data{$readname}{$readno}{"chr"} = $1;
	    
	#This uses the cigar subroutine to interpret the cigar string
        my ($cigarflag, $readend, $readstrand, $seqlength) = cigar($bitwise, $cigar, $readstart);
				
			if ($cigarflag == 0)
            {
                $data{$readname}{$readno}{"seqlength"} = $seqlength;
                $data{$readname}{$readno}{"readstart"} = $readstart;
                $data{$readname}{$readno}{"readend"} = $readend;
				$data{$readname}{$readno}{"strand"} = $readstrand;
                
                #Generates a string of the coordinates of all the split paired end reads to allow duplicates to be excluded
                my $coordinate_string = $data{$readname}{$readno}{"chr"}.":".$data{$readname}{$readno}{"readstart"}."-".$data{$readname}{$readno}{"readend"};
				push (@{$data{$readname}{"coord array"}}, $coordinate_string); #coord array method for excluding precisely matching duplicates
				
				#Hashes of arrays required for new method allowing wobble
				push (@{$coords_n2c{$readname}}, $coordinate_string); #readname pointing to an array of coordinates
				push (@{$coords_c2n{$coordinate_string}}, $readname); #coordinates pointing to an array of readnames
            }
        else{
			$counters{"00 ERROR - number cigar strings failing to parse:"}++;
			$data{$readname}{$readno}{"strand"}=$blatstrand; #Shouldnt happen - think the strand shoudl come
			next
			}
            
    #Defines the read type based on the matching by BLAT
            $data{$readname}{$readno}{"blat_type"}=$readtype; #This is the readtype as defined by BLAT
            $counters{"21d readtype: $readtype"}++;
            $data{$readname}{$readno}{"coord_start"}=$coord_start; #This is the coordinates of the oligo given in the fa oligo 'genome' file used by BLAT
            $data{$readname}{$readno}{"coord_end"}=$coord_end; #This is the coordinates of the oligo given in the fa oligo 'genome' file used by BLAT
            $data{$readname}{$readno}{"blatstrand"}=$blatstrand;
            $data{$readname}{"start hash"}{$coord_start}=$readno; #hash to pull out the read with a given start site
						$data{$readname}{"end hash"}{$coord_end}=$readno;
        
	#This part of the code defines whether the read is a capture or reporter read or whether it is proximity excluded
            my $flag=0;
	    
	#Loops through the @oligo_data array to see whether the fragment meets the criteria to be a capture or exclusion fragment
		if (exists $data{$readname}{"oligo"}){} else {$data{$readname}{"oligo"}="N"} #populates these flags if they haven't been declared
		if (exists $data{$readname}{"reporter"}){} else {$data{$readname}{"reporter"}="N"}
		
		$data{$readname}{$readno}{"type"}="R"; #Defines the type as a reporter for it to be overwritten if there is a capture
		$counters{"02 readno reporters defined"}++;
		
    foreach my $oligo(keys(%oligo_data))
    {
     # Defines if the fragment lies within the capture region = NB this version of the script requires the part of the read to hit the oligo..
      if (($data{$readname}{$readno}{"chr"} eq $oligo_data{$oligo}{"chr"}) and ($data{$readname}{$readno}{"readstart"}<$oligo_data{$oligo}{"end"} and $data{$readname}{$readno}{"readend"}>$oligo_data{$oligo}{"start"})) 
      {
		    if (exists$data{$readname}{"oligo name"}) #Checks if there is a read that has already been assigned as an oligo
					{
						if ($data{$readname}{"oligo name"} ne $oligo) #to exclude reads matching to more than one of the oligos !(might be dangerous in some circumstances with two identical targets).
						{
						 $counters{"02!_multi_match_target:"}++;
						 $counters{"22_oligo_multimatch_$oligo:"}++;
						 $data{$readname}{$readno}{"type"}="E"; #type set to E for exclude -  could be dangerous if you have two identical oligos - in which case change to ="O"
						 $counters{"02 readno reporters defined"}--;
						}
						else
						{
							$data{$readname}{$readno}{"type"}="P"; #type set to P for reads that are second Oligo matches
							$counters{"02_oligo_$oligo:"}++;
							$counters{"02_oligos_assigned:"}++;
							$counters{"02 readno reporters defined"}--;	
						}
					}
					else
					{
			if ($oligo eq $gene) #Checks that the oligo is the same as the blat match
						{ 
							$data{$readname}{$readno}{"type"}="O"; #type set to O for Oligo match
							$data{$readname}{"oligo name"}=$oligo; #Flag to show that there is an oligo matching fragment in the read
							$data{$readname}{"oligo"}="Y";
							$counters{"02_oligo_$oligo:"}++;
							$counters{"02_oligos_assigned:"}++;
							$counters{"02 readno reporters defined"}--;
						}
						else
						{
							$data{$readname}{$readno}{"type"}="E"; #Exclude this read
							$counters{"02!Blat vs oligo mismatches"}++;
							$counters{"26!Blat vs oligo mismatches $oligo"}++;
							$counters{"02 readno reporters defined"}--;
						}
					}
		  }
			
    }
		# Checks if the type has been assigned (this should only happen if it maps to the capture oligo region) Otherwise it should be defined as a reporter.
		if ($data{$readname}{$readno}{"type"} !~ /[OE]/)
		{
			$data{$readname}{"reporter"}="Y";
			$counters{"02 readname reporters defined"}++;
		}

		
        # Checks whether the name of the read has changed compared to the last one.  If there is no change it continues to loop through all the fragments of the read until it changes.
        # If the readname has changed it loads the data into the output hashes and deletes the lines from the %data hash.
        # NB this has the property of throwing data if the reads are not in order in the sam file... samtools sort -n -o file_sorted.sam file.sam 
        if ($last_read eq "first"){$last_read=$readname} #deals with the first read
        elsif($readname eq $last_read){} #checks whether the read name has changed and moves on to the next line in the sam file if the read names are the same
        else
        {
			$analysis_read = $last_read;
			my $read_reporter = $data{$analysis_read}{"reporter"};
			my $oligo_reporter = $data{$analysis_read}{"oligo"};
			
			$counters{"01a Total number of reads:"}++;
			$counters{"11 Non PCR duplicate filtered reads types reporter v oligo $read_reporter v $oligo_reporter:"}++;

					
        #This line sorts the array in the hash and converts it into a string with the sequences in ascending order 
        $data{$analysis_read}{"coord string"}= join( "_", sort {$a cmp $b} @{$data{$analysis_read}{"coord array"}}); 
		
		# This checks if the coordinate string is unique and uses this to collapse PCR duplicates - this is the first line of PCR duplicate filtering based on all coordinates being unique
		if (exists $coords_hash{$data{$analysis_read}{"coord string"}}) 
        {
    			$coords_hash{$data{$analysis_read}{"coord string"}}++;
    			$counters{"01a Duplicate reads Total:"}++;
    			$counters{"06c Duplicate reads n=".$data{$analysis_read}{"number of fragments"}.":"}++;
        } #Checks if the read coordinates have been seen before
					
        else
		{
			# Adds the coordinates string to the duplicates hash
			$coords_hash{$data{$analysis_read}{"coord string"}}++; 
			
			# Tests if there is duplication using a wobble in the start / stop coordinates of any of the slices
			if (pcr_duplicate_finder($analysis_read)==0) 
			{
			########################################################################################################################
			#Below here the reads should be PCR duplicate filtered
			########################################################################################################################
          
			$counters{"01a Number of non-duplicated reads Total:"}++;
			$counters{"11 PCR duplicate filtered reads types reporter v oligo $read_reporter v $oligo_reporter:"}++;
			$counters{"06b Number of non-duplicated reads n=".$data{$analysis_read}{"number of fragments"}}++;
		  
			#print Dumper(%data);
			
		  for (my $readno = 0; $readno< $data{$analysis_read}{"number of fragments"}; $readno++)
              {
        ########################################################################################################################
				#This part of the script looks for chimeric junctions - i.e. reads with a perfect junction (+/- an offset of 2bp) in the read and it uses between the oligo mapping reads and other reads for footprinting
				########################################################################################################################
					my $flag=0;
	          if ($data{$analysis_read}{$readno}{"type"} eq "O")
						{
							for (my $matching_readno = 0; $matching_readno< $data{$analysis_read}{"number of fragments"}; $matching_readno++)
							{									
							my $allowed_gap =2;
							my $start_offset = abs($data{$analysis_read}{$readno}{"coord_end"}-$data{$analysis_read}{$matching_readno}{"coord_start"});
							my $end_offset = abs($data{$analysis_read}{$readno}{"coord_start"}-$data{$analysis_read}{$matching_readno}{"coord_end"});
							
							if ($matching_readno == $readno){} #ignore the oligo matching read
							elsif($flag>0){} #only output one junction per read
							#### Maps the reads with the oligo match on the Left of the reporter in the original read
				      elsif ($start_offset <=$allowed_gap)
								{
								$flag=1;
								#Parameters for junction_out($readname, $chrL, $startL, $endL, $strandL, $typeL, $read_coordL,$chrR, $startR, $endR, $strandR, $typeR, $read_coordR, $hash_ref)
								junction_out (
											  $analysis_read,
											  $data{$analysis_read}{$readno}{"chr"},
											  $data{$analysis_read}{$readno}{"readstart"},
											  $data{$analysis_read}{$readno}{"readend"},
											  $data{$analysis_read}{$readno}{"strand"},
											  $data{$analysis_read}{$readno}{"type"},
											  $data{$analysis_read}{$readno}{"coord_start"},
												$data{$analysis_read}{$readno}{"coord_end"},
											  $data{$analysis_read}{$matching_readno}{"chr"},
											  $data{$analysis_read}{$matching_readno}{"readstart"},
											  $data{$analysis_read}{$matching_readno}{"readend"},
											  $data{$analysis_read}{$matching_readno}{"strand"},
											  $data{$analysis_read}{$matching_readno}{"type"},
											  $data{$analysis_read}{$matching_readno}{"coord_start"},
												$data{$analysis_read}{$matching_readno}{"coord_end"},
												"Start_offset",
												\%junction);
								print SAMOUT $data{$analysis_read}{$readno}{"whole line"}."\n";
								print SAMOUT $data{$analysis_read}{$matching_readno}{"whole line"}."\n";
								}
							#### Maps the reads with the oligo match on the right of the reporter in the original read
							elsif ($end_offset <=$allowed_gap)
								{
								$flag=1;
								junction_out (
											  $analysis_read,
											  $data{$analysis_read}{$matching_readno}{"chr"},
											  $data{$analysis_read}{$matching_readno}{"readstart"},
											  $data{$analysis_read}{$matching_readno}{"readend"},
											  $data{$analysis_read}{$matching_readno}{"strand"},
											  $data{$analysis_read}{$matching_readno}{"type"},
											  $data{$analysis_read}{$matching_readno}{"coord_start"},
												$data{$analysis_read}{$matching_readno}{"coord_end"},
											  $data{$analysis_read}{$readno}{"chr"},
											  $data{$analysis_read}{$readno}{"readstart"},
											  $data{$analysis_read}{$readno}{"readend"},
											  $data{$analysis_read}{$readno}{"strand"},
											  $data{$analysis_read}{$readno}{"type"},
											  $data{$analysis_read}{$readno}{"coord_start"},
												$data{$analysis_read}{$readno}{"coord_end"},
												"End_offset",
												\%junction);
								print SAMOUT $data{$analysis_read}{$readno}{"whole line"}."\n";
								print SAMOUT $data{$analysis_read}{$matching_readno}{"whole line"}."\n";
								}
							}
						if($flag==0){$counters{"03 non perfect junctions"}++}
						}		

										
										
						########################################################################################################################
						# This part of the script outputs reads depending on whether part of the sequence maps inside the coordinates of the capture oligo i.e. 120bp
						# This is my prefered method of analysis at the moment - it has additional safeguards - has to have both oligo mapping and reporter
						# It should also provide additional detail around the capture site
						#######################################################################################################################
						
							if (($data{$analysis_read}{$readno}{"type"} eq "R") and ($data{$analysis_read}{"oligo"} eq "Y"))
                            {
                             $counters{"02_unique_reporters_with_oligo_match:"}++;  #I think this is the best metric for the quality of the data
							 if ($data{$analysis_read}{$readno}{"chr"} eq $target_chr){$counters{"02_unique_reporters_with_oligo_match_cis:"}++}
							 else{$counters{"02_unique_reporters_with_oligo_match_cis:"}++}
                             wighash_read ($data{$analysis_read}{$readno}{"chr"}, $data{$analysis_read}{$readno}{"readstart"}, $data{$analysis_read}{$readno}{"readend"}, \%wighash_oligo_rep, \%counters );
                            }
                            elsif (($data{$analysis_read}{$readno}{"type"} eq "O") and ($data{$analysis_read}{"reporter"} eq "Y"))
                            {
                             $counters{"02_unique_oligo_matches_with_reporter:"}++;
                             wighash_read ($data{$analysis_read}{$readno}{"chr"}, $data{$analysis_read}{$readno}{"readstart"}, $data{$analysis_read}{$readno}{"readend"}, \%wighash_oligo_cap, \%counters );
                            }
							else
							{
							$counters{"02_fragments_without_oligo_map_and_reporter"}++;
							my $oligo_type = $data{$analysis_read}{$readno}{"type"};
							my $blat_type = $data{$analysis_read}{$readno}{"blat_type"};

							$counters{"22a_fragments_without_oligo_map_and_reporter_oligo_match_$oligo_type:"}++;
							$counters{"22b_fragments_without_oligo_map_and_reporter_blat_match_$blat_type:"}++;
							$counters{"22c_fragments_without_oligo_map_and_reporter_oligo_v_blat_match_$oligo_type v $blat_type:"}++;
							$counters{"22c_fragments_without_oligo_map_and_reporter_oligo_v_blat_match_$oligo_type v $blat_type:"}++;
							$counters{"22d_fragments_without_oligo_map_and_reporter_reporter_v_oligo_readtype_$read_reporter v $oligo_reporter:"}++;
							}
				}
			}
			else{$counters{"02 Wobble duplicates excluded:"}++} #Should get to this part if they are specifically wobble duplicates - hash this out if using simple filtering
		}
		if ($counters{"01c Aligning sequences:"}>10){delete $data{$analysis_read}}  #deletes the slice in the hash with the data that has been removed - the first 10 lines are kept so that they can be visualised for debugging by Data Dumper
		$last_read = $readname; 
        }	      
  }
  else
  {
    #if ($use_dump eq 1){print DUMPOUTPUT $line."\n"} #"$name\n$sequence\n+\n$qscore\n"};
    $counters{"01c Non aligning sequences:"}++; next
  }

LOOP_END:    
unless ($use_limit ==0){if ($counters{"01c Aligning sequences:"}>$use_limit){last}}   # This limits the script to the first n lines
}


###############################################################################################
#Output of deduplicated data
my $rep_hash_ref = \%wighash_oligo_rep;
my $cap_hash_ref = \%wighash_oligo_cap;

#Output of normalised data
my $norm_factor = 1;
my $norm_filename = $outputfilename."_de_norm_rep";
if (exists ($counters{"02_unique_reporters_with_oligo_match:"}))
{
	if ($counters{"02_unique_reporters_with_oligo_match:"}==0){$norm_filename .= "WARNING_NOT_NORM_ZERO"}
	else {$norm_factor = 10000/$counters{"02_unique_reporters_with_oligo_match:"}}
}
else{$norm_filename .= "WARNING_NOT_NORM_UNDEF"}

$counters{"01c Normalisation factor"}=$norm_factor;
# print Dumper(\$rep_hash_ref);
# print Dumper(\%{$junction{"ALL"}});
wigout_normalised($rep_hash_ref, $norm_filename, $norm_factor);


#Output of foot printing - this is not normalised at this stage so that the data can be combined across replicates and normalised later
wigout_normalised(\%{$junction{"ALL"}}, $outputfilename."_ALL_FP", 1);
wigout_normalised(\%{$junction{"UP"}}, $outputfilename."_UP_FP", 1);
wigout_normalised(\%{$junction{"DO"}}, $outputfilename."_DO_FP", -1);

# wigout_plusone(\%{$junction{"ALL"}}, $outputfilename."_ALL_FP_plusone", 1);
# wigout_plusone(\%{$junction{"UP"}}, $outputfilename."_UP_FP_plusone", 1);
# wigout_plusone(\%{$junction{"DO"}}, $outputfilename."_DO_FP_plusone", -1);


#Output of the report file data.
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
print REPORTFH
"#### If you don't see a this time stamp for the script finishing it has crashed!      ####
#### Script finished at: $mday/$mon/$year $hour:$min:                                 ####
\nRun Statistics\n";
output_hash_counters(\%counters, \*REPORTFH);

exit;



#########################################################################################################################################################################
#Subroutines
#########################################################################################################################################################################
sub junction_out
{
	my ($analysis_read, $chrL, $startL, $endL, $strandL, $typeL, $coord_startL, $coord_endL,$chrR, $startR, $endR, $strandR, $typeR, $coord_startR, $coord_endR, $offset_type, $hash_ref) = @_;
	$counters{"03 Unique junctions identified:"}++;
							
	my ($genome_junctionL, $genome_readend1)=junction($startL, $endL, $strandL, "L");						
	my ($genome_junctionR, $genome_readend2)=junction($startR, $endR, $strandR, "R");	
                
	$$hash_ref{"ALL"}{$chrL}{$genome_junctionL}++;
	$$hash_ref{"ALL"}{$chrR}{$genome_junctionR}++; 
	
	$$hash_ref{"L"}{$chrL}{$genome_junctionL}++;
	$$hash_ref{"R"}{$chrR}{$genome_junctionR}++;
	
	# This works out if the junction is up or downstream of the read
	#     JunctionRRRRRRRRRRRRRRRRRR  = UP    vs RRRRRRRRRRRRRRRRRJunction = DOwnstream
	my $streamL;
	my $streamR;
	
	
	if ($strandL eq "+"){wighash_plusone($chrL, $genome_junctionL, \%{$$hash_ref{"DO"}}); $streamL='DO'}
	else{wighash_plusone($chrL, $genome_junctionL, \%{$$hash_ref{"UP"}}); $streamL='UP'}
	if ($strandR eq "+"){wighash_plusone($chrR, $genome_junctionR, \%{$$hash_ref{"UP"}}); $streamR='UP'}
	else{wighash_plusone($chrR, $genome_junctionR, \%{$$hash_ref{"DO"}}); $streamR='DO'}
	
																
	if ($chrL eq $chrR){$counters{"01c cis junctions"}++;}
	else{$counters{"01c trans junctions"}++;}
	#For full sanity check of junction callling
	#print CHIMOUT "chr$chrL\t$genome_junctionL\t$strandL\t$streamL\tchr$chrR\t$genome_junctionR\t$strandR\t$streamR\tLeft:$startL-$endL|$strandL->$genome_junctionL\tRight:$startR-$endR|$strandR->$genome_junctionR\t|RType-L$typeL-R$typeR|$coord_startL-$coord_endL,$coord_startR-$coord_endR|$offset_type|ReadNumber-$analysis_read\n";
	print CHIMOUT "chr$chrL\t$genome_junctionL\t$strandL\t$streamL\tchr$chrR\t$genome_junctionR\t$strandR\t$streamR\tLeft_$startL-$endL$strandL|Right_$startR-$endR$strandR|Type_$typeL-$typeR|Coord_$coord_startL-$coord_endL,$coord_startR-$coord_endR|$offset_type\t$analysis_read\n";
}


sub junction
{
    my ($genomeL, $genomeR, $strand, $read_pos) = @_;
    if ($read_pos eq "L")
    {
        if ($strand eq "+"){return ($genomeR, $genomeL)}
        elsif ($strand eq "-"){return ($genomeL, $genomeR)}
        else {die "Strand interpretation error $strand....."}
    }
    elsif ($read_pos eq "R")
    {
        if ($strand eq "-"){return ($genomeR, $genomeL)}
        elsif ($strand eq "+"){return ($genomeL, $genomeR)}
        else {die "Strand interpretation error"}
    }
    else {die "Read position interpretation error"}
}

sub wighash
{
my ($chr, $start, $hash_ref) = @_;
$$hash_ref{$chr}{$start}++;
}

sub wighash_plusone
{
my ($chr, $start, $hash_ref) = @_;
$$hash_ref{$chr}{$start}++;
$$hash_ref{$chr}{$start-1}++;
$$hash_ref{$chr}{$start+1}++;
}

sub wighash_read
{
my ($chr, $start, $end, $read_hash_ref, $counters_ref) = @_;
for(my $i=$start; $i <=$end; $i++)
		{
        $$read_hash_ref{$chr}{$i}++;
		}					
}


# This subroutine outputs the data from a hash of dpnII fragments of format %hash{$fragtype}{$capture}{$chr}{$fragment_start}{"end"/"value"} to wig format
sub wigout_normalised
{
    my ($hashref, $filenameout, $normalisation_factor) = @_;
    unless (open(WIGOUTPUT, ">$filenameout.wig")){print "Cannot open file $filenameout.wig\n";exit;}
    foreach my $storedChr (sort keys %$hashref)  
    {
    print WIGOUTPUT "variableStep  chrom=chr$storedChr\n";
    
    foreach my $storedposition (sort {$a cmp $b} keys %{$$hashref{$storedChr}})
        {
            my $output = $$hashref{$storedChr}{$storedposition}*$normalisation_factor;
            print WIGOUTPUT "$storedposition\t".$output."\n";
            
        }
    }
    wigtobigwig($filenameout, \*REPORTFH, $public_folder, $full_url, $filenameout, $bigwig_folder );
}

# This subroutine outputs the footprinting data to wig format with -1/+1 bp
sub wigout_plusone 
{
    my ($hashref, $filenameout, $normalisation_factor) = @_;
		#print "Normalisation factor = $normalisation_factor";
    unless (open(WIGOUTPUT, ">$filenameout.wig")){print "Cannot open file $filenameout.wig\n";exit;}

	foreach my $storedChr (sort keys %$hashref)  
    {
    print WIGOUTPUT "variableStep  chrom=chr$storedChr\n";
    foreach my $storedposition (sort {$a <=> $b} keys %{$$hashref{$storedChr}})
        {   
            my $window_sum=$$hashref{$storedChr}{$storedposition};
			if (exists ($$hashref{$storedChr}{$storedposition-1})){$window_sum += $$hashref{$storedChr}{$storedposition-1}};
			if (exists ($$hashref{$storedChr}{$storedposition+1})){$window_sum += $$hashref{$storedChr}{$storedposition+1}};
			my $output = $window_sum*$normalisation_factor;
            print WIGOUTPUT "$storedposition\t".$output."\n"; 
        }
    }
    #Hash this out if you are having problems with converting from wig to bigwig
		wigtobigwig($filenameout, \*REPORTFH, $public_folder, $full_url, $filenameout, $bigwig_folder );
}


#wigtobigwig (Filename (filename without the .wig extension), Filename for report (full filename), Description for UCSC)
#This subroutine converts a wig file into big wig format and copies it into the public folder.  It also generates a small file with the line to paste into UCSC
#if this fails to work try running: module add ucsctools  at the unix command line
sub wigtobigwig 
{
    my ($filename, $report_filehandle, $public_folder, $full_url, $description, $bigwig_folder) = @_;

    if (-e "$filename.wig")
		{
			if (-d "$public_folder/")
			{
				# if (-e "$bigwig_folder/$genome\_sizes.txt")
				if (-e "$bigwig_folder")
				{
					system ("module load ucsc\nwigToBigWig -clip $filename.wig $bigwig_folder $filename.bw") == 0 or print STDERR "$filename couldn't bigwig $genome files\n";
					system ("mv $filename.bw $public_folder/") == 0 or print $report_filehandle "$filename couldn't move files\n";
					print $report_filehandle "track type=bigWig name=\"$filename\" description=\"$description\" bigDataUrl=http://$full_url/$filename.bw\n";
				}
				else{print STDERR "Unable to bigwig $filename.wig - $bigwig_folder/$genome\_sizes.txt does not exist"}
			}
			else{print STDERR "Unable to bigwig $filename.wig - public folder does not exist"}
		}
		else{print STDERR "Unable to bigwig $filename.wig - file does not exist"}
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

#This ouputs a 2column hash to a file - with the numbers in a 1000,000.00 format
sub output_hash_counters
{
    my ($hashref, $filehandleout_ref) = @_;
    foreach my $value (sort keys %$hashref)
    {
		my $csv = $$hashref{$value};
		my $fraction = $$hashref{$value} - int($$hashref{$value});
		if (($fraction ==0) or ($csv > 10000)) #Will print the number as comma separated and remove the decimal places
		{
			$csv = int($$hashref{$value});
			$csv = reverse $csv;
			my @csv = unpack("(A3)*", $csv);
			$csv = reverse join ',', @csv;
		}
    print $filehandleout_ref "$value\t".$csv."\n";
    }        
}


#This ouputs a 2column hash to a file
sub output_hash_num
{
    my ($hashref, $filehandleout_ref) = @_;
    foreach my $value (sort {$b<=>$a} keys %$hashref)
    {
    foreach my $key (sort keys %{$$hashref{$value}})
			{
				print $filehandleout_ref "$value\t".$key."\n";
			}	
    }        
}

#This reads the cigar string from bowtie and outputs the start stop coordinates
sub cigar
{
    my ($bitwise, $cigar, $readstart) = @_;
		my $strand=0;
		my $seqlength=0;
    my $readlength=0;
    my $start_clip=0;
    my $end_clip=0;
    my @cigar=();
    my $flag = 0;  
        #Assigns the strand using the bitwise code    
    if ($bitwise==0){$strand="+"}
    elsif($bitwise==16){$strand="-"}
    else{$flag = 1;$counters{"00 ERROR - bitwise non 0,16 parse - STRAND UNDEFINED:"}++}
              
        #This calculates the end of the read using the cigar - complex matches are discarded by this.
      
				# Populates the array with the different sections of the cigar
        while ($cigar !~ /^$/)
        {
        if ($cigar =~ /^([0-9]+[MIDS])/){ push @cigar, $1;} #puts each part of the cigar into the array
        else {$counters{"Unexpected cigar:$cigar"}++;}
        $cigar =~ s/$1//; #removes the part of the cigar from the cigar
        }
        
        if ($strand eq "-"){@cigar=reverse(@cigar)} #reverses the cigar if the strand is negative
        
				# Runs through the array and calculates the
				# $seqlength is the length of the sequence that maps in the genome
				# $readlength is the length of the sequence in the read itself 
				# $start_clip is the length of the clip from the left side of the read
				# $end_clip is the length of the clip from the right side of the read
        for (my $i=0; $i<scalar(@cigar); $i++)
        {
            my $cigar_part = $cigar[$i];
            if ($cigar_part =~ /(\d+)M/){$seqlength += $1;$readlength += $1;}
            elsif ($cigar_part =~ /(\d+)I/){$readlength += $1;} #seqlength end doesn't change if there's an insertion but the readlength does
            elsif ($cigar_part =~ /(\d+)D/){$seqlength += $1;} #deletions - the seqlength changes but not the readlength
            elsif ($cigar_part =~ /(\d+)X/){$seqlength += $1;$readlength += $1;} #Substitutions both increase
            elsif ($cigar_part =~ /(\d+)N/){$seqlength += $1;} #Deals with introns - this hasn't been tested - might be wrong
            elsif ($cigar_part =~ /(\d+)S/){if ($i==0){$start_clip += $1;}else{$end_clip +=$1}} #only moves the start coordinate if it is the first element in the array
            elsif ($cigar_part =~ /(\d+)H/){if ($i==0){$start_clip += $1;}else{$end_clip +=$1}}
            elsif ($cigar_part =~ /(\d+)(.*)/){$counters{"Didn't parse cigar: $2"}++;} #print TESTOUT "Non parsed cigar!! $2"}
        }
        
      
				my $readend = $readstart+$seqlength-1; #minus 1 due to 1 based sam file input. $end position in the genome
								
                # It is possible to track the start and end positions in the read rather than the genome
								#my $seqStart = $seqStart+$start_clip; #start position in the original read
                #my $seqEnd = $seqStart+$start_clip+$readlength-1; ##end position in the original read I don't think you need to subtract the end_clip 
			return ($flag, $readend, $strand, $seqlength);
}


# This works by having two searchable hashes
# %coords_n2c is a hash of all read names -> chr:start-stop
# %coords_c2n is a hash of all chr:start-stop -> read names
# If there is a perfect match for any part of the read in the c2n hash is searches all of these to see whether they match the first read within the wobble of 4bp
# If there is a wobble match this is excluded

sub pcr_duplicate_finder
{
	my ($readname) = @_;
	my @reads_to_exclude;
	my $flag =0;
	# %coords_c2n # contains a hash of arrays with the coordinates pointing to the readnames
	# %coords_n2c # contains a hash of arrays of the readnames pointing to the coordinates
	foreach my $coordinate (@{$coords_n2c{$readname}})
	{
		foreach my $read(@{$coords_c2n{$coordinate}})
		{
			push @reads_to_exclude, $read
		}
	}

	foreach my $read(@reads_to_exclude)
	{
		if ($read eq $readname){}
		else
		{
		$flag =0;
			foreach my $readname_coords (@{$coords_n2c{$readname}})
			{
				foreach my $target_coords (@{$coords_n2c{$read}})
				{
					my $pcr_test= pcr_test($readname_coords, $target_coords);
					if ($pcr_test eq "error"){$counters{"Wobble duplicate filter error"}++}
					elsif($pcr_test ==0){} #could change to $flag-- to deal with triplicates etc with some unique reads and duplicates simultaneously
					elsif($pcr_test ==1){$flag++}
				}
			}
		if ($flag >=2){return 1} #if there are more than 2 fragments that match mark as PCR duplicate (will exclude some triplets unnecessarily)
		}
	}
	if ($flag <2){return 0}
	else {return 1}
}

# Tests if two coordinates are within the wobble in the start and stop
sub pcr_test
{
	my ($coord1, $coord2) = @_;
	my $chr1; my $chr2; my $start1; my $start2; my $end1; my $end2;
	my $wobble = 4;
	if ($coord1 =~ /(.*):(\d++)-(\d++)/){$chr1=$1; $start1=$2; $end1=$3} else{return "error"}
	if ($coord2 =~ /(.*):(\d++)-(\d++)/){$chr2=$1; $start2=$2; $end2=$3} else{return "error"}
	if ($chr1 eq $chr2)
	{
		my $start_dif = abs ($start1-$start2);
		if ($start_dif <= $wobble)
		{
			my $end_dif = abs ($end1-$end2);
			if	($end_dif <= $wobble){return 1;}
		}
		else {return 0}
	}
	else {return 0}
}


# James Davies 2020