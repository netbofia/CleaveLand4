#!/usr/bin/env bash 

#####
#Created by: Bruno Costa
# 28/05/2016 
# Copyright ITQB 2016
#
# This is a wrapper for Cleaveland 4
# it allows a list of sRNA to be divided into smaller lists
# once a list has been compleated a new list is created maintaining a 
# set of n thread allways running until the query has finished

# call: Cleaveland4-TH.sh [sRNA.fasta] [density/degradome] [transcriptome] [#Threads]

# Starting with lists of 10 sRNAs
# Using an already generated density map
# Transcriptome is allways the same | Beware it will be put into memory N time
set -e

#Path to cleaveland (Hardcoded Sorry will be changed to be incorporated to path instead)
cleaveland=~/.Software/CleaveLand4-4.3/CleaveLand4.pl

# OUTPUT-COLORING
red='\e[0;31m'
blue='\e[0;34m'
green='\e[0;32m'
blink='\e[5m'
unblink='\e[25m'
invert='\e[7m'
NC='\e[0m' # No Color


while [[ $# > 0 ]]
do
  key="$1"

case $key in
  -v|--version)
  version=true
  shift # past argument
  ;;
  -q|--quiet)
  quiet=true
  shift # past argument
  ;;
  -a|--allenSort)
  allenSort=true
  shift # past argument
  ;;
  -t|--tabular)
  tabular=true
  shift #past argument
  ;;
  -r|--mfer)
  mfer="$2"
  shift #past argument
  ;;
  -o|--tplots)
  tplots="$2"
  shift # past argument
  ;;
  -d|--degDensity)
  degDensity="$2"
  shift # past argument
  ;;
  -e|--degradome)
  degradome="$2"
  shift # past argument
  ;;
  -g|--gstar)
  gstar="$2"
  shift # past argument
  ;;
  -u|--sRNA)
  sRNA_COMPLETE="$2"
  shift # past argument
  ;;
  -n|--transcriptome)
  TRANSCRIPT="$2"
  shift # past argument
  ;;
  -p|--pval)
  tplots="$2"
  shift # past argument
  ;;
  -c|--category)
  tplots="$2"
  shift # past argument
  ;;
  -s|--threads)
  THREADS="$2"
  shift # past argument
  ;;
  -h|--help)
  echo -e " 
  ${blue}-f|--lib-first
  -l|--lib-last
  -h|--help${NC}
  ---------------------
  Optional args
  ---------------------
  ${blue}-s|--step${NC} Step is an optional argument used to jump steps to start the analysis from a different point  
      ${green}Step 1${NC}: Wbench Filter
      ${green}Step 2${NC}: Filter Genome & mirbase
      ${green}Step 3${NC}: Tasi
      ${green}Step 4${NC}: Mircat
      ${green}Step 5${NC}: PAREsnip    
 ${blue}--lc${NC} Set the program to begin in lcmode instead of fs mode. The preceading substring from the lib num (Pattern) Template + Lib num, but identify only one file in the inserts_dir    
 ${blue}--fasta${NC} Set the program to start using fasta files. As an argument supply the file name that identifies the series to be used. Ex: Lib_1.fa, Lib_2.fa, .. --> argument should be Lib_
 ${blue}--fastq${NC} Set the program to start using fastq files. As an argument supply the file name that identifies the series to be used. Ex: Lib_1.fq, Lib_2.fq, .. --> argument should be Lib_ , if no .fq file is present but instead a .fastq.gz file will additionally be extracted automatically.  
  "
  exit 0
esac
shift # past argument or value
done

#Clearing forbidden options 
if [[ -n "$degradome" && -n "$degDensity" ]]; then
  echo -e "${red}Error${NC} - you can't specifiy a degradome file and a degradome density. Please submit only one."
  exit 126
fi
if [[ "$version" ]]; then
  $cleaveland -v
  exit 0       
fi

#Planed support only for modes 1 and 2 
#Modes 3 and 4 use the existing GSTAr file (Not sure it is a use case for this abroad shouldn't need paralization)
#Support for the -p,-r,-v,-q,-a,-t,-o,-g and -c options will become a reality at least most but might take some time.
#Some are simple others involve some hacking around with conditions.


#Check it is higher than the numer of threads on the machine

#####################################################################################
#Get dir
DIR=$(pwd)
#Start time for benchmarking
START=$(date +%s.%N)




#CONFIGURATIONS
#Sequences to be processed by each thread
SEQperTHREAD=2
#The starting sequence number
current_seq_num=1

#This should go to the temp folder of the sRNA generated
output_part="Cleaveland4-partial-result.tsv"
TPLOTS="tplots"
preLabel="tempCleave4_"

base_sRNA=$(basename $sRNA_COMPLETE)
if [[ -n "$degradome" ]]; then
  base_deg=$(basename $degradome)
else
  base_deg=$(basename $degDensity)
fi
base_trans=$(basename $TRANSCRIPT)
OUTPUT_final=${base_sRNA%.*}_${base_deg%.*}_${base_tras%.*}_Cleaveland4_result.tsv

#Number of sequences in fasta file
NUM_SEQ=$(grep -c ">" $sRNA_COMPLETE)


#Running in mode 1
if [[ -z $degDensity ]]; then 
    root_name="${preLabel}${current_seq_num}-"$(( $current_seq_num + $SEQperTHREAD -1))
    #Let's hide all the dirty temp folder under the rug.
    root_dir=".${root_name}"
    #make temp directory for file enter and run  
    cd $DIR
    mkdir -p  ${root_dir}
    extract_x_seqs.sh $sRNA_COMPLETE $current_seq_num $SEQperTHREAD > "${root_dir}/${root_name}.fasta" 
    fasta_part="${root_name}.fasta"
    echo "Spawning new process"
    cd ${root_dir}

    #Run cleaveland in mode 1    
    $cleaveland -e $degradome -u $fasta_part -n $TRANSCRIPT -t -o $TPLOTS > $output_part

    #Advance the current seq num
    current_seq_num=$(( $current_seq_num + $SEQperTHREAD ))	
 
    #Set the density file  
    degDensity=${degradome}_dd.txt
fi


#Continue with degradome density
while [[ "$current_seq_num" -lt "$NUM_SEQ" ]]; do
  root_name="${preLabel}${current_seq_num}-"$(( $current_seq_num + $SEQperTHREAD -1))
  #Let's hide all the dirty temp folder under the rug.
  root_dir=".${root_name}"
  #make temp directory for file enter and run  
  cd $DIR
  mkdir -p  ${root_dir}
  extract_x_seqs.sh $sRNA_COMPLETE $current_seq_num $SEQperTHREAD > "${root_dir}/${root_name}.fasta" 
  fasta_part="${root_name}.fasta"

  #Get script pid and count subprocess beining run
  bash_pid=$$
  children=`ps -eo ppid | grep -w $bash_pid`
  echo $children
  CURRENT_TH_COUNT=$( echo $children | wc -w ) 

  #Check number of running threads and keep checking until 
  #number of threads running is lower then limit
  #Num of thread is checked every 2 sec.
  while [[ $CURRENT_TH_COUNT -ge  $THREADS ]]; do
  	sleep 2
        children=`ps -eo ppid | grep -w $bash_pid`
  	#echo $children
  	CURRENT_TH_COUNT=$( echo $children | wc -w ) 
  done

  #Additional check but it should be ready to run, if at this point
  if [[ $CURRENT_TH_COUNT -le $THREADS ]]; then
   
    echo "Spawning new process"
    cd ${root_dir}
    #Run cleaveland in mode 2
    
    $cleaveland -d $degDensity -u $fasta_part -n $TRANSCRIPT -t -o $TPLOTS > $output_part &
    #Advance the current seq num
    current_seq_num=$(( $current_seq_num + $SEQperTHREAD ))	
  else
    #Shouldn't be necessary but if it happens will wait for all subprocess to end
    wait
  fi
  #For testing
  #sleep 1

done 
wait

#Unification of all bits and pieces
cd $DIR && 
mkdir -p Tplots &&
cp "."${preLabel}*/${TPLOTS}/* Tplots
cat "."${preLabel}*/$output_part > $OUTPUT_final
#get metadate for this
tempResult=tempFileCleave12345678987654321.tmp
head -8 $OUTPUT_final > $tempResult
#Remove metadate headers
sed -i '/^\#/d' $OUTPUT_final
echo "PASS1"
header=$(grep -m1 "^SiteID" $OUTPUT_final )
echo "PASS2"
sed -i '/^SiteID\tQuery/d' $OUTPUT_final
echo "PASS3"
sed -i "1i${header}" $OUTPUT_final
echo "PASS4"
cat $OUTPUT_final >> $tempResult
mv $tempResult $OUTPUT_final
echo "PASS5"
#Once again to ensure i'm in the right directory
cd $DIR &&
#clean up 
rm -r "."${preLabel}*
echo "PASS6"
END=$(date +%s.%N)
DIFF=$( echo "${END} - ${START}" | bc )

echo "Target prediction finished in ${DIFF} secs"


exit 0
