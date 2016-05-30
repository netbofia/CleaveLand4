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

#Get dir
DIR=$(pwd)
#Start time for benchmarking
START=$(date +%s.%N)

#Fasta with all the sRNA sequences in fasta format
sRNA_COMPLETE=$1
#Check to see which type was given based on the file extension
DEG=$2
#Transcriptome	
TRANSCRIPT=$3
#Check it is higher than the numer of threads on the machine
THREADS=$4

#CONFIGURATIONS
#Sequences to be processed by each thread
SEQperTHREAD=2
#The starting sequence number
current_seq_num=1
#Path to cleaveland (Hardcoded Sorry will be changed to be incorporated to path instead)
cleaveland=~/.Software/CleaveLand4-4.3/CleaveLand4.pl

#This should go to the temp folder of the sRNA generated
output_part="Cleaveland4-result.tsv"
TPLOTS="tplots"
preLabel="tempCleave4_"

base_sRNA=$(basename $sRNA_COMPLETE)
base_deg=$(basename $DEG)
base_trans=$(basename $TRANSCRIPT)
OUTPUT_final=${base_sRNA%.*}_${base_deg%.*}_${base_tras%.*}_Cleaveland4_result.tsv

#Number of sequences in fasta file
NUM_SEQ=$(grep -c ">" $sRNA_COMPLETE)

while [[ "$current_seq_num" -lt "$NUM_SEQ" ]]; do
  echo $current_seq_num
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
   
    echo "Spwaning new process"
    cd ${root_dir}
    #Run cleaveland in mode 2
    
    $cleaveland -d $DEG -u $fasta_part -n $TRANSCRIPT -t -o $TPLOTS > $output_part &
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
metadata=$(head -9  $OUTPUT_final)
#Remove metadate headers
sed -i '/^\#/d' $OUTPUT_final
header=$(grep -m "^SiteID\tQuery")
sed -i '/^SiteID\tQuery/d' $OUTPUT_final
sed -i "1i${header}" $OUTPUT_final
sed -i "1i${metadata}" $OUTPUT_final
#Backup if last sed doesn't do the job.
#temp="Cleaveland4-tempOUTPUT-12345678987654321.tsv"
#echo $header | cat - $OUTPUT_final > $temp
#mv $temp $OUTPUT_final 

#Once again to ensure i'm in the right directory
cd $DIR &&
#clean up 
rm -r "."${preLabel}*

END=$(date +%s.%N)
DIFF=$( echo "${END} - ${START}" | bc )

echo "Target prediction finished in ${DIFF} secs"


exit 0
