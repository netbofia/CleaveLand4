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
  version="true"
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
  TPLOTS="$2"
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
  pval="$2"
  shift # past argument
  ;;
  -c|--category)
  category="$2"
  shift # past argument
  ;;
  -s|--threads)
  THREADS="$2"
  shift # past argument
  ;;
  -h|--help)
  echo -e "  
  ------------------------------------------
  ${blue}Cleaveland4 threaded verion${NC}
  ------------------------------------------
  Mandatory parameters - 4 of the following:
  -d or -e both are not allowed
  ------------------------------------------

  ${blue}-u|--sRNA          ${NC}| sRNA sequence file (fasta) No spaces or strange chars in seq names
  ${blue}-d|--degDensity    ${NC}| Degradome density file
  ${blue}-e|--degradome     ${NC}| Degradome file (fasta)
  ${blue}-n|--transcriptome ${NC}| Transcriptome file (fasta) Will be stored in memory (RAM) N times
  ${blue}-s|--threads       ${NC}| Number of parallel processes to run
  ---------------------
   Optional args
  ---------------------
  ${green}-h|--help          ${NC}| This message.
  ${green}-v|--version       ${NC}| Display version of Cleaveland (Not working)
  ${green}-q|--quiet         ${NC}| Quiet mode straight to log (Not implemented yet) 
  ${green}-a|--allenSort     ${NC}| Sort by allen score instead of by MFEratio
  ${green}-t|--tabular       ${NC}| Tabular output allways on (Might be removed)
  ${green}-r|--mfer          ${NC}| Minium free energy ratio cutoff. Default: 0.65 [0-1]
  ${green}-o|--tplots        ${NC}| Directory where tplots are saved
  ${green}-g|--gstar         ${NC}| Gstar file (Not implemented)
  ${green}-p|--pval          ${NC}| P-value cutoff. Default: 1 [0-1] 		
  ${green}-c|--category      ${NC}| Category cutoff Default:4 [0-4]
  
  For more information go to: https://github.com/MikeAxtell/CleaveLand4
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
if [[ "$version" == "true" ]]; then
  $cleaveland -v
  exit 0       
fi
if [[ -z "$pval" ]]; then
  pval=1
fi
if [[ -z "$mfer" ]]; then
  mfer=0.65
fi
if [[ -z "$category" ]]; then
  category=4
fi 
if [[ -z "$THREADS" ]]; then
  THREADS=1
fi


#Planed support only for modes 1 and 2 
#Modes 3 and 4 use the existing GSTAr file (Not sure it is a use case for this abroad shouldn't need paralization)
#Support for the -p,-r,-c implemented. -v,-q,-a,-t,-g options may become a reality at least most of them but might take some time.
#Tabular format is actually the norm not sure if worth changing this 
#Some are simple others involve some hacking around with conditions. -g might not be an option


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
if [[ -z "$TPLOTS" ]]; then
  TPLOTS="Tplots"
fi
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
    $cleaveland -e $degradome -u $fasta_part -n $TRANSCRIPT -r $mfer -p $pval -c $category -t -o $TPLOTS > $output_part

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
  	CURRENT_TH_COUNT=$( echo $children | wc -w ) 
  done

  #Additional check but it should be ready to run, if at this point
  if [[ $CURRENT_TH_COUNT -le $THREADS ]]; then
   
    echo "Spawning new process"
    cd ${root_dir}
    #Run cleaveland in mode 2
    
    $cleaveland -d $degDensity -u $fasta_part -n $TRANSCRIPT -r $mfer -p $pval -c $category -t -o $TPLOTS > $output_part &
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
cp "."${preLabel}*/${TPLOTS}/* ${TPLOTS}
cat "."${preLabel}*/$output_part > $OUTPUT_final
#get metadate for this
tempResult=tempFileCleave12345678987654321.tmp
head -8 $OUTPUT_final > $tempResult
#Remove metadate headers
sed -i '/^\#/d' $OUTPUT_final
header=$(grep -m1 "^SiteID" $OUTPUT_final )
sed -i '/^SiteID\tQuery/d' $OUTPUT_final
sed -i "1i${header}" $OUTPUT_final
cat $OUTPUT_final >> $tempResult
mv $tempResult $OUTPUT_final
#Once again to ensure i'm in the right directory
cd $DIR &&
#clean up 
rm -r "."${preLabel}*
END=$(date +%s.%N)
DIFF=$( echo "${END} - ${START}" | bc )

echo "Target prediction finished in ${DIFF} secs"


exit 0
