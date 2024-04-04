#! /bin/bash

# makeCombinedIndex.sh creates a list of timesteps from multiple UDAs that can be 
# used in a single index.xml file to refer to multiple UDA directories.
#
# Run makeCombinedIndex.sh in the directory containing the UDAs of interest.

declare -a idx_arr=("CCVars" "SFCXVars" "SFCYVars" "SFCZVars")


if test $# -eq 0 || test "$1" == "-h" -o "$1" == "--h" -o "$1" == "--help" -o "$1" == "-help"; then
   echo ""
   echo "Usage: $0 <UDAs>"
   echo ""
   exit
fi


for dir in $*; do

  if test ! -d $dir; then
    echo "ERROR: '$dir' is not a directory!  Goodbye."
    exit
  fi

done

for v in "${idx_arr[@]}"
do

    idx_filename=`echo "$v.idx"`
    gidx_filename=`echo "$v.gidx"`
    
    echo '<datasets name="'$v'">' >> $gidx_filename

    testing=false
    preTimestep=-1

for dir in $*; do

   for timestep in `grep "timestep href=" $dir/index.xml | awk -F'[=/"]' '{ print $3 }'`; do

      if test -f $dir/$timestep/timestep.xml; then

         tsNum=`echo $timestep | cut -d"t" -f2`
         if test "$tsNum" -le "$preTimestep"; then
            if test "$testing" == "true"; then
               echo "ERROR: $timestep is before previous timestep: t$preTimestep"
            fi
         else

            timeline=`grep "timestep href=" $dir/index.xml | grep $timestep`

            time=`echo $timeline | cut -f3 -d" " | cut -f2 -d"="`
            oldDelt=`echo $timeline | cut -f4 -d" " | cut -f2 -d"=" | cut -f1 -d">"`
            timeNum=`echo $timeline | cut -f4 -d" " | cut -f2 -d">" | cut -f1 -d"<"`

            #echo "time is $time, oldDelt is $oldDelt, timeNum is $timeNum"
            echo "    "\<dataset name='"'$timestep'"' log_time='"'$timeNum'"' url='"file://'$PWD/$dir/$timestep/'l0'/$idx_filename'"/>' >> $gidx_filename
	    
            preTimestep=$tsNum
         fi
      else
         if test "$testing" == "true"; then
            echo "WARNING $dir/$timestep/timestep.xml does not exist"
         fi
      fi
   done
done

echo "</datasets>" >> $gidx_filename

done
