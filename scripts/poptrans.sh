#!/bin/bash
if [ $# = 2 ]
then
    :
else
    printf "Usage: poptrans <mapfile> <ba3 file>\n"
    exit
fi
if [ -e $1 ] && [ -e $2 ]
then
    :
elif [ ! -e $1 ]
then
    printf "$1 not found!\n"
    exit
elif [ ! -e $2 ]
then
    printf "$2 not found!\n"
fi
while IFS= read -r line
do
    if [ ! -z "$line" ]
    then
	regexp=`echo $line | awk '{ print $1 }'`
	popname=`echo $line | awk '{ print $2 }'`
#	re=^9[0-9_]+.*
	awk '/'$regexp'/ { print $1 "\t" p1 "\t" $3 "\t" $4 "\t" $5}' "p1=$popname" $2
    fi
done < "$1"
       
