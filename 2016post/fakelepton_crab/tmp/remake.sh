#!/bin/bash
for file in `cat job`
do
echo $file
crab remake --task=$file
dir=`echo $file | awk -F yian_ '{print $2}'`
echo $dir
crab kill -d $dir
rm -r $dir
done

