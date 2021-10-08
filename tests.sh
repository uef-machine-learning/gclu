#!/bin/bash


mdfn=data/test_md5.txt
rm $mdfn

rm tmp/s4.part
./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 100 -K 15 --algo m --costf 1 --type distance --seed 1632488925 -g 0.8
echo "S4-cf1" >> $mdfn
md5sum tmp/s4.part >> $mdfn

rm tmp/s4.part
./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 100 -K 15 --algo m --costf 2 --type distance --seed 1632488925 -g 0.8
echo "S4-cf2" >> $mdfn
md5sum tmp/s4.part >> $mdfn

rm tmp/s4.part
./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 100 -K 15 --algo m --costf 4 --type distance --seed 1632488925 -g 0.8
echo "S4-cf4" >> $mdfn
md5sum tmp/s4.part >> $mdfn

diff -w data/test_md5.txt  data/test_md5.txt_s4

