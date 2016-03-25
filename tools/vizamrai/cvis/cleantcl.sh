#/bin/sh

for i in *.tcl 
do
echo "Cleaning $i"
sed '/cvisTrace/d' $i > /tmp/temp.$$
sed '/tracervar/d' /tmp/temp.$$ > $i
done
