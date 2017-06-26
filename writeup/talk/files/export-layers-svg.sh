#!/bin/bash
# CC-BY https://answers.launchpad.net/inkscape/+question/48164
# $1 is the file to extract layers from
# $2 is not needed but can be used eg with "-d 300" to send the png resolution or any other inkscape switch through
#
# needs xmlstarlet
#
# Notes:
#  - doesn't deal with spaces in object names

type xmlstarlet >/dev/null 2>&1 || { echo >&2 "can't find xmlstarlet"; exit 1; }

TMPFILE=$(mktemp /tmp/output.XXXXXXXXX).svg
cp $1 $TMPFILE

YESLAYERS=${*:2}
ALLLAYERS=$(xmlstarlet sel -t -m "//*[@inkscape:groupmode=\"layer\"]" -v "concat(@inkscape:label,' ')" $TMPFILE|sort -Vr)
NOLAYERS=$(comm -1 -3 <(echo $YESLAYERS|tr ' ' '\n'|sort) <(echo $ALLLAYERS|tr ' ' '\n'|sort))

for layer in $NOLAYERS
do
    echo $layer
    id=$(xmlstarlet sel -t -m "//*[@inkscape:label=\"$layer\"]" -v "@id" $TMPFILE)
    xmlstarlet ed -S -L -d "//*[@id=\"$id\"]" $TMPFILE
done

inkscape --without-gui --export-pdf=/dev/stdout $TMPFILE

# [ -f $TMPFILE ] && rm $TMPFILE
