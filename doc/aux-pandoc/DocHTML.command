#!/bin/bash

make_doc()
{
    for file in "$@"
    do
	PIPPO=`basename $file .txt`
	echo "$PIPPO"
	txt2tags -t html -i $PIPPO.txt -o ../html/$PIPPO.html
#	mv $PIPPO.conv.txt $PIPPO.txt
#	HTMLDocConvert $PIPPO
    done
}

cd txt; make_doc *.txt
