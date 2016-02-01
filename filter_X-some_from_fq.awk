#!/usr/bin/awk -f

{if ($0 ~ /^@chr/) {k=0} 

if ($0 ~ /chrX/) {k=1}

if (k==1) print $0}

