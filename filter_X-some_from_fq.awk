#!/usr/bin/awk -f

{if ($0 ~ /^@[Cc]hr/) {k=0} 

if ($0 ~ /[Cc]hrX/) {k=1}

if (k==1) print $0}

