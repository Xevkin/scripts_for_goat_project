#!/usr/bin/awk -f

{if ($0 ~ /^@chr/) {k=0} 

if ($0 ~ /chr2/ || $0 ~ /chr1/ || $0 ~ /chr3/ || $0 ~ /chr4/ || $0 ~ /chr5/ || $0 ~ /chr6/ || $0 ~ /chr7/ || $0 ~ /chr8/ || $0 ~ /chr9/) {k=1}

if (k==1) print $0}

