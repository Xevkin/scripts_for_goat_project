#!/usr/bin/awk -f

{if ($0 ~ /^@[Cc]hr/) {k=0} 

if ($0 ~ /^@AJPT/) {k=0}

if ($0 ~ /^@KB0/) {k=0}

if ($0 ~ /[Cc]hr2/ || $0 ~ /[Cc]hr1/ || $0 ~ /[Cc]hr3/ || $0 ~ /[Cc]hr4/ || $0 ~ /[Cc]hr5/ || $0 ~ /[Cc]hr6/ || $0 ~ /[Cc]hr7/ || $0 ~ /[Cc]hr8/ || $0 ~ /[Cc]hr9/) {k=1}

if (k==1) print $0}

