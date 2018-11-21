#!/usr/bin/awk -f

{if ($0 ~ /^@[Cc]hr/) {k=0}

if ($0 ~ /^@AJPT/) {k=0}

if ($0 ~ /^@KB0/) {k=0}

if ($0 ~ /^@[123456789]/) {k=1}

if (k==1) print $0}

