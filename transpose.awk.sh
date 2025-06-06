#!/bin/bash
# Usage: ./transpose.awk.sh input.txt > output.txt
awk '
{
    for (i = 1; i <= NF; i++)  {
        a[i, NR] = $i
    }
    max_nf = (NF > max_nf ? NF : max_nf)
}
END {
    for (i = 1; i <= max_nf; i++) {
        for (j = 1; j <= NR; j++) {
            printf "%s%s", a[i, j], (j == NR ? ORS : OFS)
        }
    }
}
' "$1"
