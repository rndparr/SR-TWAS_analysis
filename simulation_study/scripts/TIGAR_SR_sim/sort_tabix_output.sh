#!/usr/bin/bash

temp=$1
out=$2
tabix_str=${3:--b2 -e2 -S1}

echo 'Sort/bgzip/tabix-ing output file.'
head -n1 ${temp} > ${out}

tail -n+2 ${temp} | \
sort -nk1 -nk2 >> ${out} && \
rm ${temp} || \
( tail -n+2 ${temp} | \
sort -T ${out%/*} -nk1 -nk2 >> ${out} && \
rm ${temp} )

if [ ! -f "${temp}" ] ; then
    echo 'Sort successful. Bgzip/tabix-ing.'
    bgzip -f ${out} && \
    tabix -f ${tabix_str} ${out}.gz

else
    echo 'Sort failed; Unable to bgzip/tabix output file.'
fi

