#!/bin/sh
PROG_ID=$1
echo "prog id: ${PROG_ID}"
echo ""
echo "Out:"
echo ""
head log/out.log
echo "" 
echo "Errors:"
echo ""
cat log/error.log
echo ""
echo "Benchmarks"
echo ""
sacct --jobs=${PROG_ID} -o JobId,JobName,AllocCPUs,State,TotalCPU
