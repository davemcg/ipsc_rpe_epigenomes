#!/bin/bash

in2csv data/Copy\ of\ Brooks\ retRPE\ Illumina_Sample_Submission_Form_ChIP_CUT\&RUN03.22.19.xls | grep "^Sample" -A 100 | grep -v "^," | grep -v "^\*\*\* Tubes" > data/Copy\ of\ Brooks\ retRPE\ Illumina_Sample_Submission_Form_ChIP_CUT\&RUN03.22.19.csv
