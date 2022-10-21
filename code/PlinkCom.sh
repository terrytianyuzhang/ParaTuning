#!/bin/bash
# for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
# do
#     plink2 --pfile /raid6/Ron/prs/data/bert_sample/CEU.TUNE/CEU.TUNE --chr $chr --out /raid6/Tianyu/PRS/bert_sample/CEU.TUNE/CHR/CEU.TUNE-chr${chr} --make-bed
# done

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    plink2 --pfile /raid6/Ron/prs/data/bert_sample/YRI.TRN/YRI.TRN --chr $chr --out /raid6/Tianyu/PRS/bert_sample/YRI.TRN/CHR/YRI.TRN-chr${chr} --make-bed
done

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    plink2 --pfile /raid6/Ron/prs/data/bert_sample/CEU.TRN/CEU.TRN --chr $chr --out /raid6/Tianyu/PRS/bert_sample/CEU.TRN/CHR/CEU.TRN-chr${chr} --make-bed
done



# plink2 --pfile /raid6/Ron/prs/data/bert_sample/CEU.TUNE/CEU.TUNE --chr 22 --out /raid6/Tianyu/PRS/bert_sample/CEU.TUNE/CHR/CEU.TUNE-chr22 --make-bed
