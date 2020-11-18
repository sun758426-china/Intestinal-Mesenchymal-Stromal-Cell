deeptools

multiBamSummary bins -bs 20 --bamfiles ~/CD81_1.bam  ~/CD81_2.bam ~/CD34_1.bam  ~/CD34_2.bam ~/MyoF_1.bam ~/MyoF_2.bam ~/CD138_1.bam ~/CD138_2.bam -p 28 -o ~/results.npz --maxFragmentLength 200 -e 2000 --centerReads --ignoreDuplicates -bl ~/mm10-blacklist.v2.bed

plotCorrelation -in /mnt/d/atac/results.npz -c pearson --skipZeros --removeOutliers -p heatmap --colorMap seismic --plotNumbers -min 0.4 -o /mnt/d/atac/corrplot.eps

#Identify subset specific peaks

bedtools subtract -a ~/CD81_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed -b ~/CD34_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed > ~/CD81_1.bed

bedtools subtract -a ~/CD81_1.bed -b ~/CD138_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed > ~/CD81_2.bed

bedtools subtract -a ~/CD81_2.bed -b ~/MyoF_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed > ~/CD81_2.bed

#Identify subset overlap peaks

bedtools intersect -a ~/CD81_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed -b ~/CD34_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed > ~/CD81_CD34.bed

bedtools intersect -a ~/CD81_CD34.bed -b ~/CD138_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed > /mnt/e/ATAC/CD81_CD34_CD138.bed

bedtools intersect -a ~/CD81_CD34_CD138.bed -b ~/MyoF_Specific/case/rep_concord_merge/mergedReads.PeakCallingFseq.bed > /mnt/e/ATAC/total.bed

#Merge All Peaks for Heatmap

computeMatrix reference-point --referencePoint center -b 2000 -a 2000 -R /mnt/e/ATAC/CD81_2.bed /mnt/e/ATAC/CD34_2.bed /mnt/e/ATAC/CD138_2.bed /mnt/e/ATAC/MyoF_2.bed /mnt/e/ATAC/total.bed -S /mnt/e/ATAC/CD81_1.bw /mnt/e/ATAC/CD81_2.bw /mnt/e/ATAC/CD34_1.bw /mnt/e/ATAC/CD34_2.bw /mnt/e/ATAC/CD138_1.bw /mnt/e/ATAC/CD138_2.bw /mnt/e/ATAC/MyoF_1.bw /mnt/e/ATAC/MyoF_2.bw --skipZeros -o /mnt/e/ATAC/heatmap.gz -p 28

plotHeatmap -m /mnt/e/ATAC/heatmap.gz -out /mnt/e/ATAC/heatmap.pdf --colorList blue,gold,red --sortRegions descend --zMin 0 --zMax 5 --refPointLabel center