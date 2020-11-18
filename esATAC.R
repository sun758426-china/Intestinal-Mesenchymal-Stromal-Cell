# Loading Library

library(esATAC)

#Fastq Mapping Genome

CD81vsIEC<- 
  atacRepsPipe2(genome = "mm10",caseFastqInput1=list("~/ATBS02_R1.fastq.gz","~/ATBS07_R1.fastq.gz"),
                caseFastqInput2 =list("~/ATBS02_R2.fastq.gz","~/ATBS07_R2.fastq.gz") ,
ctrlFastqInput1=list( "~/ATBS14_R1_001.fastq.gz","~/ATBS19_R1_001.fastq.gz"),
ctrlFastqInput2 = list("~/ATBS14_R2_001.fastq.gz","~/ATBS19_R2_001.fastq.gz"),
   threads = 28,
   refdir = "~/mm10",
 tmpdir = "~/New_Folder_1")  

# Bam Output

sam2bam(samInput = "~/New_Folder_1/ATBS02_R1.Bowtie2Mapping.sam",bamOutput = "~/ATBS02.bam")

