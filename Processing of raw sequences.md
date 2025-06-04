## Setup 
- create directory to work in, name this folder after the data
```bash
mkdir XXX
cd XXX
```
- Check that the following packages are installed: 
	- [sra-tools](https://github.com/ncbi/sra-tools)
	- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
	- [vsearch](https://github.com/torognes/vsearch)
	- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
	- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
	- [LULU](https://github.com/tobiasgf/lulu) (R)
	- [CREST4](https://github.com/xapple/crest4)
- check you have a version of the silvamod database downloaded. (I used silvamod138pr2.fasta)
- move fastq files, **silvamod138pr2.fasta**, **filter_seq_by_OTUtable.py**, and **LULU.R** to the working directory

#### For publicly available data:
download data from NCBI/ENA archives
`fasterq-dump SRRXXXXXXXXXX SRRXXXXXXXX ....

unzip fastq.gz files
`gunzip *.gz	
#### For raw sequences pulled from NELS (i.e. Jørgensen and Zhao, 2016) CGB5 and blank from CGB6(tag 70)
convert from bam to fastq. Need [samtools](https://www.htslib.org/) package
```bash
ls *bam | parallel -j 2 "samtools fastq {} > {.}.fastq"
```
## FastQC
- Using FastQC (Babraham Bioinformatics: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), check quality of reads in fastq format and at what length quality deteriorates to determine the length to trim. (if using downloaded data check for adapter content and add to cutadapt if additional adapters need removed). **Note 1**: Our local data was primarily run on Ion Torrent, so will have worse quality when viewed on FastQC as this is tailored for Illumina data.
`fastqc
---
# Pipeline
This pipeline is based on the alternative VSEARCH pipeline. VSEARCH was developed as a free and open source alternative USEARCH, but also claims other benefits like speed and accuracy. 
https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5075697/

________________________________________________________________________
# START OF PIPELINE
This pipeline is a bash script that I ran directly in the terminal, but can also be written into and run as a shell script.
``` bash
conda activate pipeline16s

! /bin/bash

VSEARCH=$(which vsearch)
THREADS=6

---
# Amount reads in fastq file.  
for f in *.fastq; do  
 awk '{s++}END{print s/4}' $f  
done
```
##### Record the number of reads at the beginning:
It is important to follow how the data is being trimmed, so the first thing I do is start a catalog of how many reads are lost in an excel file to track any weird samples. There is likely an efficient line of code to do this, but I just copy and paste the output of the previous chunk of code into excel file for each dataset. It will be in the order of the fastq files being read so you can do "`ls -1`"  to get the list for the rownames and then follow with the values
## Remove primers 
**Note**: In this study, we used the forward primers. Since our data was from Ion Torrent PGM, we only had forward reads and only used forward reads from all studies for consistency).
#### List of forward primers

| Study                  | 16S region | Forward Primer                     | Reverse Primer                      |
| ---------------------- | ---------- | ---------------------------------- | ----------------------------------- |
| **This Paper**         | V4         | 519F (5'-CAGCMGCCGCGGTAA-3´)       | 805R (5'- GACTACHVGGGTATCTAATCC-3´) |
| **Jørgensen and Zhao** | V4         | 519F (5'-CAGCMGCCGCGGTAA-3´)       | 805R (5'- GACTACHVGGGTATCTAATCC-3´) |
| **Zhao**               | V4         | 519F (5'-CAGCMGCCGCGGTAA-3´)       | 805R (5'- GACTACHVGGGTATCTAATCC-3´) |
| **Møller**             | V4         | 519F (5'-CAGCMGCCGCGGTAA-3´)       | 805R (5'- GACTACHVGGGTATCTAATCC-3´) |
| **Bergsten**           | V4         | 515F (5´-GTGCCAGCMGCCGCGGTAA-3´)   | 806R (5´-GGACTACHVGGGTWTCTAAT-3´)   |
| **Lee**                | V4         | 515F (5´-GTGCCAGCMGCCGCGGTAA-3´)   | 806R (5´-GGACTACHVGGGTWTCTAAT-3´)   |
| **Motamedi**           | V4         | 515F (5´-GTGCCAGCMGCCGCGGTAA-3´)   | 806R (5´-GGACTACHVGGGTWTCTAAT-3´)   |
| **Zhang**              | V4         | 515F (5´-GTGYCAGCMGCCGCGGTAA-3´)   | 806R (5´-GGACTACNVGGGTWTCTAAT-3´)   |
| **Vuillemin**          | V4         | 515F (5´-GTGYCAGCMGCCGCGGTAA-3´)   | 806R (5´-GGACTACNVGGGTWTCTAAT-3´)   |
| **Orcutt**             | V4-V5      | 515F (5´-GTGYCAGCMGCCGCGGTAA-3´)   | 926R (5´-CCGYCAATTYMTTTRAGTTT-3´)   |
| **Wee**                | V4-V5      | 515F (5´-GTGYCAGCMGCCGCGGTAA-3´)   | 926R (5´-CCGYCAATTYMTTTRAGTTT-3´)   |
| **San-Saez****         | V4-V5      | 515F (5´-GTGYCAGCMGCCGCGGTAA-3´)   | 926 R (5´-CCGYCAATTYMTTTRAGTTT-3´)  |
| **Schauberger**        | V4-V5      | 515F (5´-GTGCCAGCMGCCGCGGTAA-3´)   | 926R (5´-CCGYCAATTYMTTTRAGTTT-3´)   |
| **Hiraoka**            | V4-V5      | 530F (5´-́BTGBCAGHMGHHDCGG-3´)     | 907R: (5´-CCGYCWATTYMYTHRARTTT-3´)  |
| **Ramirez**            | V1-V3      | 28F (5´- GAGTTTGATCNTGGCTCA G -3´) | 388R (5´- TGCTGCCTCCCGTAGGAGT -3´)  |
| **Bellec**             | V3-V4      | 341F (5´-CCTACGGGNGGCWGCAG-3´)     | 785R (5´-GACTACHVGGGTATCTAATCC-3´)  |
| **Wegener****          |            | 341F (5´-CCTACGGGNGGCWGCAG-3´)     | 785R (5´-GACTACHVGGGTATCTAATCC-3´)  |
** These two studies already had the primers removed from the sequences in the public repository, thus this first step was skipped and the fastq files were renames to add "2-" to the start to continue with the pipeline. 
### Change -g to appropriate forward primer sequence
``` bash
for f in *.fastq; do
    s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'.' '{print "2-" $1}' <<< $s)
     echo $s

     # Primers used.
     # Forward primer: 519f (5-́###-3)
     # Referse primer: 805r (5-́###-3)


 cutadapt -j 3 -g XXXXXXXXXXX --no-indels --discard-untrimmed -o $s.noprimers.fastq $f
 
 done
```

Motamedi had dual-indexed primers so the script was slightly different for this dataset, this line of code was pulled from the supplementary data from the paper itself: 
```bash
cutadapt --discard-trimmed --error-rate 0.10 -g ^TTAGAWACCCVHGTAGTCCGGCTGACTGACT -o $s.noprimers.fastq $f done
```
##### Print number of reads remaining:
*It is important to follow how the data is being trimmed, so I keep a catalog of how many reads are lost in an excel file to track any weird samples. There is likely an efficient line of code to do this, but I copy and paste the output of the following code into excel file for each dataset.* 
```bash
# Amount reads in fastq file.  
for f in 2-*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done
```
---
## Trim the sequences at 220bp.  
 Based on what is seen in FastQC,  "--fastq_trunclen" can be shortened or lengthened to where the reads start to deteriorate. 220 is a good length usually and what we used for the data in this study. 
``` bash 
 for f in 2-*.fastq; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "3-" $2}' <<< $s)

     $VSEARCH --threads $THREADS \
         --fastq_filter $f \
         --fastq_maxns 0 \
         --fastq_trunclen 220 \
         --fastqout $s.trimmed.fastq
 done
```
##### Print number of reads remaining:
*Add the output to the excel file to follow how the data is being trimmed* 
``` bash
for f in 3-*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done
```
---
## Quality filtering at maxee = 2.
 This sets the maximum errors allowed. It correlates the the Qulaity/Phred score given to each base - which is the probability the base is incorrect.  https://www.drive5.com/usearch/manual/exp_errs.html default is 2, but this number can be lowered for more aggressive filtering 
``` bash
  for f in 3-*.fastq; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "4-" $2}' <<< $s)
     l=$(awk -F'-' '{print $2}' <<< $s)
     l=$(awk -F'.' '{print $1}' <<< $l)
     $VSEARCH --threads $THREADS \
         --fastq_filter $f \
         --relabel $l. \
         --fastq_maxee 2 \
         --fastaout $s.fasta \
         --fasta_width 0
 done
```

##### Print number of reads remaining:
*Add the output to the excel file to follow how the data is being trimmed* 
```bash
for f in 4-*.fasta; do  
 awk '{s++}END{print s/2}' $f  
done
```
---
## Dereplicate at sample level and relabel with sample_n
This removes identical sequences within a sample
``` bash
 for f in 4-*.fasta; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "5-" $2}' <<< $s)
     l=$(awk -F'-' '{print $2}' <<< $s)
     l=$(awk -F'.' '{print $1}' <<< $l)
     $VSEARCH --threads $THREADS \
         --derep_fulllength $f \
         --strand plus \
         --sizeout \
         --relabel $l. \
         --fasta_width 0 \
         --output $s.derep.dfa
 done 
```
##### Print number of reads remaining:
*Add the output to the excel file to follow how the data is being trimmed* 
```bash
for f in 5-*.dfa; do  
 awk '{s++}END{print s/2}' $f  
done
```
---
## Cat all fasta files
combine all files into one fasta file. At this point there is no need to print out results, but it is important to look at the reads removed and keep an eye for anything weird, record it in the excel file if you see something.
``` bash
cat 5-*.dfa > all.fasta
```
---
## Remove unneeded files.
``` bash
rm 2-*.fastq 3-*.fastq 4-*.fasta 5-*.dfa
```
---
## Dereplication of the fasta file.
Removes identical sequences. **derep.uc** is a way of mapping which sequences goes to which sample and how many times
 ``` bash 
 $VSEARCH --derep_fulllength all.fasta \
     --threads $THREADS \
     --minuniquesize 2 \
     --sizein \
     --sizeout \
     --fasta_width 0 \
    --uc derep.uc \
    --output derep.fasta
```
---
## Clustering at 97% similarity.
De Novo clustering of sequences with 97% similarity as default (similarity level can be adjusted).  See Rognes et. al (2016) for more in depth explanation
 ```bash
 $VSEARCH --cluster_size derep.fasta \
     --threads $THREADS \
     --id 0.97 \
     --strand plus \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --centroids centroids.fasta
```
---
## Sort centroids and remove singletons.
This sorts the sequences by abundance and removes the sequences found only once
 ```bash
 $VSEARCH --sortbysize centroids.fasta \
     --threads $THREADS \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --minsize 2 \
     --output sorted.fasta
```
---
## denovo chimera detection.
Looks for and removes chimeras formed from clustering
``` bash
$VSEARCH --uchime_denovo sorted.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --nonchimeras denovo.nonchimeras.fasta
```
---
## Reference chimera detection against SILVA138.
Blasts sequences against SILVA138 database. Removes chimeras
 ``` bash
 $VSEARCH --uchime_ref denovo.nonchimeras.fasta \
     --threads $THREADS \
     --db silvamod138pr2.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --dbmask none \
     --nonchimeras nonchimeras.fasta
```
---
## Relabel OTUs.
*Change "XXX" to unique file name for dataset*
 ``` bash
 $VSEARCH --fastx_filter nonchimeras.fasta \
     --threads $THREADS \
     --sizein \
     --fasta_width 0 \
     --relabel OTU_ \
     --fastaout XXX_otus.fasta
```
 ---
## map sequences to OTU
*Change "XXX" to unique file name for dataset*
 ``` bash
 $VSEARCH --usearch_global all.fasta \
     --threads $THREADS \
     --db XXX_otus.fasta \
     --id 0.97 \
     --strand plus \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --dbmask none \
     --otutabout XXX.otutab.txt 
```
#### Sort OTU table numerically 
 ``` bash
 sort -k1.5n XXX.otutab.txt > XXX.otutab.sorted.txt
 ```
---
# LULU cleanup
make blast db in terminal using last .fasta file from "relabel OTUs" step
``` bash
makeblastdb -in XXX_otus.fasta -dbtype nucl

blastn -db XXX_otus.fasta -outfmt '6 qseqid sseqid pident' -out XXX_LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query XXX_otus.fasta -num_threads 3
```

### Open R Studio
In R - Run LULU.R script and change "XXX" to unique name used in output from mapping sequences to OTUs 

#### LULU.R script
``` R
setwd("/home/hba052/Documents/University/Sequencing/data/Deep_biosphere/XXX")

require(lulu)
require(methods)

matchlist = read.delim("XXX_LULU_match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otus.all = read.delim("XXX_all.otutab.sorted.txt",row.names=1,header=T,sep="\t")
curated_result <- lulu(otus.all,matchlist, minimum_match = 97)
lulus = curated_result$curated_table
write.table(data.frame("OTU"=rownames(lulus),lulus),"XXX_table_curated.tsv", row.names=FALSE, 
            quote=F, sep="\t")

  curated_result$curated_table 
```
### Remove sequences from fasta file that got removed using LULU
(in Python)
``` bash
python filter_seq_by_OTUtable.py XXX_otus.fasta XXX_table_curated.tsv > XXX_OTUs_curated.fasta
```
________________________________________________________________________
# CREST Classifying
Using CREST4 which is installed locally on my computer in a conda environment. See CREST4 documentation on how to install and the requirements. 

#### Activate environment
``` bash
conda activate condaHRB
```
#### Run CREST4
``` bash
crest4 -f XXX_OTUs_curated.fasta
```
#### Rename the OTU files
*Change "XXX" to unique file name for dataset*
"**XXX_OTU_table.csv "
_______________________________________________________________________
# Create a merged file for assignments and OTUs

Merge the OTU table with the assignments from CREST classification into single spreadsheet. Save the table as a .csv and carefully look through for any errors. Double check sample names and that the data makes sense. 

_______________________________________________________________________
# References

