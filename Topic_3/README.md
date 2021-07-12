---
title: Topic 3
permalink: /Topic_3/
topickey: 3
topictitle: "Preprocessing"
---

### Recorded lecture 2020
<iframe src="https://monash.au.panopto.com/Panopto/Pages/Embed.aspx?id=1ed65bb6-c112-4fb2-afb9-ac8201899f60&autoplay=false&offerviewer=true&showtitle=true&showbrand=false&start=0&interactivity=all" height="405" width="720" style="border: 1px solid #464646;" allowfullscreen allow="autoplay"></iframe>

### Lecture slides
* [Slides](./quality_trimmingMonash.pdf)

### Background reading
* [Data preprocessing](./background_reading/Data_preprocessing.pdf), by Robert Schmieder.
* Del Fabbro C, Scalabrin S, Morgante M, Giorgi FM (2013) [An Extensive Evaluation of Read Trimming Effects on Illumina NGS Data Analysis](./background_reading/journal.pone.0085024.PDF). PLoS ONE 8(12): e85024. doi:10.1371/journal.pone.0085024

# Code break questions 

1) How many sequences so you have in the fasta file ~/Topic_3/data/Pine_reference_rnaseq_reduced.fa?

Hint: wc –l <file name> provides the number of lines in a file

2) How many sequences do you have in the fastq file ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R1.fastq?

Hint: for grep ^ indicates the start of the line and $ indicates the end of the line (e.g. grep ^H*?$ <filename> would find all the lines starting with H and ending in ?)

3) How many sequences contain a base with a Phred score of 2 ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R1.fastq?

To practice unix-based command line, run through some of the examples after the quality trimming section to learn how to navigate through the file system and do basic file manipulation.



# Topic 3: Preprocessing Sequence Data


## Short read data (Illumina paired end)

You will assess and filter two short read datasets today. The first is a Genotype by Sequencing (GBS) dataset from a conifer which used an enzyme called PST1. The second is a whole genome sequence dataset from a bacteria (Staphylococcus aureus).

Before we get started, we should make a directory to contain all of your output files

```bash
mkdir ~/Topic_3/out
```

## GBS dataset

We are going to use a program called Fastqc to assess the quality of our data. It will run via the GUI on your VM. Simply open the terminal and type:

```bash
fastqc
```

A window will open for the Fastqc program. Click "File" and then "Open". Navigate to the ~/Topic_3/data in your home directory and select GBS12_brds_Pi_197A2_100k_R1.fastq

Alternative, you can run Fastqc through command line by typing 

```bash
cd ~/Topic_3/out
fastqc  ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R1.fastq 
```
This will produce an html file that you can then open in your browser. I recommend using the GUI however. 

There are a number of different quality metrics for you to assess using Fastqc. Please watch this [video](https://www.youtube.com/watch?v=bz93ReOv87Y) for a detailed explanation. There are examples of "good" and "bad" data that you can review [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). There is also a written explanation of each metric [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

You can also look at the quality metrics using the program that you will use to trim and filter the data (Prinseq). Many of the metrics are similar but Prinseq supplies a few different plots that you might find helpful.

These commands make a file containing information to produced the Prinseq QC plots in the folder ~/Topic_3/out (note that the "\\" allows use to continue the command on the next line to make it more readable. You can remove the slashes but you must then place the code all on the same line.)

```bash

cd  ~/Topic_3/out
perl ~/Topic_3/scripts/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R1.fastq \
-fastq2 ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R2.fastq \
-graph_data ~/Topic_3/out/GBS12_brds_Pi_197A2_100k_graph.txt 

```

Upload the graph file (GBS12_brds_Pi_197A2_100k_graph.txt) to [here](http://edwards.sdsu.edu/cgi-bin/prinseq/prinseq.cgi) to view/download your graphs (click Get report)

Find a description of the various plot types [here](http://prinseq.sourceforge.net/manual.html#STANDALONE)

Question 1) Examine the QC results files for GBS12_brds_Pi_197A2_100k_R1.fastq and GBS12_brds_Pi_197A2_100k_R1.fastq. Are there any issues with the data? Think carefully about the data and assess if these patterns represent issues with the data that need to be rectified through filtering (hint: GBS reads will start with the same restriction site, which is part of the genomic DNA)? If so, what filters would you apply?


2. Within Prinseq (and most other programs), "trimming" will not affect the number of reads you have, but will alter the characteristics of the sequences, while "filtering" will remove reads that do not pass the criteria. This can result in un-paired reads in the two _1 and _2 files, but some programs (including Prinseq) take care of this by outputting separate files for the good (i.e. paired) reads and the unpaired (i.e. singletons) where one read has been filtered. 

Trim off bases from either end with less than a quality score of 10 (-trim_qual_left and -trim_qual_right), and filter any sequences that have less than 70 base pairs (-min_len 70):

```bash
	
cd  ~/Topic_3/out
perl ~/Topic_3/scripts/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R1.fastq \
-fastq2 ~/Topic_3/data/GBS12_brds_Pi_197A2_100k_R2.fastq \
-log gbslog1 -out_good GBS_filter1 -min_len 70 -trim_qual_left 10 -trim_qual_right 10 

```

Note that you can view the log file (gbslog1) to see the command executed and the output (including the default parameters run)

You should see 2 output files called GBS_filter1_#.fastq and two more that have a "_singletons" suffix. These "singleton" files are those that were left unpaired after the corresponding sequence was filtered from the alternate read direction file.

Experiment with other filters, such as filtering sequences with a mean quality score below some value (here, Q15):

-min_qual_mean 15

Or removing N's from the left and right hand side of the sequences:

-trim_ns_left 3 (Trim poly-N tail with a minimum length of 3 at the 5'-end)
-trim_ns_right 3 (Trim poly-N tail with a minimum length of 3 at the 3'-end)

Low complexity sequences can be filtered using:

-lc_method [dust/entropy] -lc_threshold [7

From the manual: 

	The DUST approach is adapted from the algorithm used to mask low-complexity
	regions during BLAST search preprocessing [5]. The scores are computed based on
	how often different trinucleotides occur and are scaled from 0 to 100. Higher
	scores imply lower complexity. A sequence of homopolymer repeats (e.g. TTTTTTTTT)
	has a score of 100, of dinucleotide repeats (e.g. TATATATATA) has a score around
	49, and of trinucleotide repeats (e.g. TAGTAGTAGTAG) has a score around 32.

	The Entropy approach evaluates the entropy of trinucleotides in a sequence. The 
	entropy values are scaled from 0 to 100 and lower entropy values imply lower 
	complexity. A sequence of homopolymer repeats (e.g. TTTTTTTTT) has an entropy 
	value of 0, of dinucleotide repeats (e.g. TATATATATA) has a value around 16, and
	of trinucleotide repeats (e.g. TAGTAGTAGTAG) has a value around 26.

	Sequences with a DUST score above 7 or an entropy value below 70 can be considered 
	low-complexity. An entropy value of 50 or 60 would be a more conservative choice.


NOTE: prinseq has an order of operations, which it tells you when you run it with the -h option. Be mindful of this, because it will calculate the mean quality score differently if other operations are run first.

There is no best option for trimming/filtering. The choices you make should reflect the application and the state of your data. There is such a thing as too much filtering/trimming. 

Question 2) Try a couple different filtering combinations for the GBS data  (see the [manual](http://prinseq.sourceforge.net/manual.html) for options) and examine the QC plots using Fastqc. Which options you would choose to implement if this was your data?


## WGS data

Now we are going to QC and filter some paired-end whole genome sequence data from a single bacterial sample downloaded from GAGE (Genome Assembly Gold-Standard Evaluations) [website](http://gage.cbcb.umd.edu/data/index.html). The paired end sequence files are ~/Topic_3/frag_1.fastq ~/Topic_3/frag_2.fastq.

First assess the quality using Fastqc.

Question 3) Examine the QC results for ~/Topic_3/data/frag_1.fastq  and ~/Topic_3/data/frag_2.fastq. Are there any issues with the data? How do these data compare to the GBS data above? What filters would you apply if you wanted to use these data for assembly of the bacterial genome?

Use Prinseq to the filter the data. Below is an example command but you should adjust it depending on how you think the data should be filtered:

```bash
cd  ~/Topic_3/out
perl ~/Topic_3/scripts/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ~/Topic_3/data/frag_1.fastq -fastq2 ~/Topic_3/data/frag_2.fastq \
-log wgslog1 -out_good WGS_filter1 -min_len 70 -trim_qual_left 10 -trim_qual_right 10  
	
```

Re-plot the WGS_filter1 data using FastQC. 

Question 4) What improvements have been made? How much sequence data has been removed as a result of your trimming and filtering? 



## Long read data

We going assess the quality of long reads (Nanopore) from Staphylococcus aureus genome. These long reads were taken from [NCBI](https://www.ncbi.nlm.nih.gov/sra/SRR10346110) and we will use them to assemble the bacterial genome in a couple days.

The data can be found here: 

* ~/Topic_3/data/sra_data.fastq

Assess the quality using Fastqc. 

Question 5) What do you notice about these reads compared to the Illumina short read data. Are these differences expected? 

Now we are going to remove potential adapter contamination using [Porechop](https://github.com/rrwick/Porechop), a program designed to specifically removed adapters from Nanopore data. The adapter sequences for Nanopore, were not included in this Fastqc configuration and they might not be apparent in our QC plots.

We are going to output an adapter filtered (sra_fil.fastq) fastq file in the out directory. Note, this takes about 5-10 min.

```bash
cd  ~/Topic_3/out
porechop -i ~/Topic_3/data/sra_data.fastq -o sra_fil.fastq
```

Question 6) Look at the report produced by Porechop, were there many adapter sequences removed? Were they only found at the ends of sequences?

Open the filtered file with Fastqc and observe the plots.

Question 7) Do you think we need to further trim and filter the data before we assemble the genome? If so, use Prinseq to improve the quality of the data and assess the results using Fastqc. 

Here is an example command:

```bash
cd  ~/Topic_3/out
perl ~/Topic_3/scripts/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ~/Topic_3/out/sra_fil.fastq \
-log nanolog -out_good sra_fil2 -min_len 500 -trim_left 50
```


In case you want to set these programs up on your own machine here are the instructions that I used for installation

Prinseq:
	
Download at https://sourceforge.net/projects/prinseq/files/
Unpack prinseq-lite-0.20.4.tar.gz by executing the following command:
```bash
tar -xf prinseq-lite-0.20.4.tar.gz
```

Fastqc:
```bash
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
cd FastQC 
chmod 755 fastqc 
sudo ln -s /home/ubuntu/FastQC/fastqc /usr/bin/
```
	
Porechop:
```bash
git clone https://github.com/rrwick/Porechop.git
cd Porechop
sudo python3 setup.py install
porechop -h
```
