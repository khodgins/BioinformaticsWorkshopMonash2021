---
title: "Topic 5: Assembly"
permalink: /Topic_5/
topickey: 5
topictitle: "Assembly"
---

## Recorded lecture 2020

<iframe src="https://monash.au.panopto.com/Panopto/Pages/Embed.aspx?id=835dc4b6-e73c-4892-88cc-ac820189a883&autoplay=false&offerviewer=true&showtitle=true&showbrand=false&start=0&interactivity=all" height="405" width="720" style="border: 1px solid #464646;" allowfullscreen allow="autoplay"></iframe>

## Accompanying material

Slides: [De novo Assembly](./Assembly2021.pdf)

Background reading: 
* Comparison of the two major classes of assembly algorithms: overlap-layout-consensus and de-bruijn-graph [Paper](./background_reading/Briefings in Functional Genomics-2011-Li-bfgp-elr035.pdf). Briefings in Functional Genetics 2011.
* The impact of third generation genomic technologies on plant genome assembly [Paper](./Jiao.pdf). Jiao & Schneeberger. Curr. Op. Plant Biol. 2017.
* [Velvet Manual 1.1](./background_reading/Manual.pdf)
* [Stacks paper](./background_reading/Stacks.pdf)
* [Minimap2 paper](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)
* [Miniasm paper](https://academic.oup.com/bioinformatics/article/32/14/2103/1742895)
* [Racon paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411768/)
* [A comparison of assemblers for long read assembly in prokayotes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6966772/) 

# Topic 5 Assembly

## Code break questions 

1. Write a one liner to find all the overlaps exactly 4 bp in length between CTCTAGGCC and a list of other sequences in the file ~/Topic_5/data/overlaps.fa

2. Find all the unique 9mers in a fasta sequence ~/Topic_5/data/kmer.fa

This might be a tricky one and there are likely many ways to do this. First try out these commands that might help you complete this task. It might also be helpful to determine how many characters are in the sequence (wc -c).
		
Hints: test out the following commands:

```bash
cut -c1- kmer.fa
```

```bash
cut -c1-4 kmer.fa
```
		
```bash
for num in {1..10}
    do
        echo $num >> file.txt
    done
```
What do these commands do? Can you use commands like this to find all the kmers in the sequence?

3. Sort them and keep the unique kmers

Hint: try sort (look up the options)


## Tutorial 

Your goal for today is to assemble a bacterial genome (Staphylococcus aureus), as best you can, using short reads or long reads. You will then evaluate your assemblies to decide on the optimal method.


### Short read assembly

This genome was downloaded from GAGE (Genome Assembly Gold-Standard Evaluations) website (http://gage.cbcb.umd.edu/data/index.html): 


Details for this dataset are as follows: 

Species: Staphylococcus aureus

Actual genome size: 2,860,307 bp

Type: Paired end

Read number (total - including both reads per pair): 1,294,104

Read size (each read): 101 bp 

Insert length (sd): 180 bp (+/-20 bp) 

Question 1) Given the above information, what is the average expected coverage?

The data is located in ~/Topic_5/data/illumina/. One file contains the forward read (frag_1.fastq.gz) and the other file contains the reverse read (frag_2.fastq.gz). Each read has a match from the same fragment in the other file and is in the same order in the matching file.

Step 1. Install Velvet COMPLETED 
(see Manual at http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf):

Velvet is already installed, but this was the command used to install the program:

sudo apt-get install velvet

(note I also placed the zipped program in the ~/Topic_5/scripts folder just in case, but it would have to be installed)


Step 2. The first program you need to run is velveth

Velveth takes in a number of sequence files, produces a hashtable (i.e. all the kmers), then outputs two files in an output directory (creating it if necessary). These files are called Sequences and Roadmaps, and are necessary to run velvetg (the next step). 

The syntax to run velveth is as follows:

velveth "output_directory" "hash_length"

To ensure that each k-mer cannot be its own reverse complement, k (i.e. hash length or kmer length) must be odd.

For all the options simply type velveth and refer to the velvet manual for more details.

Here is an example command line:

Make a directory in your home directory for your output

```bash
mkdir ~/Topic_5/out
```

Move into the ~/Topic_5/out directory and run velveth

```bash
cd ~/Topic_5/out

velveth sa_assembly21 21 -shortPaired -fastq  -separate ~/Topic_5/data/illumina/frag_1.fastq.gz ~/Topic_5/data/illumina/Topic_5/data/frag_2.fastq.gz

```

If this process fails because of memory issues (Segmentation fault) then copy the following directories into your output folder
```bash
cp -r ~/Topic_5/data/sa_assembly21 ~/Topic_5/out/
cp -r ~/Topic_5/data/sa_assembly31 ~/Topic_5/out/
```

21 is the kmer length
-shortParied specifies the types of reads (paired end)
-fastq is the type of sequence format
-separate tells the programs that the paired reads are in two separate files (one with read 1 (frag_1) and one with read 2 (frag_2))

A Roadmap file is produced. This file is used as input for the next stage of velvet. For each k-mer observed in the set of reads, the hash table records the ID of the first read encountered containing that k-mer and the position of its occurrence within that read. This file rewrites each read as a set of original k-mers combined with overlaps with previously hashed reads. For more information see:

http://homolog.us/blogs/blog/2011/12/06/format-of-velvet-roadmap-file/

Step 3. The next step is to make make the graph, simplify the graph, correct errors and resolve repeats. This is done by velvetg.

velvetg is run as follows:

velvetg "velveth_output_directory"

For all the options simply type velvetg and refer to the velvet manual for more details.

Here is an example command line with no options:

```bash
velvetg sa_assembly21
```

Run the above command.

Now contigs.fa appears in your output directory (sa_assembly21). This file contains your assembled genome. A log file, with information on your input parameters as well as some basic metrics of your assembly, and a stats.txt file with information on each node is present also appear. Note that node lengths are given in k-mers units (see below). 

Step 4. Assess the assembly quality

There are several basic metrics to quantify the quality of a genome assembly. N50 is a common statistic similar to a mean or median contig length, but has greater weight given to the longer contigs. It is defined as the contig length at which half the bases in the genome are in contigs that size or larger. Other metrics include the longest contig size, the total size of the assembly and total contig number. 

To assess the assembly quality use the web-based version of [Quast](http://cab.cc.spbu.ru/quast/). Simply drop you genome assembly into the box and click evaluate. 

Question 2) Can you think of other ways to assess assembly quality? What might be the trouble with only focusing on maximizing N50 or N90? 

Question 3) Quantify the assembly metrics for your first assembly that you ran without any options. Pick a couple sets of parameters to run (be sure to try exp_cov auto). Compare the resulting assemblies with one another and discuss which ones seemed to have improved the assembly and why that might be. Be prepared to share your findings with the class. 


Be careful not to use all your hard disk space or RAM. Velvet is relatively fast, but a memory hog. When you are finished the tutorial, you may want to remove files produced by Velvet from your drive to make space. 

### Recap: 

You can check how much space is on your drive by typing: 

```bash
df -h
```

This will provide you information on how large the directory you are in is:

```bash
du -sh 
```

This will tell you how much memory is being used currently by the server:

```bash
free -g 
```

To remove files you can type: 

```bash
rm <filename>
```

To remove all files in your current directory enter the directory and type the following (Be careful not to delete something you need! Every file in that directory will be GONE!)

```bash
rm *
```

To remove an empty directory type:

```bash
rmdir <directory name>
```

To move files or directories:

```bash
mv <old file name> <new file name>
```

### Back to the tutorial

If you want to modify the k-mer you must run velveth again and replace 21 with a new number. Velvet only allows k-mers up to 31 bp in length. 

Question 4) Typically, the longer the k-mer, the better the assembly, until you hit the point of too little coverage. Why might this be?

If you want to modify other parameters you can simply run velvetg without rerunning velveth. 

Some potential parameters to modify include (see the manual for details):
* min_contig_lgth
* cov_cutoff
* ins_length 
* ins_length_sd 
* exp_cov (note that for this parameter you can include an estimated expected k-mer coverage or ask velvet to estimate it from the data by typing -exp_cov auto)

Be sure to try -exp_cov auto. For a more detailed explanation of this parameter see:

http://homolog.us/blogs/blog/2012/06/08/an-explanation-of-velvet-parameter-exp_cov/


## Long reads

We going to use long reads (Nanopore) to assemble a Staphylococcus aureus genome. These long reads were taken from [NCBI](https://www.ncbi.nlm.nih.gov/sra/SRR10346110) 
The estimated genome size of this strain is 2,727,720 bp (circular) and there is a 27,265 bp circular plasmid.

The Minimap2-Miniasm-Racon pipeline is an incredibly fast and memory efficient way of assembling long read data. It is composed of three steps implemented by separate programs.

Step 1. Overlap: Use Minimap2 for fast all-vs-all overlap of raw reads

Step 2. Layout: Miniasm concatenates pieces of read sequences to generate the final assembly graph. 

Step 3. Consensus: Map the raw reads back to the assembly using Minimap2 to quickly generate a high quality consensus using Racon. 

An overview of the programs.

[Minimap2](https://github.com/lh3/minimap2)

Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. Uses include: 
(1) mapping PacBio or Oxford Nanopore genomic reads to a genome**
(2) finding overlaps between long reads (error rate <15%)**
(3) splice-aware alignment of PacBio Iso-Seq or Nanopore cDNA or Direct RNA reads against a reference genome 
(4) aligning Illumina single- or paired-end reads 
(5) assembly-to-assembly alignment
(6) full-genome alignment between two closely related species (divergence <15%)

**We will use Minimap2 for 1 & 2 today.

Minimap2 is much faster than other common long-read mappers (e.g., BLASR, BWA-MEM, NGMLR and GMAP). It even works for short reads (>100bp Illumina). In fact, Minimap2 is three times as fast as BWA-MEM and Bowtie2, and as accurate on simulated data! 

How it works in a nutshell:
1. Identifies k-mers in the reads/reference and picks out a special type called minimizers. It places them in a hash table (key=minimizer k-mer, value=position of the k-mer)
2. Finds all the matches between the query/reference of minimizers (called seeds).
3. If there are collinear stretches of seeds (i.e., k-mer matches in the same order in the query and the reference), an alignment can be conducted to extend from the ends of the seeds in the chain using dynamic programming to link the seeds

[Miniasm](https://github.com/lh3/miniasm)

Miniasm is a very fast OLC-based de novo assembler for noisy long reads. It takes all-vs-all read self-mappings (e.g., from Minimap) as input. It outputs an assembly graph in the GFA format. Unlike most assemblers, Miniasm does not have a consensus step. It simply concatenates pieces of read sequences to generate the final sequences. Consequently the per-base error rate is similar to the raw input reads.

Miniasm outputs unitigs rather than contigs. Unitigs are unique parts in the assembly graph (unitigs end where the overlap data shows multiple, mutually contradictory, paths and so unitigs end at repeats). Contigs, by contrast, may include regions with ambiguous read information, depending on the algorithm used. Consequently unitigs can be shorter than contigs and the Miniasm assemblies more fragmented.

How is works:
1. Crude read selection. For each read, find the longest contiguous region covered by three good mappings. Get an approximate estimate of read coverage.
2. Fine read selection. Use the coverage information to find the good regions again but with more stringent thresholds. Discard contained reads.
3. Generate a string graph. Prune tips, drop weak overlaps and collapse short bubbles (similar to short-read assemblers).
4. Merge unambiguous overlaps to produce unitig sequences.


[Racon](https://github.com/isovic/racon)

Racon is a consensus module for raw de novo DNA assembly of long uncorrected reads. It will take as input raw contigs generated by rapid assembly methods (e.g., Miniasm) which do not include a consensus step. The goal of Racon is to generate a consensus which is of similar or better quality compared to the output generated by assembly methods which employ both error correction of the reads and consensus steps, while providing more rapid results. Both Pacific Biosciences and Oxford Nanopore Technologies can be used.

Racon can be used as a polishing tool after the assembly with either Illumina data or data produced by third generation of sequencing. However, it can also be used as a read error-correction tool. 


Step 1. Overlap: First use Minimap2 to map the Nanopore reads against themselves. This will identify overlaps between the reads. 

First move into the directory with the data

```bash
cd ~/Topic_5/data/nanopore/
```

Run Minimap2 
```bash
minimap2 -x ava-ont sra_data.fastq.gz sra_data.fastq.gz | gzip -1 > ~/Topic_5/out/sa_minimap.paf.gz
```

The above command will compare all reads in the fastq file (sra_data.fastq) against themselves. The output will be passed to gzip (gzip -1) to compress the output and redirect (>) the output into a file called sa_minimap.paf.gz.

Tip: Always check that you produced the expected output and that the contents of that file is as expected. For example, to check you have a file called sa_minimap.paf.gz with data in it in the ~/Topic_5/out/ type:

ls -l ~/Topic_5/out/

Step 2. Layout: Use the Minimap2 output and the trimmed reads to assemble unitigs with Miniasm.

Miniasm is a very fast OLC-based de novo assembler for noisy long reads. It takes all-vs-all read self-mappings as input and outputs an assembly graph in the GFA format. 

```bash
miniasm -f sra_data.fastq.gz ~/Topic_5/out/sa_minimap.paf.gz > ~/Topic_5/out/sa_miniasm.gfa
```

Contrary to most assemblers, Miniasm does not have a consensus step. It simply concatenates pieces of read sequences to generate the final unitig sequences. Therefore we need to use another tool to provide the final consensus. 

Step 3. Consensus: Align the reads back to the layout produced by Miniasm using Minimap2 and use Racon to produce the final assembly

Step3a. GFA to fasta (produce a fasta file of the unitigs from the graph file)
```bash
awk '/^S/{print ">"$2"\n"$3}' ~/Topic_5/out/sa_miniasm.gfa > ~/Topic_5/out/sa_unitigs.fasta
```

Step 3b. Map the reads (sra_data.fastq) to the unitigs (sa_unitigs.fasta) using Minimap2
```bash
minimap2 ~/Topic_5/out/sa_unitigs.fasta sra_data.fastq.gz  > ~/Topic_5/out/sa_reads.gfa1.paf
```

Step 3c. Error correct the unitigs (sa_unitigs.fasta) using the mapping (sa_reads.gfa1.paf)
```bash
racon -t 8 sra_data.fastq.gz ~/Topic_5/out/sa_reads.gfa1.paf ~/Topic_5/out/sa_unitigs.fasta > ~/Topic_5/out/sa_ont_assembly.fasta
```

Step 4. Assess the quality of your long read assembly from error prone reads using [Quast](http://cab.cc.spbu.ru/quast/) 

Step 5. Use a program called [Bandage](https://github.com/rrwick/Bandage) to view the assembly graph files from your various assemblies using a GUI. The sequences are represented by coloured bars (nodes), and alternate paths connecting the sequences (edges) are represented by black lines.

Open the program by typing "Bandage" in the terminal at the prompt. A window will pop up. Click "File", "Open" and select your graph file of interest. Click on "draw graph"

Graph files produced in the assembly process above:

~/Topic_5/out/sa_assembly21/LastGraph (Velvet de bruijn graph) (compare with the LastGraph in sa_assembly31 folder if you conducted that assembly)
~/Topic_5/out/sa_miniasm.gfa (Miniasm)

For a more detailed description of what the assembly graphs are showing you, see the Bandage [wiki](https://github.com/rrwick/Bandage/wiki/Simple-example)

Question 5) How does your long read assembly (and assembly graph) compare to the best genome assembly derived from short reads? Why do you think there are these differences?


## Advanced:

Test out [Minipolish](https://github.com/rrwick/Minipolish#method) which performs iterative Racon polishing rounds 

usage: minipolish [-t THREADS] [--rounds ROUNDS] [--pacbio] [-h] [--version] reads assembly

reads                          Long reads for polishing (FASTA or FASTQ format)
assembly                       Miniasm assembly to be polished (GFA format)

Test if Minipolish improves the assembly:
```bash
cd ~/Topic_5/data/nanopore/
minipolish -t 8 sra_data.fastq.gz ~/Topic_5/out/sa_miniasm.gfa > ~/Topic_5/out/sa_polished.gfa
awk '/^S/{print ">"$2"\n"$3}' ~/Topic_5/out/sa_polished.gfa > ~/Topic_5/out/sa_polished.fasta
```
Upload the assembly (sa_polished.fasta) to Quast to examine the metics and examine the polished graph in Bandage (sa_polished.gfa)

Finally, you can test if the filtering of the Nanopore data (Topic_3) will improve the assembly. Instead of running each step separately, you can run the entire Minimap2+Miniasm+Minipolish pipeline using a single wrapper (miniasm_and_minipolish.sh) found in ~/Topic_5/scripts.

```bash
cd ~/Topic_5/scripts
sh miniasm_and_minipolish.sh ~/Topic_3/out/sra_fil2.fastq 8 > ~/Topic_5/out/sa_filt_polished.gfa
awk '/^S/{print ">"$2"\n"$3}' ~/Topic_5/out/sa_filt_polished.gfa > ~/Topic_5/out/sa_filt_polished.fasta
```
Upload the assembly (sa_filt_polished.gfa) to Quast to examine the metics and examine the polished graph in Bandage (sa_filt_polished.gfa).


####

How to download and install Minimap2, Miniasm, Racon & Bandage
- This was already done for you but incase you need to do it in the future, here is what was done:

```bash
git clone https://github.com/lh3/minimap2 && (cd minimap2 && sudo make)
git clone https://github.com/lh3/miniasm  && (cd miniasm  && sudo make)
```
(mv into /user/local/bin/)

```bash
git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Bandage installation
* Ensure the package lists are up-to-date: sudo apt-get update
* Install prerequisite packages: sudo apt-get install build-essential git qtbase5-dev libqt5svg5-dev
* Download the Bandage code from GitHub: git clone https://github.com/rrwick/Bandage.git
* Open a terminal in the Bandage directory.
* Set the environment variable to specify that you will be using Qt 5, not Qt 4: export QT_SELECT=5
* Run qmake to generate a Makefile: qmake
* Build the program: make
* Bandage should now be an executable file.
* Optionally, copy the program into /usr/local/bin: sudo make install. The Bandage build directory can then be deleted.
