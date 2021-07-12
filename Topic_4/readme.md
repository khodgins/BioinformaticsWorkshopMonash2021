---
title: "Topic 4: Sequence Alignment"
permalink: /Topic_4/
topickey: 4
topictitle: "Sequence Alignment"
---

### Accompanying material

* [Slides](./Topic_4.pdf)


Today we're going to align sequence data to the sunflower reference genome useing BWA and explore what a BAM file is.

The first step is to set up a directory structure so the resulting files will be organized and copy the raw data to your home directory.

```bash

#Navigate to the Topic_4 folder in your home working directory
cd ~/Topic_4

#make a red directory for your sunflower genome reference
mkdir ref

#Copy the reference directory to the ref directory
cp /mnt/workshop/data/ref/HanXRQr1.0-20151230.1mb.fa  ref/

#Copy the fastq files to this directory
cp -r /mnt/workshop/data/fastq ./

#Make a folder to put the bam
mkdir bam

```
Double check and see if the files have copied over and the files are the correct size as those in the /mnt/workshop/data/fastq folder.

Once that is done, we have to index our reference genome.

```bash
#Index the reference for BWA. 

bwa index ref/HanXRQr1.0-20151230.1mb.fa

```
Now finally we can run BWA and align our raw sequence data 
```bash
bwa mem \
  ref/HanXRQr1.0-20151230.1mb.fa \
  fastq/ANN1133.R1.fastq \
  fastq/ANN1133.R2.fastq \
  -t 2 \
  -R '@RG\tID:ANN1133\tSM:ANN1133\tPL:illumina\tPU:SBS\tLB:ANN1133_lib' \
  > bam/ANN1133.sam
  
```
Lets break this command down since it has several parts:
**bwa** <= We're calling the program _bwa_ from the directory _/usr/local/bin/_. You can call this program no matter where you are in the file system. When programs are installed in the /usr/loca/bin/ directory the path does not need to be specified.

**mem** <= This is the bwa command we are calling. It is specific to bwa and not a unix command.

**\\** <= Having this at the end of the line tells the shell that the line isn't finished and keeps going. You don't need to use this when typing commands in, but it helps break up really long commands and keeps your code more organized.

**ref/HanXRQr1.0-20151230.1mb.fa** <= This is the reference genome. We're using a relative path here so you need be in ~/Topic_4 or it won't be able to find this file.
  
**fastq/ANN1133.R1.fastq.gz** <= This is the forward read (e.g. read 1)  set for the first sunflower sequence sample. It's also a relative path and we can see that the file has been gzipped (since it has a .gz ending).

**fastq/ANN1133.R2.fastq.gz** <= This is the reverse read (e.g. read 2)  set for the first sunflower sample.
  
**-t 2** <= This is telling the program how many threads (i.e. cpus) to use. In this case we're only using two because we have limited resources on this VM.

**-R '@RG\tID:ANN1133\tSM:ANN1133\tPL:illumina\tPU:SBS\tLB:ANN1133_lib'** <= This is adding read group information to the resulting sam file. Read group information lets other programs know what the sample name along with other factors. It is necessary for GATK to run later on.

**> bam/ANN1133.sam** <= This is directing the output of the program into the file _bam/ANN1133.sam_

We now have our reads aligned to the genome in a human readable format (sam) instead of binary format (bam) which we will use later. Generally we keep our data in bam format because its more compressed but we can use this opportunity to better understand the format. 

```bash

#Lets examine that bam file
less -S bam/ANN1133.sam
#Notice the @PG line that includes the program call that created the sam file. 
#This is useful for record keeping.

```
Lets examine the sam file. It contains all the information on the reads from the fastq file, but also alignment information. 
### Questions:
1. How are reads ordered in the sam file (hint: compare the sam to the original fastq)? 
2. What does the 6th column represent? What would the string "1S93M6S" mean?
3. What are some possible reasons why mapping quality could be low for a particular read?
4. What percent of your reads mapped to the genome? Hint: <span>samtools flagstat</span>{: .spoiler}

```bash

#The next step is to convert our same file (human readable) to a 
#bam file (machine readable) and sort reads by their aligned position.
samtools view -bh bam/ANN1133.sam | samtools sort > bam/ANN1133.sort.bam 
```
With this command we're using the pipe "|" to pass data directly between commands without saving the intermediates. This makes the command faster since its not saving the intermediate file to hard disk (which is slower). It can be more risky though because if any steps fails you have to start from the beginning. 


Next we want to take a look at our aligned reads. First we index the file, then we use samtools tview.
```bash
samtools index bam/ANN1133.sort.bam  
samtools tview bam/ANN1133.sort.bam  --reference ref/HanXRQr1.0-20151230.1mb.fa
#use ? to open the help menu. Scroll left and right with H and L. 
#Try to find positions where the sample doesn't have the reference allele. 
```





We have 10 samples, so we don't want to have to type these commands out 10 times. Write a bash script to produced a sorted bam file for each sample.

Before you continue be sure that you have reviewed the [loops](http://swcarpentry.github.io/shell-novice/05-loop/index.html) and [shell scripts](http://swcarpentry.github.io/shell-novice/06-script/index.html) tutorial from Software Carpentry that we went over on the first day. Remember that they use files found in ~/shell-lesson-data.

After you complete the loop lesson try this extension. Create a list of the files in the creature directory and save it to a file.
```bash
cd ~/shell-lesson-data/creatures
ls * > file.txt
cat file.txt
```
Now create a for loop that will take the file list as input, print the file name and then the first 10 lines of the file to the screen.
```bash
for file in `cat file.txt`
do
  echo $file
  head $file
done

#remove the file with the list of files
rm file.txt
```
The `cat file.txt` takes each line of the file file.txt. Each line then becomes the variable "$file" for each iteration of the loop. The file name will be printed (echo $file) and then the first 10 lines will be printed to the screen (head $file).

A bash script is a plain text file (i.e. not rich text, nor a word doc) which contains bash commands. You can create the file on your computer and copy it over, or you can edit it directly from the server with one of the installed editors (this is covered in [topic 2, Editing](../Topic_2/#editing). The name of the file is up to you, but bash scripts are given the `.sh` extension by convention:

HINTS:
  * Write out the list of commands that you would need to align one file. What would you need to change in the command to align the second file or the third (what parts of the file name changes versus stays the same from sample to sample)?
  * You will need to use a loop 
  * Consider using variables for directory paths e.g., "ref=~/Topic_4/ref/"
  {: .spoiler}

MORE HINTS:
  * Can you make a file that contains a list of all the prefixes to the fastq files (this is the part of the name that will change from file to file e.g., ANN1133). You will need to use "ls" "grep" and "sed" (find and replace function) and redirect (>>) the output to a file (e.g., named filelist.txt)
  * Now try and take that file with the list of prefixes (e.g., filelist.txt) and print the first lines of each fastq file to the screen
{: .spoiler}
   ```bash
        for fname in `cat filelist.txt`
	do 
	  head -n1 ~/Topic_4/fastq/$fname.R1.fastq 
	done 
 ```	
 {: .spoiler}
  
  * You can do pathname manipulation with `basename` and `dirname` (see manual pages):

    ```
    dirname a/b/c          # prints a/b
    basename a/b/c.gz      # prints c.gz
    basename a/b/c.gz .gz  # prints c

    fpath=/path/to/the/file.gz
    base=$(basename "$fpath")    # assign "file.gz" to variable base
    echo "$base"                 # prints file.gz
    ```
  {: .spoiler}

<details>
<summary markdown="span">**Answer**
</summary>
```bash
  #First set up variable names
  bam=~/Topic_4/bam
  fastq=~/Topic_4/fastq
  ref_file=~/Topic_4/ref/HanXRQr1.0-20151230.1mb.fa

  #Then get a list of sample names, without suffixes and place it in the bam folder
  ls $fastq | grep R1.fastq | sed s/.R1.fastq//g > $bam/samplelist.txt

  #Then loop through the samples
  for name in `cat $bam/samplelist.txt`
  do
    bwa mem \
    -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" \
    $ref_file \
    $fastq/$name.R1.fastq \
    $fastq/$name.R2.fastq \
    -t 1 > $bam/$name.sam;

    samtools view -bh $bam/$name.sam |\
    samtools sort > $bam/$name.sort.bam;
    samtools index $bam/$name.sort.bam
  done 

```
</details>

After your final bam files are created, and you've checked that they look good, you should remove intermediate files to save space. You can build file removal into your bash scripts, but it is often helpful to only add that in once the script works. It's hard to troubleshoot a failed script if it deletes everything as it goes.
### By tomorrow, you should have created cleaned bam files for all samples so we can move on to variant calling.

## Discussion Questions
1. Is an alignment with a higher percent of mapped reads always better than one with a lower percent? Why or why not?
2. I want to reduce the percent of incorrectly mapped reads when using BWA. What setting or settings should I change in BWA?
3. What are two ways that could be used to evaluate which aligner is best?

