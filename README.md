# Bioinformatics Master

## Goal

This repository contains additional material for the introductory Bioinformatics course taught as part of the Master Life Sciences at the university of Würzburg. In this course, you will identify genes which are up-regulated in a digesting Venus Flytrap Trap and characterize their function and evolution. To reach this goal, you will learn how to work on a remote linux machine and apply a typical RNASeq workflow. We will use RNASeq data from [Bemm et al.](https://genome.cshlp.org/content/26/6/812.full.html) and the genome published in [Palvalfi et al.](https://www.sciencedirect.com/science/article/pii/S0960982220305674)

## Linux Basics

### Command Line

| Command           | Options | Example                           | Description                                               |
|-------------------|---------|-----------------------------------|-----------------------------------------------------------|
| passwd            |         |                                   | Change your password                                      |
| pwd               |         |                                   | print working directory                                   |
| ls                |         | ls /home                          | list content of directory                                 |
|                   | -l      | ls -l .                           | long; additional information for files                    |
|                   | -a      |                                   | all; also show hidden file                                |
|                   | -t      |                                   | sort by time                                              |
|                   | -r      | ls -ltr .                         | sort in reversed order                                    |
| cd                |         | cd /master/data                   | change into another directory                             |
|                   |         | cd                                | without directory, change to home directory               |
| mkdir DIRECTORY   |         | mkdir results                     | make a new directory                                      |
| rmdir DIRECTORY   |         | rmdir testDir                     | delete an empty (!) directory                             |
| touch FILE        |         | touch newFile.txt                 | if file does not exist, generate it (rarely used...)      |
| cp FROM TO        |         | cp thesis.txt backup.txt          | copy a file                                               |
| mv FROM TO        |         | mv experiment.txt newName.txt     | move a file, i.e. either rename it or move to a new place |
|                   |         | mv experiment.txt ..              |                                                           |
| less FILE         |         | less results.txt                  | show contents of file (type 'q' to quit)                  |
| head FILE         |         | head reads.fastq                  | show first 10 lines of a file                             |
|                   | -number | head -100 reads.fastq             | show the first 'number' lines of file                     |
| tail FILE         |         | tail reads.fastq                  | show the last 10 lines of a file                          |
|                   | -number | tail -20 reads.fastq              | show the last 'number' lines of a file                    |
| whoami            |         |                                   | OK, the name says it all                                  |
| groups            |         |                                   | to which groups do I belong?                              |
| chmod RIGHTS FILE |         | chmod go-r myThesis.txt           | change the access rights of a file                        |
| wget URL          |         | wget https://www.uni-wuerzburg.de | get the content of an Url                                 |
| unzip FILE        |         | unzip fastqc_v0.11.9.zip          | unzip a ziped file                                        |

### Access rights

As we are working on a multi-user system, we want to clearly define who can do what with our files. In Linux, there are three levels of access rights - (u)ser, (g)roup and (o)thers. If you want to know to which groups you belong, type `group` on the command line. In addition, there are three things a user can do with a file - (r)ead, (w)rite and e(x)ecute. For each of the three levels, we can individually set the three things that can be done.

| (u)ser                  | (g)roup | (o)thers |
|-------------------------|---------|----------|
| (r)read(w)ritee(x)ecute | rwx     | rwx      |

If you type `ls -l` in a directory, you see the current rights for all files in the directory. To change the right, use the `chmod` command. This takes as first option the rights you want to change and as second option the file name. To grant rights use `+`, to revoke them `-`. This also works for directories.

| command               | description                                                           |
|-----------------------|-----------------------------------------------------------------------|
| chmod u+x fastqc      | make the file fastqc executable for the user                          |
| chmod go-r thesis.txt | only the user can read the thesis (revoke read from group and others) |
| chmod o-x results     | assuming that results is a directory: others can't cd into results    |

### Interaction with cloud machine

You have to be in the university network, either physically or via the vpn client.

- **login**: `ssh yourCloudName@ipadress`, e.g. `ssh joesch@10.106.241.119`
- **move file from local machine to cloud**: `scp localFile yourCloudName@ipadress:targetDir`, e.g. `scp input.txt joesch@10.106.241.119:/master/home/joesch/projects/carnivores`
- **move file from cloud to local machine** (you are on your local machine): `scp yourCloudName@ipadress:filePath localFile`, e.g. `scp joesch@10.106.241.119:/master/home/joesch/projects/carnivores/input.txt .` 

## RNASeq Analysis

### Read Quality control

#### install FastQC

Go to the [FastQC Web page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and click on **Download Now**. We need the `FastQC v0.11.9 (Win/Linux zip file)` (version might have changed!). If you are using an interactive browser, the file will be stored on your local computer. Doesn't help if you want to run the program on the Cloud machine :D. There are two solutions for this:
1. Download it locally and use `scp` to [move the file to our cloud machine](#interaction-with-cloud-machine) (wastes storage on your local machine and takes time)
2. Use the command line tool `wget` on the cloud machine. To get the Url, right-click on the file name and "copy link address". on the cloud command line, you can now paste this with your right mouse button. You should get `wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip`

Next, unzip the file with `unzip fastqc_v0.11.9.zip`. Move into the new directory ( `cd FastQC/` ) and make the program file executable, i.e. give yourself execute rights on this file - `chmod u+x fastqc`.

#### run FastQC

As the fastqc program is not in your path (try `echo $PATH` to see your current path), you have to tell Linux, where to find it. So, if you are in the directory where you unpacked FastQC this is `./fastqc`. As the [documentation](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt) tells us, fastqc will write its results into the folder where the sequence files are. This is not only ugly (we mix data and results) but in our setup not even possible, as we don't have write access to the directory which contains our RNASeq data. So, we have to use the `--outdir=` parameter to write the results in a different directory. For the moment, let's write the results into the current directory:
```
./fastqc --outdir=. /master/home/data/transcriptome/DM_exp001_Tr_L1_P1.fastq
```

#### check results

FastQC will generate two result files, an `.html` and a `.zip` file. To view the result, [move](#interaction-with-cloud-machine) the `.zip` file to your local machine. Open it in the web browser of your choice. For me, the most important plots are the [Per base sequence quality](DM_exp001_Tr_L1_P1_fastqc.html#M1) and the [Per sequence GC content](DM_exp001_Tr_L1_P1_fastqc.html#M5). Here are some questions to think about:

- What can we do if the sequence quality drops at the end of the reads?
- How would the GC content plot look, if we have a massive contamination like bacteria?

### Read Mapping

#### Install STAR

Check the STAR [GitHub repository](https://github.com/alexdobin/STAR). You can clone this repository to your local computer using the URL given by the green '<> Code' button. Make sure you selected 'https' and copy the URL. Now, in a directory of your choice (may I suggest `src`;-), clone the git repository.
```
git clone https://github.com/alexdobin/STAR.git
```
(Just as a note: You can clone this git repository just the same way)

The manual (check the GitHub start page or the README.md file in your local git repository) tells us, that we already have pre-compiled binaries in the bin directory and that the `static` executables are the easiest to use. Let's check this:

```
cd STAR/bin
ls
cd Linux_x86_64_static
./STAR
```

In case you are a freak (in this course this is an honor - same for geek and nerd), you might want to compile the source code yourself :o . To do this, we follow the documentation and 'run `make` in  the source directory'. This will take some time. If you want to get some idea about what's going on here, there is an [awesome C++ course on YouTube](https://www.youtube.com/playlist?list=PLlrATfBNZ98dudnM48yfGUldqGD0S4FFb). Most relevant here is the [How the C++ Works](https://www.youtube.com/watch?v=SfGuIVzE_Os&list=PLlrATfBNZ98dudnM48yfGUldqGD0S4FFb&index=5). And of course, if you want to learn how to program C++, just work through the course (don't :D). But I digress....

#### Running STAR

##### Generating Genome Index
The manual tells us, that the basic STAR workflow consists of two steps - generating a genome index and mapping the reads. So, let's start with generating a genome index. Just as a reminder, the genome fasta file is in `/master/home/data/genome/Dm_genome_assembly.fa`. Here, you can also find the annotation of the genome (`Dm_annotation.gff`). Be aware that this is not a **GTF** but a **GFF** file! Check the manual how to handle these! Also, remember that we have a big genome (3 gigaBases) and the genome assembly is very fragmented with about 100.000 scaffolds. For the number of threads, remember that we have 16 cores and there might be more than one user....
<!--
```
~/src/STAR/STAR/bin/Linux_x86_64_static/STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ./genome_index --genomeFastaFiles /master/home/data/genome/Dm_genome_assembly.fa --sjdbGTFfile /master/home/data/genome/Dm_annotation.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --genomeChrBinNbits 15
```
-->

##### Mapping Reads to Genome

Now that we have generated the genome index, we can mapp the reads onto the genome. Remember that we have paired-end reads, i.e. two read files per experiment. This step will take some time, so maybe only map one experiment right now. Plus, you might want to take a smaller experiment ;-).

<!--
```
~/src/STAR/STAR/bin/Linux_x86_64_static/STAR --runThreadN 8 --genomeDir ./genome_index --readFilesIn /master/home/data/transcriptome/DM_exp001_Tr_L1_P1.fastq /master/home/data/transcriptome/DM_exp001_Tr_L1_P2.fastq
```
-->