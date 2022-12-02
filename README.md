# Bioinformatics Master

This repository contains additional material for the introductory Bioinformatics course taught as part of the Master Life Sciences at the university of Würzburg. 

## Linux Basics

### Command Line

| Command         | Options | Example                       | Description                                               |
|-----------------|---------|-------------------------------|-----------------------------------------------------------|
| pwd             |         |                               | print working directory                                   |
| ls              |         | ls /home                      | list content of directory                                 |
|                 | -l      |                               | long; additional information for files                    |
|                 | -a      |                               | all; also show hidden file                                |
|                 | -t      |                               | sort by time                                              |
|                 | -r      |                               | sort in reversed order                                    |
| cd              |         | cd /master/data               | change into another directory                             |
|                 |         | cd                            | without directory, change to home directory               |
| mkdir DIRECTORY |         | mkdir results                 | make a new directory                                      |
| rmdir DIRECTORY |         | rmdir testDir                 | delete an empty (!) directory                             |
| touch FILE      |         | touch newFile.txt             | if file does not exist, generate it (rarely used...)      |
| cp FROM TO      |         | cp thesis.txt backup.txt      | copy a file                                               |
| mv FROM TO      |         | mv experiment.txt newName.txt | move a file, i.e. either rename it or move to a new place |
|                 |         | mv experiment.txt ..          |                                                           |
| less FILE       |         | less results.txt              | show contents of file (type 'q' to quit)                  |
| head FILE       |         | head reads.fastq              | show first 10 lines of a file                             |
|                 | -number | head -100 reads.fastq         | show the first 'number' lines of file                     |
| tail FILE       |         | tail reads.fastq              | show the last 10 lines of a file                          |
|                 | -number | tail -20 reads.fastq          | show the last 'number' lines of a file                    |
