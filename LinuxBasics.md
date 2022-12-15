# Linux Basics

## Command Line

| Command           | Options | Example                           | Description                                               |
|-------------------|---------|-----------------------------------|-----------------------------------------------------------|
| passwd            |         |                                   | Change your password                                      |
| pwd               |         |                                   | print working directory                                   |
| ls                |         | ls /home                          | list content of directory                                 |
|                   | -l      | ls -l .                           | long; additional information for files                    |
|                   | -a      |                                   | all; also show hidden file                                |
|                   | -t      |                                   | sort by time                                              |
|                   | -r      | ls -ltr .                         | sort in reversed order                                    |
|                   | -h      | ls -h .                           | human readable                                            |
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
| grep              |         |                                   | 'grep' Content from file                                  |
| unzip FILE        |         | unzip fastqc_v0.11.9.zip          | unzip a zipped file                                       |
| wc                |         | wc file.txt                       | count the number of lines, words and bytes                |
| cut               |         |                                   | get columns from rows                                     |
| sort              |         |                                   | sort a file alphanumerically                              |
| uniq              |         |                                   | merge identical lines (which are next to each other)      |
|                   | -c      | uniq -c                           | count number of merged lines                              |
|                   |         |                                   |                                                           |


## Access rights

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


## Pipes

- Leite den Output eines Programs in ein anderes. `program1 | program1`.
- Leite den Output eines Programs in eine neue Datei um. `program1 > neueDatei.txt`

## Interaction with cloud machine

You have to be in the university network, either physically or via the vpn client.

- **login**: `ssh yourCloudName@ipadress`, e.g. `ssh joesch@10.106.241.119`
- **move file from local machine to cloud**: `scp localFile yourCloudName@ipadress:targetDir`, e.g. `scp input.txt joesch@10.106.241.119:/master/home/joesch/projects/carnivores`
- **move file from cloud to local machine** (you are on your local machine): `scp yourCloudName@ipadress:filePath localFile`, e.g. `scp joesch@10.106.241.119:/master/home/joesch/projects/carnivores/input.txt .` 
