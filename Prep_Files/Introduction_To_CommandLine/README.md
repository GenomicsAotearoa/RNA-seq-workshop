# Introduction to the command line

---

### Objectives

* Navigating your file system
* Copying, Moving, Renaming and Removing files
* Examining file contents
* Redirection, manipulation and extraction
* Loops
* Shell scripting

### Navigating your file system

Find the current location by running the command `pwd` which stands for *print working directory*. At any given time, **current working directory** is current default directory.

```bash
$ cd
$ pwd
$ /home/UserName/
```

We can see what files and subdirectories are in this directory by running `ls`, which stands for "listing".

```bash
$ ls
```

Navigating to the `MGSS_Intro/` directory can be done with `cd` command which stands for *change directory*.

```bash
$ cd RNA_seq/
```

Run the `ls` command to list the content of the current directory. Check whether there is a directory named **Intro**. Go inside this directory and run `ls` command. There should 2 x *.fastq* files



```bash
$ cd Intro

$ ls 

SRR1.fastq  SRR2.fastq
```
The `mkdir` command (*make directory*) is used to make a directory. Enter `mkdir` followed by a space, then the directory name you want to create

```bash
$ mkdir backup
```
### Copying, Moving, Renaming and Removing files

Make a second copy of `SRR1.fastq` and rename it as `backup1.fastq`.  Then move that file to `backup/` directory.

```bash
$ cp SRR1.fastq backup1.fastq

$ mv backup1.fastq backup
```

Take a look at the content of `backup` directory with `ls` command.

```bash
$ ls backup/
```

See whether you can remove the `backup/` directory by using the `rm` command . 

```bash
$ rm backup/
# rm : can not remove 'backup/': Is a directory
```

By default, `rm` will not delete directories. This can be done by using `-r` (recursive) option.

```bash
$ rm -r backup
```

### Examining file contents

There are a number of ways to examine the content of a file. `cat` and `less` are two commonly used programs for a quick look. Check the content of `SRR097977.fastq` by using these commands. Take a note of the differences. 

```
$ cat SRR1.fastq
$ less SRR1.fastq
```

A few useful shortcuts for navigating in `less`

<img src="./img/ex1_less_shortcuts.jpg" alt="drawing" width="400"/></p>

There are ways to take a look at parts of a file. For an example the `head` and `tail` commands will scan the beginning and end of a file, respectively. 

```bash
$ head SRR1.fastq

$ tail SRR1.fastq
```

Adding `-n` option to either of these commands will print the first or last *n* lines of a file.

```bash
$ head -n 1 SRR1.fastq
# @SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
```
### Redirection, manipulation and extraction

Although using `cat` and `less` commands will allow us to view the content of the whole file, most of the time we are in search of particular characters (strings) of interest, rather than the full content of the file. One of the most commonly used command-line utilities to search for strings is `grep`. Let's use this command to search for the string `NNNNNNNNNN` in `SRR2.fastq` file.

```bash
$ grep NNNNNNNNNN SRR2.fastq
```

Retrieve and discuss the output you get when `grep` was executed with the `-B1` and `-A1` flags.

```bash
$ grep -B1 -A2 NNNNNNNNNN SRR2.fastq
```

In both occasions, outputs were printed to the terminal where they can not be reproduced without the execution of the same command. In order for "string" of interest to be used for other operations, this has to be "redirected" (captured and written into a file). The command for redirecting output to a file is `>`. Redirecting the string of bad reads that was searched using the `grep` command to a file named `bad_reads.txt` can be done with

```bash
$ grep -B1 -A2 NNNNNNNNNN SRR2.fastq > bad_reads.txt
```

Use the `wc` command to count the number of words, lines and characters in the `bad_reads.txt` file.

```bash
$ wc bad_reads.txt
```

Add `-l` flag to `wc` command and compare the number with the above output

```bash
$ wc -l bad_reads.txt
```

#### basename

`basename` is a function in UNIX that is helpful for removing a uniform part of a name from a list of files. In this case, we will use `basename` to remove the .fastq extension from the files that we've been working with.

```bash
$ basename SRR1.fastq .fastq
```
### Loops

Loops are a common concept in most programming languages which allow us to execute commands repeatedly with ease. There are three basic loop constructs in `bash` scripting,

* **for** - iterates over a list of items and performs the given set of commands

```
for item in [LIST]
do
    [COMMANDS]
done
```
For most of our uses, a `for loop` is sufficient for our needs, so that is what we will be focusing on for this exercise.

Shell identifies the `for` command and repeats a block of commands once for each item in a list. The for loop will take each item in the list (in order, one after the other), assign that item as the value of a variable, execute the commands between the `do` and `done` keywords, then proceed to the next item in the list and repeat over. The value of a variable is accessed by placing the `$` character in front of the variable name. This will tell the interpreter to access the data stored within the variable, rather than the variable name. For example

```bash
$ i="rna-seq"

$ echo i
# i
$ echo $i
# rna-seq
$ echo ${i}
# rna-seq
```

This prevents the shell interpreter from treating `i` as a string or a command. The process is known as *expanding* the variable. We will now wrtite a for loop to print the first two lines of our *fastQ* files:

```
$ for filename in *.fastq
>do
>    head -n 2 ${filename}
>done
```

`basename` is rather a powerful tool when used in a for loop. It enables the user to access just the file prefix which can be use to name things

```
$ for filename in *.fastq
>do
>    name=$(basename ${filename} .fastq)
>    echo ${name}
>done
```

### Scripts

Executing operations that contain multiple lines/tasks or steps such as for loops via command line is rather inconvenient. For an example, imagine fixing a simple spelling mistake made somethwhere in the middle of a for loop that was directly executed on the terminal.

The solution for this is the use of shell scripts, which are essentially a set of commands that you write into a text file and then run as a single command. In UNIX-like operating systems, in built text editors such as `nano`, `emacs`, and `vi` provide the platforms to write scripts. For this workshop we will use `nano` to create a file named `ForLoop.sh`.

```
$ nano forLoop.sh
```

Add the following for loop to the script (note the header `#!/bin/bash`).

```bash
#!/bin/bash

for filename in *.fastq
do
    head -n 2 ${filename}
done
```

Because `nano` is designed to work without a mouse for input, all commands you pass into the editor are done via keyboard shortcuts. You can save your changes by pressing `Ctrl + O`, then exit `nano` using `Ctrl + X`. If you try to exit without saving changes, you will get a prompt confirming whether or not you want to save before exiting, just like you would if you were working in **Notepad** or **Word**.

Now that you have saved your file, see if you can run the file by just typing the name of it (as you would for any command run off the terminal). You will notice the command written in the file will not be executed. The solution for this is to tell the machine what program to use to run the script. 

```bash
$ bash ForLoop.sh
```