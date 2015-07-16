===================
pMLST
===================

This project documents MLST service


Documentation
=============

## What is it?

The MLST service contains one perl script *pMLST-1.4.pl* which is the script of the lates
version of the pMLST service. The method enables investigators to determine the ST based on WGS data.

## Usage

To use the service some data needs to be pre-installed: *database*, *blast-2.2.26* and *Makefile*

The folder *database* includes all the MLST and pMLST schemes and needs to be updataed to get the best results.
The datasets are extracted from the http://pubmlst.org/ webside weakly and can be downloaded from
http://cge.cbs.dtu.dk/services/data.php. 

The folder *blast-2.2.26* includes blastall and formatdb which are used by the *pMLST-1.4.pl* script

The file *Makefile* installs the nesserary perl modules to run the *pMLST-1.4.pl* script. It is used by writing:
    make install

The program can be invoked with the -h option to get help and more information of the service.

```bash
Usage: perl pMLST-1.4.pl [options]

Options:

    -h HELP
                    Prints a message with options and information to the screen
    -d DATABASE
                    The path to where you have located the database folder
    -b BLAST
                    The path to the location of blast-2.2.26
    -i INFILE
                    Your input file which needs to be preassembled partial or complete genomes in fasta format
    -o OUTFOLDER
                    The folder you want to have your output files places
    -s SPECIES
                    The pMLST scheme you want to use. The options can be found in the *pmlst_schemes* file
```

## Example of use with the *database* and *blast-2.2.26* folder loacted in the current directory
    
    perl pMLST-1.4.pl -i test.fsa -o OUTFOLDER -s incf 

## Example of use with the *database* and *blast-2.2.26* folder loacted in antoher directory

    perl pMLST-1.4.pl -d path/to/database -b path/to/blast-2.2.26 -i INFILE.fasta -o OUTFOLDER -s incf 
    

## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/pMLST/


## The Latest Version


The latest version can be found at
https://bitbucket.org/genomicepidemiology/pmlst/overview

## Documentation


The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/pmlst/overview.

Installation
=======

The scripts are self contained. You just have to copy them to where they should
be used. Only the *database* folder needs to be updataed mannually. 

Remember to add the program to your system path if you want to be able to invoke the program without calling the full path.
If you don't do that you have to write the full path to the program when using it.

Citation
=======

When using the method please cite:

PlasmidFinder and pMLST: in silico detection and typing of plasmids.
Carattoli A, Zankari E, Garcia-Fernandez A, Volby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.
Antimicrob. Agents Chemother. 2014. April 28th.
[Epub ahead of print]


License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
