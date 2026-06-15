pMLST
===================

Plasmid Multi-Locus Sequence Typing

pMLST determines plasmid sequence types from assembled genomes or sequencing
reads using pMLST schemes from the CGE pMLST database.


## Installation

Bioconda installation is preferred because it installs pMLST together with
BLAST and KMA runtime dependencies:

```bash
conda install -c conda-forge -c bioconda pmlst
```

Install or update the pMLST database after installing the tool. With conda,
this installs into the current environment's default pMLST database location:

```bash
pmlst-download-db
```

To install the database at a specific location instead:

```bash
pmlst-download-db /path/to/pmlst_db
```

Then pass a custom database path with `-p` or export it:

```bash
export PMLST_DB=/path/to/pmlst_db
```

pMLST looks for the database in this order:

1. `-p` / `--database`
2. `PMLST_DB`
3. `$CONDA_PREFIX/share/pmlst/db`
4. the platform-specific user data directory from platformdirs

Install from PyPI. This installs the Python package only; `blastn` is required
for FASTA input and `kma` is required for FASTQ input. They must be available on
`PATH`, or their executable paths must be supplied with `-mp`.

```bash
python3 -m pip install pmlst
```

Install from source. This also installs the Python package only; it does not
install BLAST or KMA.

```bash
git clone https://bitbucket.org/genomicepidemiology/pmlst.git
cd pmlst
python3 -m pip install .
```

The Docker image includes the packaged tool and a bundled pMLST database.
Rebuild the image to fetch the current upstream database into the image:

```bash
docker build -t pmlst .
docker run --rm -v "$PWD:/workdir" pmlst pmlst -i /workdir/test_data/test.fsa -s incf -x
```

To keep an updated database outside the image, mount a host directory and run
the database downloader in the container:

```bash
docker run --rm -v "$PWD/pmlst_db:/db" pmlst pmlst-download-db /db
```

## Usage

Show command help and version:

```bash
pmlst --help
pmlst --version
pmlst-download-db --help
```

Run pMLST on assembled contigs with BLAST:

```bash
pmlst -i test_data/test.fsa -s incf -p /path/to/pmlst_db -x
```

Run with an explicit method executable:

```bash
pmlst -i test_data/test.fsa -s incf -p /path/to/pmlst_db -mp /usr/bin/blastn
```

Run with custom output and temporary directories:

```bash
pmlst -i sample.fasta -s incf -p /path/to/pmlst_db -o results -t tmp -x
```

Common options:

- `-i INPUTFILE`: FASTA or FASTQ input file.
- `-s SCHEME`: scheme name, comma-separated scheme list, or `all`.
- `-p DATABASE`: pMLST database directory.
- `-o OUTDIR`: output directory.
- `-t TMP_DIR`: temporary directory for external tool output.
- `-x`: write extended text and FASTA outputs.
- `-mp METHOD_PATH`: path to `blastn` or `kma`.
- `-c COVERAGE`: minimum coverage.
- `-id IDENTITY`: minimum identity.

## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/pMLST/

Citation
=======

When using the method please cite:

PlasmidFinder and pMLST: in silico detection and typing of plasmids.
Carattoli A, Zankari E, Garcia-Fernandez A, Volby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.
Antimicrob. Agents Chemother. 2014. April 28th.
[Epub ahead of print]

References
=======

1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics 2009; 10:421. 
2. Clausen PTLC, Aarestrup FM, Lund O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 2018; 19:307. 

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
