pMLST (Plasmid Multi-Locus Sequence Typing)
===================

pMLST determines plasmid sequence types from assembled genomes or sequencing
reads using pMLST schemes from the CGE pMLST database.


## Installation

```bash
# Bioconda is preferred; it installs pMLST, BLAST, and KMA.
conda install -c conda-forge -c bioconda pmlst_ssi

# PyPI installs the Python package only. Install blastn/kma separately.
python3 -m pip install pmlst_ssi

# Source install also requires blastn/kma from conda or the system.
git clone https://github.com/ssi-dk/pmlst.git
cd pmlst
python3 -m pip install .
```

## Database

Install or update the pMLST database after installing the tool. With conda, the
default destination is the current environment's pMLST database location.

```bash
# Default location: conda env DB if conda is active, otherwise user data dir.
pmlst-download-db

# Custom location.
pmlst-download-db /path/to/pmlst_db

# Use a custom DB without passing -p every time.
export PMLST_DB=/path/to/pmlst_db
```

pMLST looks for the database in this order:

1. `-p` / `--database`
2. `PMLST_DB`
3. `$CONDA_PREFIX/share/pmlst/db`
4. the platform-specific user data directory from platformdirs

## Docker

The Docker image includes pMLST, runtime dependencies, and a bundled database.

```bash
docker build -t pmlst .
docker run --rm -v "$PWD:/workdir" pmlst pmlst -i /workdir/test_data/test.fsa -s incf -x
```

Rebuild the image to update the bundled database:

```bash
docker build --no-cache -t pmlst .
```

## Usage

Common command examples:

```bash
# Run pMLST on assembled contigs with BLAST.
pmlst -i test_data/test.fsa -s incf -x

# Run with an explicit method executable.
pmlst -i test_data/test.fsa -s incf -p /path/to/pmlst_db -mp /usr/bin/blastn

# Run with custom output and temporary directories.
pmlst -i sample.fasta -s incf -p /path/to/pmlst_db -o results -t tmp -x
```

Common command-line options:

```text
Required:
  -i INPUTFILE [INPUTFILE ...]   One or more FASTA or FASTQ input files.
  -s SCHEME                      Scheme, comma-separated schemes, or all.
  -pf PF_RESULTS                 Precomputed PlasmidFinder results.

  Use either -s or -pf.
  Available schemes: incac, incf, inchi1, inchi2, inci1, incn,
                     pbssb1-family, shigella.

Output options:
  -o OUTDIR                      Output directory.
  -x                             Write extended text and FASTA outputs.
  -xm                            Like -x, but prefix extended output filenames
                                 with the scheme name for multi-scheme runs.
  -so                            Write simple tab-separated allele profile output.
  -q                             Quiet command-line output.
  -t TMP_DIR                     Temporary directory for external tool output.

Other options:
  -p DATABASE                    pMLST database directory.
  -mp METHOD_PATH                Path to blastn or kma.
  -c COVERAGE                    Minimum coverage.
  -id IDENTITY                   Minimum identity.
```

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
