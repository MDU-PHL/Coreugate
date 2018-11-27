[![Build Status](https://travis-ci.org/kristyhoran/Coreugate.svg?branch=master)](https://travis-ci.org/kristyhoran/Coreugate)

# COREugate - A pipeline for cgMLST
## From reads to cgMLST profile.

This is a simple pipeline that allows the user to input paired-end reads or assemblies and using a user-defined cgMLST schema and will output a cgMLST profile for the isolates as well as statistics for allele calling and a matrix of pairwise allelic distances.

1. Assemble - default using [Shovill](https://github.com/tseemann/shovill) (implementing the latest version of Spades)
2. PrepSchema (if necessary) and Call alleles using [chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki).
3. Combine profiles and statisitics for the whole dataset.
4. Calculate pairwise allelic distances (missing data is ignored)

### Dependencies
```
Python >=3.6
Biopython >=1.70
Snakemake >=5.3.0
chewBBAC >=2.0.16
```
COREugate uses Shovill (with SPAdes, Skesa or Velvet or the assemblers alone) to perform assemblies. There is a singularity container of these tools available. If you do not want to use the singularity container you will need to ensure installation of these assemblers. The easiest way to do so is to install `shovill`. Detailed intructions can be found [here](https://github.com/tseemann/shovill).

```
brew install brewsci/bio/shovill
```

OR 

```
conda install -c bioconda shovill
```

### Biopython
Biopython is used here for quality control of assemblies, information about biopython can be found [here](https://biopython.org)
```
pip3 install biopython
```

### Snakemake
Ensure that you have Snakemake installed. Detailed instructions can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

```
pip3 install snakemake
```

### chewBBACA
[chewBBACA](https://github.com/B-UMMI/chewBBACA/) is used here to prepare the schema, by selecting exemplar alleles for comparison and to call allele profiles. More information about chewBBACA and how it is works can be found [here](https://github.com/B-UMMI/chewBBACA/wiki). COREugate can use a singularity version of chewBBACA, however if you want to install the latest version (>=2.0.16)

```
pip3 install chewbbaca
```

### Run COREugate

#### Get COREugate
```
pip3 install git+https://github.com/kristyhoran/Coreugate
```

If you are installing COREugate on a server using `--user` please ensure that your `~/.local/bin` is part of your PATH
```
export PATH=$PATH:/path/to/.local/bin
```

#### Running COREugate
```
coreugate run -h

Usage:
  run [options]

Options:
  -f, --input_file=INPUT_FILE                Input file - three columns <isolate_id> <path/to/R1> <path/to/R2> [default: ""]
  -d, --schema=SCHEMA                        Path to species schema/allele db [default: ""]
  -w, --workdir=WORKDIR                      Working directory, defaults to CWD [default: "/home/khhor/cgMLST/listeria/test"]
  -r, --template_dir=TEMPLATE_DIR            Path to Snakefile and config.yaml templates. [default: "/home/khhor/dev/Coreugate/templates"]
  -t, --threads=THREADS                      Number of threads for chewBBACA [default: 36]
  -i, --id=ID                                Name of the job [default: ""]
  -s, --singularity=SINGULARITY              Use singularity containers for assembly and chewBBACA rather than local install [default: "N"]
  -c, --min_contig_size=MIN_CONTIG_SIZE      Minumum contig size required for QC [default: 500]
  -m, --min_contigs=MIN_CONTIGS              Minumum number of contigs required for QC [default: 0]
  -a, --assembler=ASSEMBLER                  Assembler to be used (options are: shovill-spades, shovil-skesa, shovill-velvet, skesa, spades) [default: "shovill-spades"]
  -p, --prodigal_training=PRODIGAL_TRAINING  Prodigal file to be used in allele calling. See https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files for options [default: ""]
  -h, --help                                 Display this help message
```
**Input**
You can run COREugate with paired-end reads or assemblies. For the sake of consistency it is recomended that if you are using assemblies you use assemblies constructed with the same tool, since assembler and assembler settings can greatly impact allele calling.

###### Sample data
Sample/isolate data is in the form of a tab-delimited file containing a sample identifier and the path to reads or assemblies

*Reads*
```
isolate_name	path/to/reads/R1.fq.gz	path/to/reads/R2.fq.gz
```

*Assemblies*
```
isolate_name	path/to/assembly.fa	
```
###### Species cgMLST schema
COREugate requires an exisiting cgMLST schema, this can be a schema generated by the user or downloaded from one of the publically available databases. These schema should be in the format of a `fasta` file for each loci, each file should contain the different alleles for each loci. It should be noted that during allele calling, chewBBACA (implemented by COREugate) will add inferred alleles ([more information](https://github.com/B-UMMI/chewBBACA/wiki)) to your schema, so it is recommended that the schema path be fixed, that is that the schema is kept in a central location and a single version is used for each species/study.

###### Other optional arguments
* `job_id` this is compulsory and should be unique
* `min_contig_size` defaults to 500bp, but can be user-defined.
* `min_contigs`, user defined, although defaults to no minimum. It is recommended to determine a minimum number of contigs, since low number of contigs improves allele calling.
* `assembler` defaults to shovill implementation of SPAdes, other options are `spades, skesa, shovill-skesa ` and `shovill-velvet`.
* `singularity` defaults to system software versions (no singularity). To use singularity containers use `--singularity/-s Y`
* `prodigal_training` a prodigal training file for allele calling. Recommended by chewBBACA developers, a list of default training files and further information can be found [here](https://github.com/B-UMMI/chewBBACA/wiki).


**Example**
Run with singularity and prodigal training file
`coreugate run -f test.tab -d path/to/schema -i job_id -s y -p species_prodigal_training.trn`

### Limitations of the pipeline
* Coreugate is only able to work with pre-exisiting schemas that have been prep as described above, to derive profiles for isolates.
* Possibly more, I just haven't found them yet!!

