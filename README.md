# tRFtarget-pipeline
This repository holds pipeline used in **tRFtarget database** (http://trftarget.net/). It can be used to find interaction sites on target transcripts by [*RNAhybrid*](https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid) and [*IntaRNA*](http://rna.informatik.uni-freiburg.de/IntaRNA/Input.jsp) for transfer RNA-derived fragments (tRFs) which are not indexed in the database . It can also be used to find target genes for other small RNAs such as miRNAs

## Enclosed Package Version

* ***RNAhybrid***: 2.1.2
* ***IntaRNA***: 3.1.3 with *Vienna RNA* 2.4.14 and *boost* 1.70.0

Please refer the [Prediction Tools Setting](http://trftarget.net/method) section in tRFtarget database for detailed description of option settings of *RNAhybrid* and *IntaRNA* in this tRFtarget-pipeline, but generally speaking, both prediction tools are tuned to provide binding sites with different prediction mechanism:

* ***RNAhybrid***: **minimum free energy**
* ***IntaRNA***: **minimum free energy** + **seed match** + **accessibility** feature

## Installation

We provided **Docker** or **Singularity** images for immediate using of tRFtarget-pipeline

```bash
# For Docker
docker pull az7jh2/trftarget:0.2.1
# For Singularity
singularity build trftarget-0.2.1.sif docker://az7jh2/trftarget:0.2.1
```

To test the installation (should print the version of tRFtarget-pipeline)

```bash
# For Docker
docker run -it --rm az7jh2/trftarget:0.2.1 tRFtarget -v
# For Singularity
singularity exec trftarget-0.2.1.sif tRFtarget -v
```

## Usage

### Command

The command to run tRFtarget-pipeline with default setting is:

```bash
# For Docker
docker run -it --rm -v <path>:/data az7jh2/trftarget:0.2.1 tRFtarget -q <query_fasta_file_name> -t <target_fasta_file_name> -n 1 --e_rnahybrid -15 --e_intarna 0 -b 1 -s 6
# For Singularity
singularity exec -B <path>:/data trftarget-0.2.1.sif tRFtarget -q <query_fasta_file_name> -t <target_fasta_file_name> -n 1 --e_rnahybrid -15 --e_intarna 0 -b 1 -s 6
```

`<path>` is the valid and **absolute** path of the folder in the host machine to be mounted in the Docker/Singularity image for data exchanging (`readlink -f` can be used to get the absolute path of one folder)

 `<query_fasta_file_name>` and `<target_fasta_file_name>` are the file names (**Without path**) of *FASTA* files of query small RNAs and target RNAs respectively. both of them are required to be located in the `<path>` folder

### Input

* A *FASTA* file of query small RNAs. Compressed file (such as `.gz`) currently not supported
* A *FASTA* file of target RNAs (optional). Compressed file (such as `.gz`) currently not supported. If not provided, use 100,218 Protein-coding transcript sequences (GRCh38.p13) as target RNAs instead

### Output

The output of tRFtarget-pipeline are 6 *CSV* files located in the `<path>` folder:

1. `trfs_info.csv` : show tRF ID, sequence and sequence length
2. `transcripts_info.csv` : show transcript ID, sequence and length
3. `rnahybrid_results.csv` : predicted RNA-RNA interactions by *RNAHybrid*
4. `intarna_results.csv` : predicted RNA-RNA interactions by *IntaRNA*
5. `consensus_results.csv` : a consensus predictions between *RNAHybrid* and *IntaRNA*. The definition of consensus can be referred in the [Definition of consensus entries in *RNAhybrid* and *IntaRNA* predictions](http://trftarget.net/manual) section in tRFtarget database
6. `tRF_level_consensus_stats.csv` : a summary of numbers of interactions predicted by *RNAHybrid* and *IntaRNA*, as well as the number of consensus predictions. It also includes the percentage of consensus predictions in *RNAHybrid* and *IntaRNA* predictions, respectively

### Binding sites in *CSV* files

The *CSV* files containing predicted binding sites (`rnahybrid_results.csv`, `intarna_results.csv` and `consensus_results.csv`) have the unified format. The total 14 columns are shown as below:

| Column          | Description                                                  |
| --------------- | ------------------------------------------------------------ |
| `tRF_ID`        | ID of query sequence, corresponding to the sequence ID in *FASTA* file of query small RNAs. |
| `Transcript_ID` | ID of target sequence, corresponding to the sequence ID in *FASTA* file of target RNAs. |
| `Demo`          | An ASCII RNA-RNA interaction illustration.                   |
| `Start_tRF`     | Start index of RNA hybrid in query sequence. Index starts at 5', and index number starts from 1. |
| `End_tRF`       | End index of RNA hybrid in query sequence. Index starts at 5', and index number starts from 1. |
| `Start_Target`  | Start index of RNA hybrid in target sequence. Index starts at 5', and index number starts from 1. |
| `End_Target`    | End index of RNA hybrid in target sequence. Index starts at 5', and index number starts from 1. |
| `MFE`           | Calculated Free Energy of that binding site.                 |
| `HybridDP`      | VRNA dot-bracket notation for RNA hybrid (interaction sites only). |
| `SubseqDP`      | Hybrid subsequences compatible with `HybridDP`.              |
| `Max_Hit_DP`    | Hybrid subsequences in maximum complementary region. Please refer the [Definition of maximum complementary length](http://trftarget.net/manual) section in tRFtarget database for detailed description of maximum complementary region and maximum complementary length. |
| `Max_Hit_Len`   | Sequence length of maximum complementary region. Please refer the [Definition of maximum complementary length](http://trftarget.net/manual) section in tRFtarget database for detailed description of maximum complementary region and maximum complementary length. |
| `Tool`          | Tool used to predict that binding site (*RNAhybrid* or *IntaRNA*). |
| `Consensus`     | Indicate the consensus entries (`=1`). Non-consensus entries are labelled as `0`. |

### All Options

| Option                 | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| `-q` or`--query`       | *FASTA* file of query small RNAs. Required.                  |
| `-t` or `--target`     | *FASTA* file of target RNAs. If not provided, use 100,218 **Human** Protein-coding transcript sequences (GRCh38.p13) as target RNAs. |
| `-n` or `--n_cores`    | Number of CPU cores used for parallel computing. Default value is 1 (no parallel computing). |
| `--e_rnahybrid`        | Free energy threshold for *RNAhybrid*, used for *RNAhybrid* `-e` option. Default value is -15. |
| `--e_intarna`          | Free energy threshold for *IntaRNA*, used for *IntaRNA* `--outMaxE` option. Default value is 0. |
| `-b` or `--suboptimal` | Reported number of interaction sites on each target RNA, used for *RNAhybrid* `-b` option and *IntaRNA* `-n` option. Default value is 1. |
| `-s` or `--seed_len`   | For *RNAhybrid*, threshold of maximum complementary length interactions with maximum complementary length less than it are filtered out. <br/>For *IntaRNA*, threshold of the number of base pairs within the seed sequences, used for *IntaRNA* `-seedBP` option.<br/>Default value is 6 |

### Elapsed Time & Output File Size

Take **1** tRF (*tRF-1001*) and the default target RNAs (**100,218** Protein-coding transcript sequences) for example. All options are leaving as default (No parallel computing)

Elapsed time for whole pipeline: 48.26 hours

* running *RNAhybrid*: 0.72 hours
* running *IntaRNA*: 47.50 hours
* consensus target predictions: 0.02 hours

File size of output CSV files:

* `rnahybrid_results.csv`: 60 MB; including 90,398 target site entries
* `intarna_results.csv`: 58 MB; including 100,141 target site entries
* `consensus_results.csv`: 28 MB; including 22,492 *RNAhybrid* entries and 22,492 *IntaRNA* entries

It's recommended to **turn on the parallel computing by specifying `-n` or `--n_cores` option, which will significantly reduce the running time of *IntaRNA***

## Citation

If you use tRFtarget-pipeline, please cite:

Ningshan Li, Nayang Shan, Lingeng Lu, Zuoheng Wang. tRFtarget: a database for transfer RNA-derived fragment targets. *Nucleic Acids Research*, 2020, gkaa831. [https://doi.org/10.1093/nar/gkaa831](https://doi.org/10.1093/nar/gkaa831)

## Contact Us

Users are welcome to send feedbacks, suggestions or comments related to the **tRFtarget database** through our GitHub repository [tRFtarget](https://github.com/ZWang-Lab/tRFtarget).

For issues in using **tRFtarget-pipeline**, please report to this GitHub repository [tRFtarget-pipeline](https://github.com/ZWang-Lab/tRFtarget-pipeline).
