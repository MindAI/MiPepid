# MiPepid

`MiPepid` is a software specifically for predicting the coding capabilities of sORFs.

The corresponding paper "MiPepid: *Mi*cro*pep*tide *id*entification tool using machine learning"[*] is now submitted to BMC Bioinformatics.

sORFs / smORFs are short open reading frames with length <= 303 bp (including the stop codon), and if translated, they encode micropeptides that are <= 100 amino acids. 

Micropeptides were traditionally ignored due to their size but are now gaining increasing attentions because they have been shown to play critical roles in many vital biological activities.

### What does MiPepid do?
Given a fasta file containing DNA fasta sequences, for each sequence, `MiPepid` will find all the sORFs (length <= 303 bp) present in all the 3 translation frames of the sequence, and for each sORF it will return the predicted class label (coding or noncoding) as well as the probability of being in that class. All the results will be written in an output .csv file.  


### Dependencies: 

Language dependency:
**`Python 3`** (Please do *not* use Python 2 to run the code.)

Library dependency:
* `Numpy`
* `pandas`
* `Bio` (the Biopython package: https://biopython.org)


### How to use:
```sh
cd MiPepid
python3 ./src/mipepid.py input_fasta_file_path output_fasta_file_path
```

Note:
* The `input_fasta_file_path` is the path of your input fasta file. 
  * Please make sure it is a fasta file with only DNA sequences (not RNA sequences or protein sequences). 
  * Please make sure that each record in your fasta file has a *unique* ID, as the IDs will be used to name the found ORFs. 
  * Please also make sure that your sequences *only* contain the four letters: `A`, `T`, `C`, `G`. Other DNA letters, such as `N`, `R`, `Y`, etc. are currently *not*  supported.
* The `output_fasta_file_path` is the path of the output `.csv` file. 
  * It is *optional*. 
  * If you choose to specify it, it is suggested that you name this output file with `.csv` extension. 
  * If you do not specify the output file, a file with name `Mipepid_results.csv` will be automatically created under the same `./MiPepid` directory as the output file. 

The output `.csv` file contains the following columns:
* `sORF_ID`: the ID of an sORF. 
* `sORF_seq`: the DNA sequence of the sORF (containing the stop codon with the start codon as `ATG`).
* `transcript_DNA_sequence_ID`: the ID of the original DNA sequence in your input fasta file where this sORF is extracted.
* `start_at`: the 1-based starting position of this sORF on the original DNA sequence. 
* `end_at`: the 1-based ending position of this sORF on the original DNA sequence. 
* `classification`: the predicted class of the sORF, either `coding` or `noncoding`. 
* `probability`: the probability that this sORF is in the assigned class.

### How to run a demo
There is a sample DNA sequence file `sample_seqs.fasta` under the directory `./demo_files/`. You can try to run `MiPepid` on this file:

```sh
cd MiPepid
python3 ./src/mipepid.py ./demo_files/sample_seqs.fasta ./demo_files/MiPepid_results_on_sample_seqs.csv
```

This will output a file `MiPepid_results_on_sample_seqs.csv` under the same directory (`./demo_files/`).

### Regarding the datasets
`datasets.tar.gz` contains all major datasets used in the paper. 

<br>
<br>

[*]: Mengmeng Zhu, Michael Gribskov. MiPepid: Micropeptide identification tool using machine learning. *BMC Bioinformatics* 20, 559 (2019) doi:10.1186/s12859-019-3033-9
