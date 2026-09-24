# fuzzy_k-mer_counter

Counts mers with possible "N" characters from fasta files.

## Usage

Counts number of occurrences of ${ACGTN}^k$ -mers in fasta file(s).
Reads fasta file paths and mers from parameters, batch files
or standard input (if no paths are given) for piping.

For batch operation, maximum $k$ is 32, unless run with the `-L`
option for max $k = 64$. For non-batch operation max $k$ is 64.

Usage: `fkm_count [fasta path] [-b file_path] [-r] [-s] [-L] [-t N] [-p path] [mers ...]`

| Parameter | Description |
|-----------|-------------|
| `fasta_file` | Path to fasta file (possibly compressed). Needs to be placed before any mer parameters. |
| `mers ...` | List of mers to search for. |
| `-b file_path` | Reads file containing one fasta_file and at least one mer per line. Overrides `fasta_path` and `mers ...`. |
| `-r` | Do not include reverse complements for non-palindromic mers. |
| `-s` | Report at most one match per fasta record. |
| `-t N` | Set the maximum number of threads. Only runs multithreaded on batch jobs. |
| `-p path` | Prefix path to add to (each) fasta_file. |
| `-L` | Extend max $k$ to 64, by sacrificing some performance. |

## Example

```bash
$ fkm_count -p ex_seq/ -b test.txt
```

With `test.txt`:

```
BARHL2_TCCAGT40NGAC_AI_1.fastq AATTGNNNNNNNNNNNNNTAAACG TAATTGNNNNNNNNNNTGTAA AATTACNNNNNNNNNNNNNNNTTAAA AAAGCNNNNNNNNNNCGTTTA ...
BARHL2_TCCAGT40NGAC_AI_3.fastq AATTGNNNNNNNNNNNNNTAAACG TAATTGNNNNNNNNNNTGTAA AATTACNNNNNNNNNNNNNNNTTAAA AAAGCNNNNNNNNNNCGTTTA ...
```

Will count the given mers from `ex_seq/BARHL2_TCCAGT40NGAC_AI_1.fastq` and `ex_seq/BARHL2_TCCAGT40NGAC_AI_3.fastq` in parallel.

## Building

Requires a modern C++ compiler (tested with GCC) and `zlib`.

### Conda

```bash
git clone https://github.com/saskeli/fuzzy_k-mer_counter.git
cd fuzzy_k-mer_counter
conda create -n fkm python -y # or use an existing environment
conda activate fkm
conda install -c conda-forge cmake make zlib llvm-openmp -y
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/fkm_count -h
```

### GNU Make

```bash
make count
./count -h
```