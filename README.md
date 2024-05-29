# utriobin: trio binning using kmers

## Sep 1. Compute kmers using any kmer counter, in our case jellyfish

```bash
zcat ./samples/HG06807/parents/pan010/blood/*gz ./samples/HG06807/parents/pan010/cell1/*gz ./samples/HG06807/parents/pan010/cell2/*gz | jellyfish count -m 23 -s 5G -t 100 -o pan010.23.jf2 -

jellyfish count -m 23 -s 5G -C -t 100 -o pan011.23.jf2 <(zcat ./samples/HG06807/parents/pan011/blood/*gz ./samples/HG06807/parents/pan011/cell1/*gz ./samples/HG06807/parents/pan011/cell2/*gz)

jellyfish count -m 23 -s 5G -C -t 100 -o HG06807.23.jf2 HG06807.fasta
```

And dump them to kmer/tf files:

```bash
jellyfish dump -c -t -L 2 pan010.23.jf2 > pan010.23.kmers
jellyfish dump -c -t -L 2 pan011.23.jf2 > pan011.23.kmers
```

```bash

## Sep 2. From child genome take only kmers that have frequency == 1

```bash
jellyfish dump -c -t -U 1 HG06807.23.jf2 > HG06807.23.L1.kmers
```

## Step 3. From the parents, take only kmers that are in the child

```bash
./intersect HG06807.23.L1.kmers pan010.23.kmers pan010.23.L1.kmers
./intersect HG06807.23.L1.kmers pan011.23.kmers pan011.23.L1.kmers
```

## Step 4. Aggregate data to the one file

```bash
./aggregate HG06807.23.L1.kmers pan010.23.L1.kmers pan011.23.L1.kmers freq_vectors.tsv
```

Fields:

| Kmer | Freq1 (File1) | Freq2 (File1) | Freq1 (File2) | Freq2 (File2) |
|------|---------------|---------------|---------------|---------------|

## Step 5. Classify contigs from the fasta file

```python
from trseeker.seqio.fasta_file import sc_iter_fasta_brute
from trseeker.tools.sequence_tools import get_revcomp
from tqdm import tqdm

child_kmers_file = "./HG06807.23.L1.kmers"
parents_data_file = "./freq_vectors.tsv"

child_kmers = {}
with open(child_kmers_file) as fh:
    for line in tqdm(fh):
        kmer, tf = line.strip().split()
        child_kmers[kmer] = int(tf)
        

parents_kmers = {}
with open(parents_data_file) as fh:
    fh.readline()
    for line in tqdm(fh):
        kmer, tf1, tf2, tf3, tf4 = line.strip().split()
        parents_kmers[kmer] = (int(tf2), int(tf4))
        
genome = {}
for header, seq in sc_iter_fasta_brute(child_genome_fasta_file):
    genome[header.strip()[1:]] = seq.upper()

classifications = {}

k = 23
for header in genome:
    seq = genome[header]
    print(header, len(seq))
    p = 0
    m = 0
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        if not kmer in child_kmers:
            kmer = get_revcomp(kmer)
        if child_kmers[kmer] == 1 and parents_kmers[kmer].count(0) == 1:
            if parents_kmers[kmer][0] == 0:
                p += 1
            else:
                m += 1
    classifications[header] = (p, m)
```