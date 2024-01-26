# April
[![Open Source Starter Files](https://github.com/nrminor/April/actions/workflows/open-source-starter.yml/badge.svg)](https://github.com/nrminor/April/actions/workflows/open-source-starter.yml) [![Rust CI](https://github.com/nrminor/April/actions/workflows/rust-ci.yml/badge.svg)](https://github.com/nrminor/April/actions/workflows/rust-ci.yml) 

### Backstory
First off, a disclaimer: April is not built yet and is in its very early stages. That said, I can still provide a little backstory.

For the past year, I've been working on a tool for finding noteworthy SARS-CoV-2 lineages called [ALPINE](https://github.com/nrminor/ALPINE). One part of the pipeline makes a pairwise nucleotide distance matrix of sequences from each month of the Covid-19 pandemic, with the aim being to find sequences that are highly mutated. However, we kept getting the weird result that hardly any sequences were springing up as high-distance, and those that were springing up were hardly impressive. Among other possible explanations, I decided to start exploring whether a) sequences contaminated with primer sequences, and b) sequences bioinformatically contaminated with the SARS-CoV-2 reference sequence, were artificially inflating the pairwise distances for most sequences each month, thereby making it harder for noteworthy sequences to stand out.

So, the idea of April is to allow noteworthy samples to "spring up" from a background of potentially contaminated consensus sequences--because in most parts of the Northern hemisphere, most plants start springing up in April.

April could also be an acronym for the following (I just haven't chosen one yet):
- **A**pproximate **P**ermutations of **R**eference **I**nfiltration e**L**iminated
- **A**ccuracy **P**rofiler and **R**eference **I**nfluence **L**imiter
- **A**ccurate **P**urification of **R**eference **I**ntegration **L**ayers

April is being built in Rust with an emphasis on a rich and configurable command line interface, parallel processing and filtering of FASTA records, fast kmer containment computations, and informative report generation.
