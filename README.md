# ![nf-core/phylogenetictree](docs/images/nf-core-phylogenetictree_logo.png)

**Phylogenetic Tree**.

[![GitHub Actions CI Status](https://github.com/nf-core/phylogenetictree/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/phylogenetictree/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/phylogenetictree/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/phylogenetictree/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/phylogenetictree.svg)](https://hub.docker.com/r/nfcore/phylogenetictree)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23phylogenetictree-4A154B?logo=slack)](https://nfcore.slack.com/channels/phylogenetictree)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/phylogenetictree** is a bioinformatics pipeline for sars-cov2 phylogenetic analysis, given a consensus sequences the workflow will output phylogenetic tree and SNP information. The pipeline also allows to filter and find the most related sequences in GISAID

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Install nhhaidee/phylogenetic

Nextflow will automatically download the latest version of pipeline. You can show the pipeline help message with usage information with:

```bash
nextflow run nhhaidee/phylogenetictree --help
```

4. Start running your own analysis!

<!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

```bash
nextflow run nhhaidee/phylogenetictree -with-docker nhhaidee/phylogenetic:dev1.0 --filter_gisaid false --reference_name  'MN908947.3' --reference_fasta '/path/to/nCoV-2019.reference.fasta' --input '/path/to/consensus/*.fasta'
OR
nextflow run nhhaidee/phylogenetictree -with-docker nhhaidee/phylogenetic:dev1.0 --filter_gisaid true --gisiad_sequences /path/to/seq.fasta --gisiad_metadata /path/to/metadata.tsv --sample_lineage B.1.1.306 --region 'North America' --country 'Canada
```

##  Usage

Options for nuilding phylogenetic tree of consensus sequences, typical command to run is as follow:
    nextflow run nhhaidee/phylogenetictree -with-docker nhhaidee/phylogenetic:dev1.0 --filter_gisaid false --reference_name  'MN908947.3' --reference_fasta '/path/to/nCoV-2019.reference.fasta' --input '/path/to/consensus/*.fasta'

    --input                      The direroty path to consensus sequences
    --reference_name             Name of reference sequenc (MN908947.3)
    --reference_fasta            Directory path to reference fasta file

Options for filtering sequences against GISIAD, typical command is as follow:
    nextflow run nhhaidee/phylogenetictree -with-docker nhhaidee/phylogenetic:dev1.0 --filter_gisaid true --gisiad_sequences /path/to/seq.fasta --gisiad_metadata /path/to/metadata.tsv --sample_lineage B.1.1.306 --region 'North America' --country 'Canada

    --filter_gisaid              Filter against GISIAD sequences or not (Default is false)
    --gisiad_sequences           Directory path to GISIADS sequences (Download GISIAD sequences form gisaid.org), this is mandotory of filter_gisaid is true
    --gisiad_metadata            Direcorty path to metadata file of GISIAD Sequences
    --country                    Find sequences belong to country (Canada)
    --region                     Find sequences belong to region (North America)
    --sample_lineage             Lineage of sample that we want to filters
    --lmin                       Remove sequences that lenght < lmin
    --lmax                       Remove sequences that lenght > lmax
    --xambig                     Remove sequences that have number of ambiguous sequences > xambig


## Credits

nf-core/phylogenetictree was originally written by Hai Nguyen.

We thank the following people for their extensive assistance in the development
of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#phylogenetictree` channel](https://nfcore.slack.com/channels/phylogenetictree) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/phylogenetictree for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

In addition, references of tools and data used in this pipeline are as follows:

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
