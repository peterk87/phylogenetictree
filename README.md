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

3. Download the pipeline and test it on a minimal dataset with a single command: the enviroment contains R package so running with Conda profile will be very slow, it is recommended to run with Docker profile in which mamba is used to install packages

    ```bash
    nextflow run nhhaidee/phylogenetictree -with-docker nhhaidee/phylogenetic:dev1.0 --input '/path/to/consensus_sequences/*.fasta' --reference_name ='MN908947.3' --reference_fasta = '/path/to/ref_seq/nCoV-2019.reference.fasta'
    ```

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nhhaidee/phylogenetictree -with-docker nhhaidee/phylogenetic:dev1.0 --filter_gisaid true
    ```

See [usage docs](https://nf-co.re/phylogenetictree/usage) for all of the available options when running the pipeline.


## Documentation

The nf-core/phylogenetictree pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/phylogenetictree/usage) and [output](https://nf-co.re/phylogenetictree/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

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
