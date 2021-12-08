title: 'Good news everyone! Bender is a set of Genomic tools for deployment on Slurm managed systems'
tags:

- Go
- genomics
- command line
  authors:
- name: Salvador Daniel Rivas-Carrillo ^[first author, corresponding, author]
  orcid: 0000-0002-0143-2687
  affiliation: 1
- name: Deparment of Medical Biochemistry and Microbiology, Uppsala University
  index: 1
  date: 1 December
  bibliography: paper.bib

---

# Summary

The everincreasing complexity of analysis requiring computational resources is overhelming. Computational tools and analysis pipelines are being made public at an incresing pace in parallel with large amounts of data become available. And with these, the hardous task of managing resources, often on remote clusters. To meet these chanllenges, I developed Bender, a command line application to integrate and control genomic pipelines locally or on SLURM-controlled systems.

# Description

Bender is written in Go, a statically typed programming language that compiles to an cross-platform executable, which makes deployment on different hosts natural. Bender is built on the libraries Cobra and Viper. Thus, it uses a hierchical organization to access commands and subcommands, each with different arguments, functionality and extensive documentation that can be access directly. Bender accepts default arguments, parameters passed on a configuration file or parameters directly on the command line, with this order of priority.

# Statement of need

# Software Repository

The software is available as source code from https://github.com/DanielRivasMD/Bender/ under the General Public License GLP-3.0.

# Acknowledgements

<!-- TODO: acknowledge not author contributors -->

# References
