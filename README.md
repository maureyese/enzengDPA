# Enzimatic engineering of DPA synthase from DPAHelix project

This repository corresponds to one of the five computational biology analyses of the DPAHelix project in the GOGEC Competition of 2026.

> **Objective**: Improve the catalytic activity of DPA synthase through multiple alignment, molecular docking and directed mutagenesis.

We worked using VS Code, along with the WSL extension (Ubuntu 22.04) as our development environmnet (Only-Windows users).

1. WSL Installation: https://learn.microsoft.com/en-us/windows/wsl/install
2. VS Code Installation: https://code.visualstudio.com/docs/setup/windows#_install-vs-code-on-windows
3. WSL extension on VS Code: https://code.visualstudio.com/docs/remote/wsl

## Activity 1: Identify available functional and annotated metadata

We accessed to enzymatic-related databases to retrieve all available information regarding DPA synthase. We first identified general metadata of the enzyme for later search:

- **EC Nomenclature**: 1.3.1.- (through literature)
- **RHEA Number**: RHEA:47092 (through RHEA database)
- **Alternative names**: _spovF_, _dpaA_, _dpaB_, _spoVFA_, _spoVFB_, dipicolinate synthase (through literature and UniProt database)

_DPA Synthase is an oxidoreductase that catalyzes a dehydrogenation reaction on a C-CH group of its substrate, using NAD⁺ or NADP⁺ as an electron acceptor, to form dipicolinic acid (DPA)._

It is worth noting the dash as a fourth number in the EC nomenclature. There might be some uncertainty about its precise substrate specificity or mechanistic details for a final classification.

Considerations:

- Since there isn't a specific EC number for DPA synthase, we were unable to identify any quantitative, catalytic information in BRENDA database.
- Using the Rhea number of the reaction, we identified 155 DPA synthase-related data. Nonetheless, most of all annotated data is subunit-specific, either A or B. Two distinct groups (subunit-A, and subunit-B) were formed among the 155 sequences.
- There are no DPA synthase structures bounded with DPA. There are only 2 PDB structure per each specific subunit.

## Activity 2: Multiple Sequence Alignment (MSA)

_(under development...)_