# Enzimatic engineering of DPA synthase from DPAHelix project

This repository corresponds to one of the five computational biology analyses of the DPAHelix project in the GOGEC Competition of 2026.

- **DPAHelix Members involved (No specific order)**: Andrea Catalina Lozano Garza, Miguel Angel Cerino Espinosa, Naomi Lizeth Renovato Castorena.
- **Advisor**: LBG Mauricio Reyes-Elizondo
- **PI**: Dr J Claudio Moreno-Rocha

> **Initial objective (Later modified)**: Improve the catalytic activity of DPA synthase through multiple alignment, molecular docking and directed mutagenesis.

The analysis was mostly conducted on databases or web platforms. Nonetheless, a code version of this workflow was created for reproducibility, run on Google Colab and Ubuntu through WSL.

## Workflow

1. Retrieve similar sequences using BLAST (100 per subunit)
2. Multiple Sequence Alignment MSA for A and B
3. Identify Conserved Residues & Functional Domains using Clustal-O & InterProScan
4. Activity 4: Predict A-B Complex Structure using Structure Predictive Algorithms
5. Prepare HTPA/Substrate Ligand Structure
6. Targeted Protein-Ligand Binding Affinity using Docking Algorithms
7. Analyze Docking Poses & Conserved Residues in Binding Site
8. Final Mutagenesis Target List (Conserved, Proximal, Catalytically Relevant)

## Activity 1: Identify available functional and annotated metadata

We accessed to enzymatic-related databases to retrieve all available information regarding DPA synthase (IUBMB, Rhea, BRENDA, UniProt, PDB). We first identified general metadata of the enzyme for later search:

- **EC Nomenclature**: 1.3.1.- (through UniProt database)
- **RHEA Number**: RHEA:47092 (through RHEA database)
- **Alternative names**: _spovF_, _dpaA_, _dpaB_, _spoVFA_, _spoVFB_, dipicolinate synthase, dipicolinic acid synthetase (through UniProt database)

_DPA Synthase is an oxidoreductase that catalyzes a dehydrogenation reaction on a C-CH group of its substrate, using NAD⁺ or NADP⁺ as an electron acceptor, to form dipicolinic acid (DPA)._

It is worth noting the dash as a fourth number in the EC nomenclature. There might be some uncertainty about its precise substrate specificity or mechanistic details for a final classification.

Considerations:

- Since there isn't a specific EC number for DPA synthase, we were unable to identify any quantitative, catalytic information in BRENDA database.
- Using the Rhea number of the reaction, we identified 155 DPA synthase-related entries in UniProt. These sequences naturally separated into two distinct groups corresponding to subunit-A and subunit-B. Notably, only one sequence corresponded to subunit A, with the remaining 154 corresponding to subunit B, suggesting a significant data imbalance.
- There are only two PDB structures available for each specific subunit (e.g., two for A and two for B), and none are bound to DPA or its precursor.

## Activity 2: Retrieve similar sequences using NCBI-BLAST

Given the lack of available data regarding quantitative catalytic information or a list of active site amino acids of the DPA synthase in the selected databases, we modified our objective **to identify potential catalytically relevant amino acids using comparative bioinformatics.** The identified highly conserved residues will serve as the candidate list for targeted mutagenesis experiments to improve catalytic activity.

We accessed to NCBI and identified template sequences of subunit A (GenBank: SPY12701.1; ```msa_results/SPY12701_1_result.xml```) and subunit B (GenBank: QRQ54142.1; ```msa_results/QRQ54142_1_result.xml```) from _Bacillus subtilis_. Each template was queried in Protein-BLAST to identify similar sequences. The results provided a table of 100 sequences per each subunit. For biosafety explanation, we strictly used the 100 sequences only for referencing potential conserved regions in the original sequence (MSA and InterPro Analysis).

Results from Subunit A:

- BLAST XML File: ```msa_results/PAXKDZ9J016_results.xml```
- BLAST DataFrame File: ```msa_results/PAXKDZ9J016_results.csv```
- BLAST FASTA File: ```msa_results/PAXKDZ9J016_blast.fasta```

Results from Subunit B:

- BLAST XML File: ```msa_results/PAMZAJN8016_results.xml```
- BLAST DataFrame File: ```msa_results/PAMZAJN8016_results.csv```
- BLAST FASTA File: ```msa_results/PAMZAJN8016_blast.fasta```

## Activity 3: Multiple Sequence Alignment (MSA) and InterProScan

A MSA through Clustal-Omega was performed separately for each subunit by aligning the _B. subtilis_ template sequence with the respective 100 homologous sequences retrieved via BLAST. This analysis aims to identify potential conserved regions and functionally critical residues. We installed Clustal-Omega using the following commands:

```shell
sudo apt-get update
sudo apt-get install clustalo
```

Results from Subunit A:

- ```msa_results/PAXKDZ9J016_msa_10_sequences.png``` (n sequences = 10, 25, 50, 75, 101)

Results from Subunit A:

- ```msa_results/PAMZAJN8016_msa_10_sequences.png``` (n sequences = 10, 25, 50, 75, 101)

Additionally, we ran the InterProScan tool to identify conserved functional domains within DPA synthase. For this analysis, we queried our _B. subtilis_ subunit-A and subunit-B sequences separately, using the top 5 most homologous sequences retrieved from the previous BLAST searches as supplementary context for domain identification.

Results from Subunit A:

- ```msa_results/SPY12701_1_uniprot_principal.png``` (Our query sequence)

For strictly reference-only:

- ```msa_results/WP_333517570_uniprot_second.png``` (Bacillus sp. 0102A)
- ```msa_results/WP_277780810_uniprot_third.png```(Bacillus atrophaeus)
- fourth file was not exported properly (UniProt Status: Failure)
- ```msa_results/WP_243940340_uniprot_fifth.png``` (MULTISPECIES Bacillus)

Results from Subunit B:

- ```msa_results/QRQ54142_1_uniprot_principal.png``` (Our query sequence)

For strictly reference-only:

- ```msa_results/WP_187002230_uniprot_second.png``` (Bacillus subtilis)
- ```msa_results/WP_306012148_uniprot_third.png``` (Bacillus pumilus)
- ```msa_results/WP_309331194_uniprot_fourth.png``` (Bacillus atrophaeus)
- ```msa_results/WP_254003724_uniprot_fifth.png``` (Bacillus velezensis)

These procedures were replicated on ```subunit_a_workflow.ipynb``` and ```subunit_b_workflow.ipynb```.

## Activity 4: Predict A-B Complex Structure using Structure Predictive Algorithms

The selected sequences of subunit A (GenBank: SPY12701.1; ```msa_results/SPY12701_1_result.xml```) and subunit B (GenBank: QRQ54142.1; ```msa_results/QRQ54142_1_result.xml```) were submitted to predict the three-dimensional complex structure using several algorithms hosted on the Tamarind Bio Platform: AlphaFold, ESMFold, Chai‑1, and Boltz‑2. The resulting complex structure was intended to be used for predicting a potential ligand-binding pocket. However, none of the algorithms successfully identified a stable complex structure. The output files are listed below:

AlphaFold:
- ```tamarind_results/complex_alphafold.zip``` (Overall files)
- ```tamarind_results/complex_alphafold_metrics.csv``` (Metrics results from algorithm)
- ```tamarind_results/complex_alphafold_model.pdb``` (PDB file)

ESMFold:
- ```tamarind_results/complex_esmmfold.zip```
- ```tamarind_results/complex_esmmfold_metrics.csv```
- ```tamarind_results/complex_esmmfold_model.pdb```

Chai-1:
- ```tamarind_results/complex_chai_1.zip```
- ```tamarind_results/complex_chai_1_metrics.csv```
- ```tamarind_results/complex_chai_1.pdb```

Boltz-2:
- ```tamarind_results/complex_boltz_2.zip```
- ```tamarind_results/complex_boltz_2_metrics.csv```
- ```tamarind_results/complex_boltz_2_model.pdb```

Subsequently, we identified a relevant article by Yang et al. (2025), which analyzed the catalytic activity of subunits A and B individually and in complex. The authors reported that both subunits exhibit catalytic activity, though subunit A appears to have greater activity and plays a more significant role in DPA production in *B. subtilis*. We continued our analysis predicting the binding pockets of both subunit A and B.

## Activity 5: Targeted Protein-Ligand Binding Affinity using Docking Algorithms

We predicted the tridimensional structures of subunit A and subunit B using AlphaFold.

Subunit A:
- ```tamarind_results/subunit_a_alphafold.zip```
- ```tamarind_results/subunit_a_alphafold_metrics.csv```
- ```tamarind_results/subunit_a_alphafold_model.pdb```

Subunit B:
- ```tamarind_results/subunit_b_alphafold.zip```
- ```tamarind_results/subunit_b_alphafold_metrics.csv```
- ```tamarind_results/subunit_b_alphafold_model.pdb```

To identify a potential binding pocket, the SMILES structure of the precursor of DPA, (2S,4S)-4-hydroxy-2,3,4,5-tetrahydrodipicolinate (HTPA) was retrieved. The predicted structures were separately added to different algorithms along with the SMILES of ligand: Boltz-2, DiffDock, and Autodock Vina. The output files are listed below:

Here are the corrected file names for each docking algorithm:

**Autodock Vina**

(Subunit A):
- `tamarind_results/subunit_a_autodock_vina.zip` (Overall files)
- `tamarind_results/subunit_a_autodock_vina_log.txt` (Log report)
- `tamarind_results/subunit_a_complex_autodock_vina.pdb` (Complex PDB file)
- `subunit_a_ligand_autodock_vina.sdf` (Ligand file)

(Subunit B):
- `tamarind_results/subunit_b_autodock_vina.zip`
- `tamarind_results/subunit_b_autodock_vina_log.txt`
- `tamarind_results/subunit_b_complex_autodock_vina.pdb`
- `subunit_b_ligand_autodock_vina.sdf`

**DiffDock**

(Subunit A):
- `tamarind_results/subunit_a_diffdock.zip`
- `tamarind_results/subunit_a_diffdock_log.log`
- `tamarind_results/subunit_a_complex_diffdock.pdb`
- `subunit_a_ligand_diffdock.sdf`

(Subunit B):
- `tamarind_results/subunit_b_diffdock.zip`
- `tamarind_results/subunit_b_diffdock_log.txt`
- `tamarind_results/subunit_b_complex_diffdock.pdb`
- `subunit_b_ligand_diffdock.sdf`

**Boltz-2** (Ligand files were retrieved from complex PDB file)

(Subunit A):
- `tamarind_results/subunit_a_boltz.zip`
- `tamarind_results/subunit_a_boltz_log.log`
- `tamarind_results/subunit_a_complex_boltz.pdb`
- `subunit_a_ligand_boltz.sdf`

(Subunit B):
- `tamarind_results/subunit_b_boltz.zip`
- `tamarind_results/subunit_b_boltz_log.txt`
- `tamarind_results/subunit_b_complex_boltz.pdb`
- `subunit_b_ligand_boltz.sdf`

The binding poses from all docking algorithms (Autodock Vina, DiffDock, and Boltz) for both subunits were analyzed using the ProLIF (Protein-Ligand Interaction Fingerprints) Python Library to retrieve detailed molecular interaction data.The output results of interactions were stored in `interaction_results/`.

## END