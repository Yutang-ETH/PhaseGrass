# PhaseGrass
## What is PhaseGrass
PhaseGrass is a single-sample reference-based phasing workflow, compatible with both PacBio HiFi and ONT data, generating chromosome-level haplotype-resolved diploid assembly containing balanced haplomes for highly heterozygous genomes. 

## Motivation, Novelty and Highlight
A high level of heterozygosity might cause graph-based phasing methods to generate unbalanced haplomes. This has motivated us to develop PhaseGrass, a reference-based phasing workflow generating chromosome-level haplotype-resolved diploid assembly with balanced haplomes for highly heterozygous genomes. PhaseGrass was initially designed and tested for highly heterozygous forage grass species, such as _L. multiflorum_ and _L. perenne_, which possess 1 SNP per 20 bp and extensive structural variations between haplotypes through chromosomes. For other diploid species with a similar level of heterozygosity, PhaseGrass should also be applicable. Compared to other reference-based phasing methods (DipAsm, WhatsHap), the novelty of PhaseGrass is that we combine reference-based phasing with k-mer-based sequence binning instead of alignment-based read binning. With haplotype-specific k-mers, PhaseGrass can bin not only long reads but also unitigs to different haplotypes. More importantly, binning sequences with haplotype-specific k-mers can mitigate reference bias, especially when sequence divergence between haplotypes is high. PhaseGrass is compatible with both ONT and PacBio HiFi data and can use both technologies jointly with PacBio HiFi data to accomplish assembly and ONT data to complete phasing. This allows users to flexibly choose the sequencing technology and take full advantage of both technologies. 





Fig. 1 | PhaseGrass workflow overview
![phasegrass_workflow_new](https://github.com/Yutang-ETH/PhaseGrass/assets/84848653/c4c43412-66dc-4c9f-ab58-5b4330431206)






To generate a chromosome-level haplotype-resolved diploid assembly, PhaseGrass includes the following steps (Fig. 1): 

Step 1: generating a primary assembly with either ONT or PacBio HiFi reads. For PacBio HiFi reads, primary contigs are generated with Hifiasm. Unitigs generated by Hifiasm are kept for haplotype partition in step 4. For ONT reads, NextDenovo is chosen to generate the primary contigs, which are then polished with NextPolish.

Step 2: generating a chromosome-level unphased haploid assembly. PurgeGrass removes redundant allelic contigs from the primary assembly. Tritex scaffolds the deduplicated primary contigs to pseudo-chromosomes. With the pseudo-chromosomes, JuiceGrass generates a Hi-C contact map for manual curation.

Step 3: reference-based phasing with the chromosome-level unphased haploid assembly. With whole-genome sequencing (WGS) short reads, long reads and Hi-C data, DipGrass phases SNPs through pseudo-chromosomes resulting in one chromosome-level phase block per pseudo-chromosome.

Step 4: Partitioning long reads or unitigs to different haplotypes. SortGrass employs haplotype-specific k-mers extracted from the chromosome-level phase blocks to bin ONT reads or unitigs to different haplotypes, resulting in two sets of ONT reads or unitigs.

Step 5: assembling each set of ONT reads independently and haplotype-aware polishing with PolishGrass for the resulting two sets of contigs. This step is only intended for ONT reads.

Step 6: final scaffolding and manual curation. ScaffoldGrass scaffolds the two ONT-derived contig assemblies or the two PacBio HiFi-derived unitig assemblies independently, resulting in two chromosome-level assemblies. JuiceGrass concatenates these two chromosome-level assemblies to construct a diploid Hi-C contact map for manual curation, which returns the final chromosome-level haplotype-resolved diploid assembly containing two haplomes.
