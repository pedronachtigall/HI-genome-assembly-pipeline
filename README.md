<!---
![genome_assembly_pipeline](/genome_assembly_workflow.png)
--->

# Genome-assembly-pipeline
Bioinformatics pipeline and tutorial for performing genome assembly using HiFi and HiC data.

:construction:	**Under construction!** :construction:	

## Dependencies
 - [Python](https://www.python.org/) and [biopython](https://biopython.org/)
 - [cutadapt](https://github.com/marcelm/cutadapt)
 - [Trim_Galore](https://github.com/FelixKrueger/TrimGalore)
 - [MITGARD](https://github.com/pedronachtigall/MITGARD)
 - [kraken2](https://github.com/DerrickWood/kraken2)
 - [hifiasm](https://github.com/chhylp123/hifiasm)
 - [YaHS](https://github.com/c-zhou/yahs)
 - [chromap](https://github.com/haowenz/chromap)
 - [tidk](https://github.com/tolkit/telomeric-identifier)
 - [minimap2](https://github.com/lh3/minimap2)
 - [samtools](https://github.com/samtools/samtools)
 - bigWigToBedGraph
 - [deepTools](https://github.com/deeptools/deepTools)
 - [SeqKit](https://github.com/shenwei356/seqkit)
 - [PretextMap](https://github.com/sanger-tol/PretextMap)
 - [PretextView](https://github.com/sanger-tol/PretextView)

Ensure that all dependencies are installed and working properly.

## Summary
 - [Model species](#model-species)
 - [Trim adapters](#trim-adapters)
 - [Remove of contaminants](#remove-of-contaminants)
 - [Mitochondrial genome assembly](#mitochondrial-genome-assembly)
 - [Draft genome assembly](#draft-genome-assembly)
 - [References](#references)

## Model species
We will use the Golden lancehead (*Bothrops insularis*) as a model for this tutorial.

The genomic data is linked to the manuscript ["Unveiling the toxin repertoire of the golden lancehead: insights into the genomic evolution of a critically endangered species"](https://doi.org/10.1093/molbev/msaf058) published in *in prep* and it is available in NCBI under the under the project number [PRJNA679826](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA679826/).

The raw data is listed below:
| Sample ID | Data type | NCBI accession |
| :-------- | :-------: | :------------: | 
| SB1851    | HiFi | SRR32358152, SRR32358153 |
| SB1851    | HiC | SRR32358142, SRR32358143 |

## Trim adapters
We used [cutadapt](https://github.com/marcelm/cutadapt) to remove reads with adapters from HiFi data and [Trim_Galore](https://github.com/FelixKrueger/TrimGalore) to trim adapters and low-quality reads from HiC data.

```
#trim and merge hifi reads
cutadapt -j 20 -n 3 -O 35 --revcomp --discard-trimmed --anywhere="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" --anywhere="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" SRR32358152.hifi_reads.fastq.gz --output SRR32358152_HiFi_cutadapt.fastq
cutadapt -j 20 -n 3 -O 35 --revcomp --discard-trimmed --anywhere="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" --anywhere="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" SRR32358153.hifi_reads.fastq.gz --output SRR32358153_HiFi_cutadapt.fastq
cat SRR32358152_HiFi_cutadapt.fastq SRR32358153_HiFi_cutadapt.fastq > Binsu.hifi.trimmed.fastq

#trim hic reads
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR32358142_hic_tg SRR32358142_R1.fastq.gz SRR32358142_R2.fastq.gz
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR32358143_hic_tg SRR32358143_R1.fastq.gz SRR32358143_R2.fastq.gz
cat SRR32358142_hic_tg/SRR32358142_R1_val_1.fq.gz SRR32358143_hic_tg/SRR32358143_R1_val_1.fq.gz > Binsu.hic.R1.fastq.gz
cat SRR32358142_hic_tg/SRR32358142_R2_val_2.fq.gz SRR32358143_hic_tg/SRR32358143_R2_val_2.fq.gz > Binsu.hic.R2.fastq.gz
```

## Remove of contaminants
We checked and removed any bacterial and/or human reads (if needed) in the hifi data using [kraken2](https://github.com/DerrickWood/kraken2). TaxonomyID of Sauria 32561 was used to retrieve reads not matching to bacteria or human. We used a custom kraken2 database comprising the standard libraries and other squamata genomes. It is important to design a custom and reliable kraken2 database to ensure you have high quality reads for further steps. You can refer to the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual) to learn how to design a custom database.
```
mkdir kraken2 && cd kraken2
kraken2 --threads 20 --db krakendb ../Binsu.hifi.trimmed.fastq --report krakendb_report.txt --output krakendb_output.txt
extract_kraken_reads.py -k krakendb_output.txt --report krakendb_report.txt -s ../Binsu.hifi.trimmed.fastq -t 32561 -o ../Binsu.hifi.fastq --include-children
```

## Mitochondrial genome assembly
We used [MITGARD](https://github.com/pedronachtigall/MITGARD) to perform the mitogenome assembly. The available mitochondrial genome of *B. jararaca* (NC_030760.1) was used as the reference.
```
mkdir MITGARD_output && cd MITGARD_output
MITGARD-LR.py -s Binsu_mitogenome -m pacbio_hifi -r ../Binsu.hifi.fastq -R NC_030760.1.fasta
```

## Draft genome assembly
We used [hifiasm](https://github.com/chhylp123/hifiasm) to perform the draft genome assembly using HiFi and HiC reads to acquire the primary and both resolved haplotypes.
```
hifiasm -o Binsu -t32 --h1 Binsu.hic.R1.fastq.gz --h2 Binsu.hic.R2.fastq.gz Binsu.hifi.fastq
awk '/^S/{print ">"$2;print $3}' Binsu.hic.p_ctg.gfa > Binsu.hic.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' Binsu.hic.hap1.p_ctg.gfa > Binsu.hic.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' Binsu.hic.hap2.p_ctg.gfa > Binsu.hic.hap2.p_ctg.fasta
```

## Scaffold draft genome
We used [YaHS](https://github.com/c-zhou/yahs) to scaffold the primary genome assembled by hifiasm. To map reads against the draft genome, we used [chromap](https://github.com/haowenz/chromap).
```
mkdir YAHS_primary && cd YAHS_primary
ln -s ../Binsu.hic.p_ctg.fasta .
ln -s ../Binsu.hic.R1.fastq.gz hic_R1.fastq.gz
ln -s ../Binsu.hic.R2.fastq.gz hic_R2.fastq.gz

chromap -i -r Binsu.hic.p_ctg.fasta -o Binsu.hic.p_ctg.fasta.index
samtools faidx Binsu.hic.p_ctg.fasta

chromap --preset hic -r Binsu.hic.p_ctg.fasta -x Binsu.hic.p_ctg.fasta.index --remove-pcr-duplicates -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz --SAM -o aligned.sam -t 32
samtools view -@ 32 -bh aligned.sam | samtools sort -@ 32 -n > aligned.bam
rm aligned.sam

yahs Binsu.hic.p_ctg.fasta aligned.bam
```

## Review of scaffolded genome
We used [PretextMap](https://github.com/sanger-tol/PretextMap) and the [PretextView](https://github.com/sanger-tol/PretextView) to manually review the scaffolded genome following the [Rapid curation guide](https://gitlab.com/wtsi-grit/rapid-curation/-/tree/main). In addition, we used the read coveerage of hifi reads and telomeric repeats to help in the assembly review.

<!---

```
#map hic to generate the contact map
chromap -i -r yahs.out_scaffolds_final.fa -o yahs.out_scaffolds_final.fa.index
samtools faidx yahs.out_scaffolds_final.fa
chromap --preset hic -r yahs.out_scaffolds_final.fa -x yahs.out_scaffolds_final.fa.index --remove-pcr-duplicates -1 hic_R1.fastq -2 hic_R2.fastq --SAM -o aligned.sam -t 32
samtools view -@ 32 -bh aligned.sam | samtools sort -@ 32 -n > aligned.bam
rm aligned.sam

samtools view -h aligned.bam | PretextMap -o hic_map.pretext --sortby length --sortorder descend --mapq 0

#track - genomeCoverage - hifi
minimap2 -ax map-hifi -t 16 yahs.out_scaffolds_final.fa ../Binsu.hifi.fastq | samtools sort -@16 -O BAM -o hifi_unfilt.bam
samtools view -b -F 256 hifi_unfilt.bam > hifi.bam
rm hifi_unfilt.bam
samtools index hifi.bam
bamCoverage -b hifi.bam -o hifi.bw
bigWigToBedGraph hifi.bw  /dev/stdout | PretextGraph -i hic_map.pretext -n "hifi_cov" -o hifi_cov.pretext

#track - telomeric repeats
tidk search -s TTAGG --dir tidk_out --output TTAGG --fasta yahs.out_scaffolds_final.fa --extension bedgraph
PretextGraph -i hic_map.pretext -n "telomer" -o telomer.pretext < TTAGG_telomeric_repeat_windows.bedgraph
```

## References
If you use this tutorial or any of the resources/scripts, please consider citing: [Nachtigall et al., in prep](https://doi.org/10.1093/molbev/msaf058).

Please, cite the original manuscript of each tool used in this tutorial.
--->
