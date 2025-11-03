# Nextflow-GWAS-PopGen

## หัวข้อ
1. [บทนำ](#1-บทนำ)
2. [การใช้งาน Nextflow-GWAS-PopGen](#2-การใช้งาน-Nextflow-GWAS-PopGen)
3. [การเตรียมเครื่องมือและข้อมูลสำหรับ Nextflow-GWAS-PopGen](#3-การเตรียมเครื่องมือและข้อมูลสำหรับ-Nextflow-GWAS-PopGen)
4. [รายละเอียดขั้นตอนใน Nextflow-GWAS-PopGen](#4-รายละเอียดขั้นตอนใน-Nextflow-GWAS-PopGen)
5. [Output](#5-Output)
   
---
## 1. บทนำ
 Nextflow-GWAS-PopGen เป็น bioinformatics pipline ที่พัฒนาขึ้นสำหรับการเคราะห์ GWAS และ PopGen โดยจะมีขั้นตอนดังต่อไปนี้ 
1. การทำ Variant Annotations
2. การสร้าง Database (BuildCustomDB)
3. การเปรียบเทียบข้อมูล Variant ที่ซ้ำกับข้อมูล Variant ที่มีอยู่ (Comapare_VCF)
4. การดึงข้อมูล Variant Annotations ที่ซ้ำกับข้อมูล Variant ที่มีอยู่ (Call_ANN)
5. การรวมไฟล์ (Combine_VCF)
6. การใช้ SnpSift (ANN_SnpSift)

## 2-การใช้งาน-Nextflow-GWAS-PopGen 
### For vcf file
#### GWAS
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --analyze GWAS \
     --output result
```
### PopGen
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze PopGen \
     --PopGenTools phylo,faststructure,ipcaps \
     --output result
```
### BOTH
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze BOTH \
     --PopGenTools phylo,admixture,ipcaps \
     --output result
```
## For hmp file
### GWAS
```bash
nextflow run -profile gb main.nf --hmpFile data/inputonebase.hmp.txt --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --analyze GWAS \
     --output result
```
### PopGen
```bash
nextflow run -profile gb main.nf --hmpFile data/inputonebase.hmp.txt --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze PopGen \
     --PopGenTools phylo,faststructure,ipcaps \
     --output result
```
### BOTH
```bash
nextflow run -profile gb main.nf --hmpFile data/inputonebase.hmp.txt --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze BOTH \
     --PopGenTools phylo,admixture,ipcaps \
     --output result
```
## 6. Output
### ภาพรวม Output
```bash
Annotations
└── ANN_snpEff
     ├── {samples}.ann.vcf.gz 
     ├── {samples}_summary.genes.txt       
     └── {samples}_summary.html
```

```bash
Annotations_custom
├── BuildCustomDB
│    ├──snpeff_build.log
└── ANN_snpEff
     ├── {samples}.ann.vcf.gz 
     ├── {samples}_summary.genes.txt       
     └── {samples}_summary.html
```

```bash
Annotations
├── ANN_SnpSift
│    ├──{samples}_SnpSift.vcf.gz
└── ANN_snpEff
     ├── {samples}.ann.vcf.gz 
     ├── {samples}_summary.genes.txt       
     └── {samples}_summary.html
```

```bash
Annotations
├── Call_ANN
│    ├── {samples}_overlap_shared.vcf.gz
├── Combine_VCF
│    ├── {samples}_combine.vcf.gz
├── Compare_results
│    ├── {samples}_overlap.vcf.gz
│    └── {samples}_unique.vcf.gz
└── ANN_snpEff
     ├── {samples}.ann.vcf.gz 
     ├── {samples}_summary.genes.txt       
     └── {samples}_summary.html
```
