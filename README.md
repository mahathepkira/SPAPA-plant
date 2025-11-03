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
1. Selected data type
2. Preprocess Data
3. GWAS Analysis
4. Selected PruneLD Tools and Population anlysis
5. Visualize QC
![ภาพ nextflow](GWASPopGen.drawio.png)

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
#### PopGen
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze PopGen \
     --PopGenTools phylo,faststructure,ipcaps \
     --output result
```
#### BOTH
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze BOTH \
     --PopGenTools phylo,admixture,ipcaps \
     --output result
```
### For hmp file
#### GWAS
```bash
nextflow run -profile gb main.nf --hmpFile data/inputonebase.hmp.txt --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --analyze GWAS \
     --output result
```
#### PopGen
```bash
nextflow run -profile gb main.nf --hmpFile data/inputonebase.hmp.txt --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze PopGen \
     --PopGenTools phylo,faststructure,ipcaps \
     --output result
```
#### BOTH
```bash
nextflow run -profile gb main.nf --hmpFile data/inputonebase.hmp.txt --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze BOTH \
     --PopGenTools phylo,admixture,ipcaps \
     --output result
```
## 3. การเตรียมเครื่องมือและข้อมูลสำหรับ Nextflow-GWAS-PopGen
### เครื่องมือ
1. Nextflow: version 19 
### การเตรียม config
ผู้ใช้งานสามารปรับแต่งเครื่องมือที่ใช้งานในไฟล์ gb.config ให้เหมาะสมกับทรัพยากรในเครื่อง โดย gb.config จะทำงานรวมกับ nextflow.config โดยจะใช้ตัวเลือก `-profile` เพื่อเลือก config ที่จะใช้งาน
```bash
```
## 4. รายละเอียดขั้นตอนใน Nextflow-GWAS-PopGen
### Selected data type
### Preprocess Data
### GWAS Analysis
### Selected PruneLD Tools and Population anlysis
### Visualize QC
## 5. Output
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
