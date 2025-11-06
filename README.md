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
### สำหรับ input vcf file
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
### สำหรับ input hmp file
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
### Options
- `--vcfgzFile` = ไฟล์ input vcf (จำเป็น:กรณี input เป็น vcf)
- `--hmpFile` = ไฟล์ input hmp (จำเป็น:กรณี input เป็น hmp)
- `--traitsFile` = ไฟล์ traits สำหรับการทำ GWAS และ PopGen (จำเป็น)
- `--QCTools` = เครื่องมือในการทำ QC ข้อมูล (PLINK หรือ BCFTools|ค่าเริ่มต้น:PLINK)
- `--homozygous` = การ filetr homozygous (homozygous,null|ค่าเริ่มต้น:null)
- `--maf` = ค่าในการกรอง Minor Allele Frequency (ค่าเริ่มต้น:0.05)
- `--geno` =  ค่าในการกรอง Genotypeing rate (missing rate) (ค่าเริ่มต้น:0.02)
- `--method` = model ในการทำ GWAS (CMLM, MLM, GLM, MLMM, FarmCPU, FASTLMM|ค่าเริ่มต้น:CMLM)   
- `--PruneLDTools` = เครื่องมือในการทำ Prune LD ข้อมูล (PLINK หรือ BCFTools|ค่าเริ่มต้น:PLINK)
- `--analyze` = เลือกการวิเคราะห์ (GWAS, PopGen, BOTH|ค่าเริ่มต้น:GWAS)
- `--PopGenTools` = เลือกการวิเคราะห์ PopGen (phylo, admixture, ipcaps, faststructure|ค่าเริ่มต้น:admixture)
- `--output` = โฟลเดอร์หรือไฟล์ output (จำเป็น)
## 3. การเตรียมเครื่องมือและข้อมูลสำหรับ Nextflow-GWAS-PopGen
### เครื่องมือ
1. Nextflow: version 19
2. Preprocess Data : 
3. GWAS Analysis :
4. Selected PruneLD Tools and Population anlysis :
5. Visualize QC :
### การเตรียม config
ผู้ใช้งานสามารปรับแต่งเครื่องมือที่ใช้งานในไฟล์ gb.config ให้เหมาะสมกับทรัพยากรในเครื่อง โดย gb.config จะทำงานรวมกับ nextflow.config โดยจะใช้ตัวเลือก `-profile` เพื่อเลือก config ที่จะใช้งาน
```bash
process {
  executor = 'slurm'
  queue = 'memory'
  cache = 'lenient'

  withName: TrimHmpAndTraits {
    module = 'Python/3.7.4-GCCcore-8.3.0'
    cpus = 8
  }

  withName: ConvertHapMapToMapPed {
    module = 'TASSEL/5.2.59'
    cpus = 12
    memory = "128 GB"
  }

  withName: CheckTypeHapmap {
    module = 'TASSEL/5.2.59'
    cpus = 2
    memory = "8 GB"
  }

  withName: ConvertBedBimFamToMapPed {
    module = 'PLINK/1.9b_4.1-x86_64'
    cpus = 2
    memory = "16 GB"
  }

  withName: ConvertBedBimFamToVCF {
    module = 'PLINK/1.9b_4.1-x86_64:BCFtools/1.17-GCC-12.2.0'
    cpus = 2
    memory = "16 GB"
  }

  withName: ConvertMapPedToVCF {
    module = 'PLINK/1.9b_4.1-x86_64:BCFtools/1.17-GCC-12.2.0'
    cpus = 2
    memory = "8 GB"
  }

  withName: ConvertVCFToHapMap {
    module = 'TASSEL/5.2.59'
    cpus = 2
    memory = "16 GB"
  }

  withName: ConvertHapMapToVCF {
    module = 'TASSEL/5.2.59:BCFtools/1.17-GCC-12.2.0'
    cpus = 2
    memory = "16 GB"
  }

  withName: ConvertVCFToMapPed {
    module = 'PLINK/1.9b_4.1-x86_64'
    cpus = 12
    memory = "128 GB"
  }

  withName: ConvertVCFToBedBimFam {
    module = 'PLINK/1.9b_4.1-x86_64'
    cpus = 12
    memory = "128 GB"
  } 

  withName: ConvertMapPedToHapMap {
    module = 'TASSEL/5.2.59'
    cpus = 12
    memory = "128 GB"
  }

  withName: VCFstats_GWAS {
    module = 'VCFtools/0.1.16-GCC-11.3.0:Python/3.10.4-GCCcore-11.3.0'
    cpus = 4
    memory = "8 GB"
  }

  withName: VCFstats_PopGen {
    module = 'VCFtools/0.1.16-GCC-11.3.0:Python/3.10.4-GCCcore-11.3.0'
    cpus = 4
    memory = "8 GB"
  }  

  withName: IPCAPs {
    container = '/nbt_main/share/singularity/r_ldheatmap_ipcaps.sif'
    cpus = 12
    memory = "128 GB"
  }

  withName: TestGAPIT {
    container = '/nbt_main/share/singularity/r_ldheatmap_ipcaps.sif'
  }

  withName: Admixture {
    module = 'ADMIXTURE/1.3.0'
    cpus = 12
    memory = "128 GB"
  }

  withName: AdmixtureIter1 {
    module = 'ADMIXTURE/1.3.0'
    cpus = 12
    memory = "128 GB"
  }

  withName: AdmixtureIter2 {
    module = 'ADMIXTURE/1.3.0'
    cpus = 12
    memory = "128 GB"
  }

  withName: AdmixtureIter3 {
    module = 'ADMIXTURE/1.3.0'
    cpus = 12
    memory = "128 GB"
  }

  withName: SummarizeAdmixture {
    cpus = 1
  }

  withName: Phylogenetic {
    //module = 'Python/2.7.15-foss-2018b:RAxML/8.2.9-foss-2018b-hybrid-avx2'
    module = 'PLINK/2.00a3.6-GCC-11.3.0:Python/2.7.15-foss-2018b'
    cpus = 12
    memory = "256 GB"
  }

  withName: StructureByFastStructure {
    container = '/nbt_main/share/singularity/faststructure:1.0--py27hfaf7806_1'
    //module = 'Python/2.7.15-foss-2018b:fastStructure/1.0-foss-2019a-Python-2.7.15'
    cpus = 12
    memory = "128 GB"
  }

  withName: QualityControlByPLINK {
    module = 'PLINK/1.9b_4.1-x86_64'
    cpus = 4
    memory = "8 GB"
  }

  withName: PruneLDByPLINK {
    module = 'PLINK/1.9b_4.1-x86_64'
    cpus = 2
    memory = "16 GB"
  }

  withName: QualityControlByBCFTools {
    module = 'BCFtools/1.17-GCC-12.2.0'
    cpus = 4
    memory = "8 GB"
  }

  withName: PruneLDByBCFTools {
    module = 'BCFtools/1.17-GCC-12.2.0'
    cpus = 2
    memory = "16 GB"
  }

  withName: GWASAnalysis {
    container = '/nbt_main/share/singularity/r_ldheatmap_ipcaps.sif'
    cpus = 12
    memory = "256 GB"

  }

  withName: WEKAConvertFromCSV {
    module = 'WEKA/custom-Java-1.8.0_201'
    cpus = 2
    memory = "16 GB"
  }

  withName: WekaFilterNumericAttributesToNominal {
    module = 'WEKA/custom-Java-1.8.0_201'
    cpus = 2
    memory = "16 GB"
  }

  withName: WekaInfoGainAttributeEval {
    module = 'WEKA/custom-Java-1.8.0_201'
    cpus = 2
    memory = "16 GB"
  }

  withName: ParseFilterAndBayesClassifyByTopN {
    module = 'WEKA/custom-Java-1.8.0_201'
    cpus = 1
    memory = "4 GB"
  }

}

singularity {
    enabled = true
    autoMounts = true
}
```
## 4. รายละเอียดขั้นตอนใน Nextflow-GWAS-PopGen
### Preprocess Data
```bash
process QualityControlByPLINK {

  tag { key }

  input:
  tuple key, file(map), file(ped)

  output:
  tuple key, file("${prefix}.clean.map"), file("${prefix}.clean.ped")

  script:
  prefix=map.baseName
  """
  plink --file $prefix \
    --no-fid \
    --no-parents \
    --no-sex \
    --no-pheno \
    --geno ${params.geno} \
    --recode \
    --allow-extra-chr \
    --out ${prefix}.clean \
    --maf ${params.maf} \
  """
}
```
```bash
process QualityControlByBCFTools {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${prefix}.clean.vcf.gz")

  script:
  prefix=vcfgz.baseName
  """
  bcftools view -i 'MAF >= ${params.geno} && F_MISSING <= ${params.maf}' \
    -Oz -o ${prefix}.clean.vcf.gz ${vcfgz}
  """
}
```
### GWAS Analysis
```bash
process GWASAnalysis {

	tag { "$key:$method" }

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

	input:
	tuple key, traits, file(traitTxt), file(hmp)
	val method

	output:
	tuple key, file("*.csv"), file("*.pdf")

	script:
	cmds="""
  cat $traitTxt | cut -f1,3- > traits.trimmed.txt
	"""
	if (method == "FASTLMM")
		cmds+"""
		GAPIT BICmodelSelection traits.trimmed.txt $hmp
		for trait in ${traits.join(" ")}; do
			optimalPCsVal=\$(optimalPCs GAPIT.MLM.\${trait}.BIC.Model.Selection.Results.csv)
			if (\$((optimalPCsVal != "0")));
			then
				cut -d "," -f 1-\$((optimalPCsVal+1)) GAPIT.PCA.csv > GAPIT.PCA.\${optimalPCsVal}.\${trait}.csv
				GAPIT $method traits.trimmed.txt $hmp GAPIT.PCA.\${optimalPCsVal}.\${trait}.csv
			else
				GAPIT $method traits.trimmed.txt $hmp
			fi
		done
		"""
	else if (method != "FASTLMM")
		cmds+"""
		GAPIT $method traits.trimmed.txt $hmp
		"""
	//else if (method == "cMLM")
	//	cmds+"""
	//	GAPIT $method traits.trimmed.txt $hmp
	//	"""
	//else if (method == "GLM")
	//	cmds+"""
	//	GAPIT $method traits.trimmed.txt $hmp
	//	"""
	else 
		throw new Exception("Unknown method...")
}
```
### Selected PruneLD Tools and Population anlysis
```bash
process PruneLDByPLINK {

  tag { key }

  input:
  tuple key, file(map), file(ped)

  output:
  tuple key, file("${prefix}.final.bed"), file("${prefix}.final.bim"), file("${prefix}.final.fam")

  script:
  prefix=map.baseName
  """
  plink --file ${prefix} --indep-pairwise 100 5 0.2 --allow-extra-chr --out ${prefix}.pruneLD
  plink --file ${prefix} --extract ${prefix}.pruneLD.prune.in --make-bed --allow-extra-chr --out ${prefix}.final
  """
}
```
```bash
process PruneLDByBCFTools {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${prefix}.final.vcf.gz")

  script:
  prefix=vcfgz.baseName
  """
  bcftools +prune ${vcfgz} -m 0.2 -w 100 -N maxAF -Oz -o ${prefix}.final.vcf.gz
  """
}
```
```bash
process IPCAPs {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  errorStrategy "ignore"

  tag { key }

  input:
  tuple key, file(bed), file(bim), file(fam), traits, file(traitTxt)

  output:
  tuple key,
    file("node*.txt"),
    file("cluster_output/images/*"),
    file("cluster_output/RData/*"),
    file("cluster_output/*.html"),
    file("cluster_output/groups.txt")

  script:
  """
  cat $traitTxt | awk 'NR > 1 { if (\$2 == "NA") print \$1" "\$1; else print \$1" "\$2; }' > sampleLabel.txt
  IPCAPS $bed sampleLabel.txt || true
  for rdataFile in \$(find . -name "node*.RData"); do
    extractIPCAPSnode $fam \$rdataFile
  done
  """
}

```
```bash
process Admixture {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  tag { "$key:$k" }

  input:
  tuple key, file(bed), file(bim), file(fam)
  each k

  output:
  tuple k, key, file(bed), file(bim), file(fam), file("*.Q"), file("${key}.${k}.cv.errors.csv")

  script:
  """
  admixture $bed $k --cv=10 -j${task.cpus}
  perl -n -e'/^CV error \\(K=(\\d+)\\): ([\\d\\.]*)\$/ && print \$1,",",\$2' .command.log > ${key}.${k}.cv.errors.csv
  """
}
```
```bash
process SummarizeAdmixture {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  tag { key }

  input:
  tuple optimalK, key, file(bed), file(bim), file(fam), file(minCVq), minCVerror
  file(minCVerrors)

  output:
  tuple optimalK, key, file("cv.errors.csv"), file("optimalK.csv")

  script:
  """
  echo "K,CrossValidationError" > cv.errors.csv
  for eachCVErrorByK in \$(find . -name "*.cv.errors.csv" | sort -k3,3g -t '.'); do
    cat \$eachCVErrorByK >> cv.errors.csv
    echo "" >> cv.errors.csv
  done
  echo "$optimalK,$minCVerror" >> optimalK.csv
  """
}

```
```bash
process GatherAdmixtureResults {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple file(cvQs), file(cvErrors)

  output:
  tuple file(cvQs), file(cvErrors)

  script:
  """
  echo "Just Gathering results to same s3 path..."
  """
}
```
```bash
process StructureByFastStructure {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  errorStrategy "ignore"

  tag { "$key:$k" }

  input:
  tuple key, file(bed), file(bim), file(fam)
  each k

  output:
  tuple k, key, file("*.meanQ")

  script:
  prefix=bed.baseName

  """
  structure.py -K $k --input=$prefix --output=$key

  """
}
```
```bash
process Phylogenetic {

  errorStrategy 'ignore'

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  tag { key }

  input:
  tuple key, file(vcf)

  output:
  tuple key, file("*.min4.phy.treefile")

  script:
  prefix = vcf.simpleName.replaceFirst(/\.vcf$/, "").replaceFirst(/\.gz$/, "")
  """
  vcf2phylip.py -i ${vcf} -o ${prefix}.min4.phy

  iqtree2 -s ${prefix}.min4.phy -B 1000 -T Auto

  """
}
```
### Visualize QC
```bash
process VCFstats {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(vcf)

  output:
  file("*.frq")
  file("*.lmiss")
  file("*.TsTv.summary")
  file("*.summary")
  file("*.csv")
  script:
  prefix=vcf.baseName

  """
  vcftools --gzvcf "$vcf" --freq --out "${prefix}"
  vcftools --gzvcf "$vcf" --missing-site --out "${prefix}"
  vcftools --gzvcf "$vcf" --TsTv-summary --out "${prefix}"
  zcat "$vcf" | vcf-annotate --fill-type | grep -oP "TYPE=\\w+" | sort | uniq -c > "${prefix}.summary"
  python /nbt_main/home/lattapol/nextflow-Callvariants/bin/create_AF_his.py ${prefix}.frq
  python /nbt_main/home/lattapol/nextflow-Callvariants/bin/create_lmiss_his.py ${prefix}.lmiss  
  """
}
```
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
