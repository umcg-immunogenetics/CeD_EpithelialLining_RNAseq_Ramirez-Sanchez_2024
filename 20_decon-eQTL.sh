
srun --cpus-per-task=1 --mem=16gb --nodes=1 --qos=priority --time=00:59:00 --pty bash -i

cat shared.samplenames.new.line.txt > genotypes_decon-eQTL.txt
echo "" >> genotypes_decon-eQTL.txt
awk -F '\t' '(NR>1)' intersected.SNPs.sharedSamples.genotypes.matrix.dosage.txt >> genotypes_decon-eQTL.txt


#probes
awk -F '\t' '(NR>1) {print $1}' RNA_decon_mcc.table > only_genes.txt 

awk 'FNR==NR { a[$NF]; next } ($NF in a)' only_genes.txt \
/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/references/f5.bed > Expressed.Homo_sapiens.GRCh37.75.genes.f5.bed


awk -F "\t" '{print $1"\t"int($2+(($3-$2+1)/2)-1)"\t"int($3-(($3-$2+1)/2)+1)"\t"$4"\t"$5}' Expressed.Homo_sapiens.GRCh37.75.genes.f5.bed > Center.Expressed.Homo_sapiens.GRCh37.75.genes.f5.bed

ml BEDTools

bedtools window \
-a /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/references/Top_CeD_SNPs_May2019.bed \
-b Expressed.Homo_sapiens.GRCh37.75.genes.f5.bed \
-w 250000 > Expressed_CeD_SNPs_Genes_0.5MB.bed

#zcat intersection_SNPs_CeDgenes_center.vcf.gz | awk -F "\t" '{ print $5"\t"$8}' > snpsToTest_CeD_Center.txt

awk -F "\t" '{print $11"\t"$4}' Expressed_CeD_SNPs_Genes_0.5MB.bed > CeDSNPsGenesExpressed.txt

awk -F "\t" '!seen[$2]++' CeDSNPsGenesExpressed.txt > Only.CeDSNPsGenesExpressed.txt
awk -F "\t" '(NR>1)' Only.CeDSNPsGenesExpressed.txt > Only.CeDSNPsGenesExpressed2.txt
awk -F "\t" '{print $2}' Only.CeDSNPsGenesExpressed2.txt > Only.CeDSNPsGenesExpressed3.txt
#intersecting those unqiue snps with dosage 
working_dir=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/2020-Epithelial_Oslo_Deconvolution/ongoing

awk -F "\t" 'FNR==NR { a[$1]; next } ($1 in a)' Only.CeDSNPsGenesExpressed3.txt \
${working_dir}/genotypes_for_deconeQTL/sharedSamples.genotypes.matrix.dosage > intersected.SNPs.sharedSamples.genotypes.matrix.dosage 
#getting list of SNPs in vcf and probe
awk -F "\t" '{print $1}' intersected.SNPs.sharedSamples.genotypes.matrix.dosage > Only.CeDSNPs.inVCF.inRNA.txt

awk -F "\t" 'FNR==NR { a[$NF]; next } ($NF in a)' Only.CeDSNPs.inVCF.inRNA.txt \
CeDSNPsGenesExpressed.txt > CeDSNPsGenesExpressedGenotypes.txt
sed  -i '1i gene\tsnp' CeDSNPsGenesExpressedGenotypes.txt 





#!/bin/bash
#SBATCH --job-name=DeconeQTL_models
#SBATCH --output=../job_log/%x.out
#SBATCH --error=../job_log/%x.err
#SBATCH --time=1:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --nodes=1
ml Java

working_dir=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/2020-Epithelial_Oslo_Deconvolution/ongoing/eQTL
#cell counts:
FACS2=${working_dir}/input_decon-eQTL/FACS/cell_proportions_scaled.2groups.txt
FACS6=${working_dir}/input_decon-eQTL/FACS/cell_proportions_scaled.6groups.txt

#RNA expression
RNAm=${working_dir}/input_decon-eQTL/gene_expression/RNA_decon_marsh.table
RNAcr=${working_dir}/input_decon-eQTL/gene_expression/RNA_decon_crypt.table
RNAco=${working_dir}/input_decon-eQTL/gene_expression/RNA_decon_condition.table
RNAmcc=${working_dir}/input_decon-eQTL/gene_expression/RNA_decon_mcc.table

#genotypes:
GENOTYPES=${working_dir}/genotypes_for_deconeQTL/genotypes_decon-eQTL.txt

#probes
PROBES=${working_dir}/input_decon-eQTL/gene_expression/CeDSNPsGenesExpressedGenotypes.txt 


#output:
FACS2_RNAm=${working_dir}/eQTL/decon-eQTL/FACS2_RNAm
FACS2_RNAcr=${working_dir}/eQTL/decon-eQTL/FACS2_RNAcr
FACS2_RNAco=${working_dir}/eQTL/decon-eQTL/FACS2_RNAco
FACS2_RNAmcc=${working_dir}/eQTL/decon-eQTL/FACS2_RNAmcc

FACS6_RNAm=${working_dir}/eQTL/decon-eQTL/FACS6_RNAm
FACS6_RNAcr=${working_dir}/eQTL/decon-eQTL/FACS6_RNAcr
FACS6_RNAco=${working_dir}/eQTL/decon-eQTL/FACS6_RNAco
FACS6_RNAmcc=${working_dir}/eQTL/decon-eQTL/FACS6_RNAmcc


mkdir ${FACS2_RNAm}
mkdir ${FACS2_RNAcr}
mkdir ${FACS2_RNAco}
mkdir ${FACS2_RNAmcc}

mkdir ${FACS6_RNAm}
mkdir ${FACS6_RNAcr}
mkdir ${FACS6_RNAco}
mkdir ${FACS6_RNAmcc}

java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS2} \
                        -e ${RNAm} \
                        -g ${GENOTYPES} \
                        -o ${FACS2_RNAm} \
                        -sn ${PROBES} \
                        -r

java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS2} \
                        -e ${RNAcr} \
                        -g ${GENOTYPES} \
                        -o ${FACS2_RNAcr} \
                        -sn ${PROBES} \
                        -r
java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS2} \
                        -e ${RNAco} \
                        -g ${GENOTYPES} \
                        -o ${FACS2_RNAco} \
                        -sn ${PROBES} \
                        -r
java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS2} \
                        -e ${RNAmcc} \
                        -g ${GENOTYPES} \
                        -o ${FACS2_RNAmcc} \
                        -sn ${PROBES} \
                        -r


java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS6} \
                        -e ${RNAm} \
                        -g ${GENOTYPES} \
                        -o ${FACS6_RNAm} \
                        -sn ${PROBES} \
                        -r
java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS6} \
                        -e ${RNAcr} \
                        -g ${GENOTYPES} \
                        -o ${FACS6_RNAcr} \
                        -sn ${PROBES} \
                        -r
java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS6} \
                        -e ${RNAco} \
                        -g ${GENOTYPES} \
                        -o ${FACS6_RNAco} \
                        -sn ${PROBES} \
                        -r
java -jar /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/tools/decon-eQTL/Decon-eQTL-v1.4-jar-with-dependencies.jar \
                        -c ${FACS6} \
                        -e ${RNAmcc} \
                        -g ${GENOTYPES} \
                        -o ${FACS6_RNAmcc} \
                        -sn ${PROBES} \
                        -r
