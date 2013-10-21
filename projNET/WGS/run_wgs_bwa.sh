#success!
#running script: qsub -N name -e logs -o logs script
mkdir -p /ifs/scratch/c2b2/ys_lab/jh3283/wgstmp_ref
# index ref genome
ref_hm="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
echo $ref_hm
bwa index -a bwtsw $ref_hm

