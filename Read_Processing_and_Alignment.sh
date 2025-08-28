#Download Endo RNA-seq Data
#mkdir /media/storageA/hani/Apo_LUTI_Project/RNA_Endoderm_1_2
#cd /media/storageA/hani/Apo_LUTI_Project/RNA_Endoderm_1_2
#aws s3 sync s3://apo.endorep.1and2.rnaseq.pairedend.data/JK12_fastqfiles/ ./
#mkdir raw_fastq
#mv *gz ./raw_fastq
#sudo mkdir /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/RNA_Endoderm_1_2
#sudo rsync -avh --progress --update --remove-source-files ./ /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/RNA_Endoderm_1_2
#ln -s /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/RNA_Endoderm_1_2/* ./

#Download Endo Ribo-seq Data: Replicate 1
#mkdir /media/storageA/hani/Apo_LUTI_Project/Ribo_Endoderm_1
#cd /media/storageA/hani/Apo_LUTI_Project/Ribo_Endoderm_1
#mkdir raw_fastq
#cd raw_fastq/
#aws s3 sync s3://apo.endo.rep1.riboseq.data/JK09_JD04_fastqfiles/ ./
#mv Endoderm_Rep1_Riboseq_SampleInfo.xlsx ../
#sudo mkdir /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/Ribo_Endoderm_1
#sudo rsync -avh --progress --update --remove-source-files ./ /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/Ribo_Endoderm_1/
#ln -s /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/Ribo_Endoderm_1/* ./

#Download Endo Ribo-seq Data: Replicate 2
#mkdir /media/storageA/hani/Apo_LUTI_Project/Ribo_Endoderm_2
#mkdir /media/storageA/hani/Apo_LUTI_Project/Ribo_Endoderm_2/raw_fastq
#cd /media/storageA/hani/Apo_LUTI_Project/Ribo_Endoderm_2/raw_fastq
#aws s3 sync s3://apo.endo.rep2.riboseq.data/JK07_08_fastqfiles/ ./
#mv Endoderm_Rep2_Riboseq_SampleInfo.xlsx ../
#sudo mkdir /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/Ribo_Endoderm_2
#ln -s /mnt/Engram_Locker/fastq/hani/Apo_LUTI_Project/Ribo_Endoderm_2/* ./

#Check that the files look ok on Engram & the links ok on storage A. Then Align everything


#STAR RNA Alignment: Endoderm
#Go back and align all the paired-end RNA-seq for Ribosome Profiling using STAR
#LibType = ISR which is the same as the Takara Kit I use. 
#Align read pairs, then extract later from Bam 
cd /media/storageA/hani/Apo_LUTI_Project/RNA_Endoderm_1_2
array=()
while IFS=  read -r -d $'\0'; do
    array+=("$REPLY")
done < <(find ./raw_fastq -regex ".*R1.fastq.gz" -print0) #cant use -type f with symlinks.

for i in ${array[@]}; do
  echo Identified File $i
  my_path=$(echo $i | grep -Po ".*(?=_R1.fastq.gz)")
  my_name=$(echo $i | grep -Po "(?<=raw_fastq/JK12_)BC[0-9]+(?=_R1.fastq.gz)")
  echo Extracted name $my_name
  echo Created directory ./STARout_RNA_$my_name 
  mkdir ./STARout_RNA_$my_name
  echo Initated Star Mapping to ./STARout_RNA_$my_name"/RNA_"$my_name"_"
  STAR --runMode alignReads --runThreadN 8 --genomeDir /media/storageA/hani/alignment/hg38_annotations/STARDir_hg38_Coding75nt --outFileNamePrefix ./STARout_RNA_$my_name"/RNA_"$my_name"_" --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMstrandField intronMotif --readFilesCommand zcat --outFilterMultimapNmax 10 --readFilesIn $i $my_path"_R2.fastq.gz"
done

rna_bams=()
while IFS=  read -r -d $'\0'; do
    rna_bams+=("$REPLY")
done < <(find ./ -iregex .*STARout_.*bam\$ -type f -print0)

for i in ${rna_bams[@]}; do
    echo Working on file $i
    my_path=$(echo $i | grep -Po ".*(?=_Aligned.sortedByCoord.out.bam)")
    samtools view -b -q 30 -f 129 --threads 20 -o $my_path"_R2only_q30.bam" $i
    echo Indexing file $my_path"_R2only_q30.bam"
    samtools index $my_path"_R2only_q30.bam"
    echo Generating wiggle files with output $my_path"_R2only_q30"
    make_wiggle --count_files $my_path"_R2only_q30.bam" --countfile_format BAM --min_length 17 --fiveprime -o $my_path"_R2only_q30"
    echo Generating bigwigs...
    wigToBigWig $my_path"_R2only_q30_fw.wig" /media/storageA/hani/alignment/hg38_annotations/hg38.chrom.sizes $my_path"_R2only_q30_fw.wig.bw"
    wigToBigWig $my_path"_R2only_q30_rc.wig" /media/storageA/hani/alignment/hg38_annotations/hg38.chrom.sizes $my_path"_R2only_q30_rc.wig.bw"
    echo Removing wiggles...
    rm $my_path"_R2only_q30_fw.wig"
    rm $my_path"_R2only_q30_rc.wig"
    echo Running featurecounts...
    custom_featurecounts.py --run_mode rnaseq --transcript_annotations /media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly.obj --transcript_format Obj --path_to_rnaseq $my_path"_R2only_q30.bam" --rnamapfactory fiveprime --nthreads 22 --outfile_prefix $my_path"_rnacounts"
    echo '...done!'
done

#running in screen, then finish up!

#####################################
##### Ribo Profiling Alignments #####
#####################################
#Basic Pipeline. Clipped first 3 reads and adapter trim polyA tails. NOT requiring polyA tail for now. Bowtie and star still allowing 0 mismatches. 
function riboprof_basic_align {
    echo Received Input Fastq: $1
    prefix=$(echo $1 | grep -Po ".*(?=raw_fastq)")
    name=$(echo $1 | grep -Po "(?<=raw_fastq/).*(?=.fastq.gz)")
    echo Parsed name $name
    echo Pre-processing input...
    cutadapt_log_out=$prefix"cutadaptlog_"$name".log"
    unalig_out=$prefix$name"_unalig.fq"
    bowtie_log_out=$prefix"bowtielog_"$name".log"
    echo "Aligning..."
    echo "Cutadapt Log Out: "$cutadapt_log_out
    echo "Unaligned Fq: "$unalig_out
    echo "Bowtie Log Out: "$bowtie_log_out
    zcat $1 | cutadapt --cut 3 -a "A{15}" -j 8 --nextseq-trim=20 --minimum-length 17 - 2>$cutadapt_log_out | bowtie -p $2 -v 0 --un $unalig_out /media/storageA/hani/alignment/hg38_45srRNA_bowtie/hg38_rRNA - >/dev/null 2>$bowtie_log_out
    echo "...done!"
    
    echo Compressing $unalig_out using $2 cores...
    pigz $unalig_out -p $2
    echo "...done!"
    
    star_out_dir=$prefix"STARout_Ribo_"$name
    echo "Starting Mapping..."
    echo "STAR Out Dir: "$star_out_dir
    star_out_prefix=$star_out_dir"/Ribo_"$name"_"
    mkdir $star_out_dir
    echo "STAR Out Prefix: "$star_out_prefix
    STAR --runMode alignReads --runThreadN $2 --genomeDir /media/storageA/hani/alignment/hg38_annotations/STARDir_hg38_Coding30nt --outFileNamePrefix $star_out_prefix --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 400 --outFilterMismatchNmax 0 --outFilterMatchNmin 15 --outFilterMultimapNmax 1 --readFilesCommand zcat --readFilesIn $unalig_out".gz"
    echo "...done!"
    echo "........................................."
}

#Align all the Ribo
array=()
while IFS=  read -r -d $'\0'; do
    array+=("$REPLY")
done < <(find ./Ribo_Endoderm_1/raw_fastq ./Ribo_Endoderm_2/raw_fastq ./Ribo_Mesoderm_1/raw_fastq ./Ribo_Mesoderm_2/raw_fastq -regex ".*.fastq.gz" -print0)

iter_=1
x=$(echo ${#array[@]})
for i in ${array[@]}
do
    echo "Iter: "$iter_"/"$x
    riboprof_basic_align $i 8
    iter_=$((iter_+1))
done

#for each sub folder:
#1) Moved all logs to ./logs
#2) Moved unalig files to proper folder on Engram & symlinked back to storageA

############################################
##### QC/Processing all Endo Ribo #####
############################################

#Search list
for_postproc=()
while IFS=  read -r -d $'\0'; do
    for_postproc+=("$REPLY")
done < <(find ./Ribo_Endoderm_1 ./Ribo_Endoderm_2 ./Ribo_Mesoderm_1 ./Ribo_Mesoderm_2 -type f -iname "*Aligned.sortedByCoord.out.bam" -print0)

#Pipeline. did first manually to test it, hence the array slice here.
for i in ${for_postproc[@]}
do
    echo "Found File: $i"
    echo "Indexing..."
    samtools index $i
    my_path=$(echo $i | grep -Po ".*(?=_Aligned.sortedByCoord.out.bam)")
    echo "Running psite script. Outfile Name Set as $my_path"
    psite --count_files $i --countfile_format BAM --min_counts 50 --require_upstream --min_length 17 --max_length 35 /media/storageA/hani/alignment/hg38_annotations/hg38_gencode_plastid_metagene/hg38_CDSStart_50bp_rois.txt $my_path
    read -p "Please edit the P Offset file located in folder $my_path, then press enter to continue"
    echo "Running Metagene script. Outfile Name Set as $my_path"
    metagene count --count_files $i --fiveprime_variable --offset $my_path"_p_offsets.txt" --normalize_over 20 50 --min_counts 50 /media/storageA/hani/alignment/hg38_annotations/hg38_gencode_plastid_metagene/hg38_CDSStart_50bp_rois.txt $my_path
    echo "Running Phase by Size Script. Outfile Name Set as $my_path"
    phase_by_size --count_files $i --countfile_format BAM --fiveprime_variable --offset $my_path"_p_offsets.txt" /media/storageA/hani/alignment/hg38_annotations/hg38_gencode_plastid_metagene/hg38_CDSStart_50bp_rois.txt $my_path
    echo "...done!"
done

#######################################
##### FeatureCounts all Endo #####
#######################################

#Run TE/FeatureCounts
path_tx="/media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly.obj"
path_uORF="/media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly_2codon_uorfs.obj"

{
  while IFS=$'\t' read -ra line
  do
    echo Working on ${line[0]}:
    echo "(RNA): "${line[1]}
    echo "(Rib): "${line[2]}
    echo "(Psi): "${line[3]}
    echo '.....................'
    
    
    custom_featurecounts.py --run_mode both --transcript_annotations $path_tx --transcript_format Obj --uORF_annotations $path_uORF --uORF_format Obj --path_to_ribprof_bam ${line[2]} --path_to_psite ${line[3]} --path_to_rnaseq ${line[1]} --rnamapfactory fiveprime --nthreads 8 --outfile_prefix /media/storageA/hani/Apo_LUTI_Project/FeatureCountsOut_Endo_Meso/${line[0]}"_TEcounts"
    echo ...done
  done
} < /media/storageA/hani/Apo_LUTI_Project/FeatureCountsOut_Endo_Meso/featurects_params.tsv
