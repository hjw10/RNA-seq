dir1='/public/workspace/stu21230110/myNGS/bigwork'

star_index='/public/workspace/lincs/lab7/ref/star_271a'
hisat2_index='/public/workspace/stu21230110/myNGS/bigwork/grch37_tran'
gtf_dir='/public/workspace/lincs/lab7/ref/Homo_sapiens.GRCh38.98.gtf'
genom_dir='/public/workspace/lincs/lab7/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
salmon_dir='/public/workspace/lincs/lab7/ref/grch38_salmon_index'
rsem_dir='/public/workspace/lincs/lab7/ref/rsem-star-index/rsem-star-index'

#GEO数据库找到项目编号，再到ENA数据库找数据下载
#将SRR_Acc_List.txt放进来
prefetch -O /public/workspace/stu21230110/myNGS/bigwork/data/ --option-file SRR_Acc_List.txt
#sra转fastq
#fastq-dump --gzip --split-files ../data/SRR14112694/SRR14112694.sra
#......
#或者批量处理 
fastq-dump --gzip --split-files ../data/SRR*/*
#处理前做一次fastqc
for file in *.fastq.gz
do
fastpc $file &
done
wait
fastqc .

#fastp做预处理，这里还是批量操作
mkdir -p 2.cleandata
in_dir=1.fastqfile
out_dir=2.cleandata
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
	echo $name
	echo `date`
	echo "** Starting to clean $name **"
	fastp -i $in_dir/$name'_1.fastq.gz' -o $out_dir/$name'_1.fastq.gz' \
          -I $in_dir/$name'_2.fastq.gz' -O $out_dir/$name'_2.fastq.gz' \
		  -w 10
done
#处理后看一下fastqc
cd 2.cleandata
for file in *.fastq.gz
do
    fastqc $file &
done
wait
multiqc --filename clean.html .
cd ..


#进行比对
#文献中用Tophat +htseq 太慢，我的流程一：HISAT2 + HTSeq
#数据是hg19所以先下索引文件http://daehwankimlab.github.io/hisat2/download/

mkdir 3.hisat2_sam
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
	echo $name
	echo `date`
	echo "** Starting to HISAT2 $name **"
	singularity exec /public/workspace/lincs/lab7/soft/hisat2_2.1.0.sif hisat2 \
	--dta -p 8 -x $hisat2_index/genome_tran \
	-1 $dir1/2.cleandata/$name'_1.fastq.gz' \
	-2 $dir1/2.cleandata/$name'_2.fastq.gz' \
	-S $dir1/3.hisat2_sam/$name'_hisat2.sam'
done


#sam转bam并且排序
mkdir -p 4.hisat2_bam
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
	echo $name
	echo `date`
	echo "** Starting to sort $name **"
	singularity exec /public/workspace/lincs/lab7/soft/samtools_latest.sif \
	samtools sort -@ 8 -o $dir1/4.hisat2_bam/$name'_hisat2.bam' \
    $dir1/3.hisat2_sam/$name'_hisat2.sam'
done


#HTSeq 计数
mkdir -p 5.hisat2_htseq
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
	echo $name
	echo `date`
	echo "** Starting to HTSeq $name **"
	singularity exec /public/workspace/lincs/lab7/soft/htseq_count_0.11.2.sif \
	htseq-count -f bam -r pos -s no -a 16 -t exon -i gene_id -m intersection-nonempty \
	$dir1/4.hisat2_bam/$name'_hisat2.bam' \
	$dir1/Homo_sapines.GRCh38.98.gtf > $dir1/5.hisat2_htseq/$name'_hisat2_counts.txt'
done

#流程2：STAR + FeatureCounts
#STAR比对
mkdir -p 6.star
#因为数据是GRCh37的，所以先自己构建索引
#从encode 网站下载参考文件
#注释文件: wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz
#参考基因组: wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
#构建索引但是后面不用
mkdir star_271
time singularity exec /public/workspace/lincs/lab7/soft/star_2.7.1a.sif \
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $dir1/star_271 \
--genomeFastaFiles $genom_dir \
--sjdbGTFfile $gtf_dir

#star比对
cd $dir1/6.star
cp /public/workspace/lincs/lab7/ref/star_271a/* /star_271a/
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
	echo $name
	echo `date`
	echo "** Starting to STAR $name **"
	singularity exec /public/workspace/lincs/lab7/soft/star_2.7.1a.sif \
	STAR --runThreadN 10 --runMode alignReads --genomeDir star_271a/ \
	--readFilesCommand gunzip -c \
	--readFilesIn $dir1/2.cleandata/${name}'_1.fastq.gz' $dir1/2.cleandata/${name}'_2.fastq.gz' \
	--outFileNamePrefix ${name}'_2nd_quant_chim' \
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode TranscriptomeSAM --chimSegmentMin 25
done

#salmon
mkdir 8.salmon
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
mkdir $name
done
cd ..
for name in SRR14112694 SRR14112695 SRR14112696 SRR14112697 SRR14112698 SRR14112699 SRR14112700 SRR14112701 SRR14112702 SRR14112703 SRR14112704 SRR14112705
do
	echo $name
	echo `date`
	echo "** Starting to salmon $name **"
	singularity exec /public/workspace/lincs/lab7/soft/salmon_1.5.0.sif \
	salmon quant -i $dir1/grch38_salmon_index \
	-l A -p 20 \
	-1 $dir1/2.cleandata/$name'_1.fastq.gz' \
	-2 $dir1/2.cleandata/$name'_2.fastq.gz' \
	-o $dir1/8.salmon/$name/$name'_salmon'
done










