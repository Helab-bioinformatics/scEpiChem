for i in `ls *_R1.fastq.gz`
do 
base=$(basename $i "_R1.fastq.gz")
#$base,base,${base}
#/media/helab/data1/00_public/database/genomes/E_coli_bowtie2/E_coli
#/media/helab/data1/00_public/database/genomes
bowtie2 -p 3 -N 1 --dovetail --very-sensitive-local --no-unal --no-mixed --no-discordant -X 2000 -x /media/helab/data1/00_public/database/genomes/hg19_bowtie2/hg19 -1 ${base}_R1.fastq.gz -2 ${base}_R2.fastq.gz -S ${base}.sam 2> ${base}.align.log
done

for i in `ls *.sam`
do 
base=$(basename $i ".sam")

samtools view -hbS -q 30 ${base}.sam | samtools sort -T ${base}  - >${base}.bam
done

for i in ./*bam
do 
samtools view -c $i >./${i%.*}.totolreads.log
done

rm *sam

for j in *bam
do 
base=$(basename $j ".bam")
java -jar /media/helab/data1/00_public/software/picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${base}.bam o=${base}.rmdup.bam M=${base}.picard.txt
done

for i in *.rmdup.bam
do 
samtools index $i
done

for i in ./*.rmdup.bam
do
base=$(basename $i ".rmdup.bam")
num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
bamCoverage --scaleFactor  $res -b  $i   -o   ./${base}.10M.bw -p 3
bamCoverage --scaleFactor  $res -b  $i  -e 300  --smoothLength 500 -o  ./${base}.ext300.smo500.bw -p 5
macs2 callpeak -t $i -f BAM -g hs -n ${base} --nolambda --nomodel --broad
done


mkdir raw
mkdir 02_align
mkdir 03_bam
mkdir 04_rmdup
mkdir 05_bw

mv *.gz raw
mv *.align.log 02_align
mv *.rmdup.bam 04_rmdup
mv *.bai 04_rmdup
mv *.picard.txt 04_rmdup
mv *.bam 03_bam
mv *totolreads.log 03_bam
mv *.bw 05_bw











