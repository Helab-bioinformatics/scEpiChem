mkdir raw
mv *fq.gz raw/

cd raw/
function umi_tools_whitelist(){
	file="$1"
	base=$(basename $file "_1.fq.gz")
	umi_tools whitelist --stdin $file --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8}).*'  --extract-method=regex --set-cell-number 1000 --ed-above-threshold=correct  --error-correct-threshold=2 --log2stderr --knee-method=distance > ${base}_whitelist.txt --plot-prefix ${base}
}
export -f umi_tools_whitelist
parallel -j30 umi_tools_whitelist ::: ./*_1.fq.gz
cd ..
mkdir 01_extract


cd raw/

function umi_tools_extract () {
    file="$1"
    base=$(basename $file "_whitelist.txt")
    umi_tools extract --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8}).*' --stdin ${base}_1.fq.gz  --stdout ../01_extract/${base}_1.extract.fq.gz  --read2-in ${base}_2.fq.gz  --read2-out ../01_extract/${base}_2.extract.fq.gz --error-correct-cell --extract-method=regex --whitelist ${base}_whitelist.txt
}
export -f umi_tools_extract
parallel -j30 umi_tools_extract ::: ./*_whitelist.txt


cd ../01_extract/
function trim_cutadapt() {
    file="$1"
    base=$(basename $file "_2.extract.fq.gz")
    cutadapt -u 52 -o ${base}.trimmed.2.extract.fq.gz ${base}_2.extract.fq.gz
}
function apply_cutadapt() {
    file="$1"
    base=$(basename $file "_2.extract.fq.gz")
    cutadapt -q 20 -O 10 -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 10 --max-n 0.1 --trim-n -o ${base}.cut.1.fq.gz -p ${base}.cut.2.fq.gz ${base}_1.extract.fq.gz ${base}.trimmed.2.extract.fq.gz
}
function extract_1() {
    file="$1"
    base=$(basename $file "_1.extract.fq.gz")
    zcat ${base}.cut.1.fq.gz | awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' - > ${base}_1.extract.fa
}
function extract_2() {
    file="$1"
    base=$(basename $file "_2.extract.fq.gz")
    zcat ${base}.cut.2.fq.gz | awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' - > ${base}_2.extract.fa
}
export -f trim_cutadapt apply_cutadapt extract_1 extract_2

parallel -j30 trim_cutadapt ::: *_2.extract.fq.gz
parallel -j30 apply_cutadapt ::: *_2.extract.fq.gz
parallel -j30 extract_1 ::: *_1.extract.fq.gz
parallel -j30 extract_2 ::: *_2.extract.fq.gz


cd ..
mkdir 02_align
function align_to_hg19() {
    file="$1"
    base=$(basename $file "_2.extract.fa")
    bowtie2 -p 10 -f -N 1 --very-sensitive-local --no-unal -x /media/helab/data1/00_public/database/genomes/hg19_bowtie2/hg19 -1 ./01_extract/${base}_1.extract.fa -2 ./01_extract/${base}_2.extract.fa -S 02_align/${base}.hg19.sam 2> 02_align/${base}.hg19.align.log
}
export -f align_to_hg19
parallel -j30 align_to_hg19 ::: ./01_extract/*_2.extract.fa


mkdir 03_bam
function process_sam() {
    file="$1"
    base=$(basename $file ".sam")
    samtools view -hbS -q 20 $file | samtools sort -T ${base} - > 03_bam/${base}.bam
}

export -f process_sam
parallel -j30 process_sam ::: 02_align/*.sam


cd 02_align
rm *sam
cd 01_extract
rm *fa
cd ..

mkdir 04_rmdup
function mark_duplicates() {
    i="$1"
    base=$(basename $i ".bam")
    java -jar /media/helab/data1/00_public/software/picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=./03_bam/${base}.bam o=./04_rmdup/${base}_rmdup.bam M=./04_rmdup/${base}_rmdup_picard.txt
}
export -f mark_duplicates
parallel -j30 mark_duplicates ::: 03_bam/*.bam

function process_rmdup() {
    i="$1"
    base=$(basename $i "_rmdup.bam")
    samtools view ./04_rmdup/${base}_rmdup.bam | sed 's/_/\tCB:Z:/'  - | sed 's/_/\tTB:Z:/' - |awk '{print$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$2"\t"$3}' -> ./04_rmdup/${base}_CB_TB.bam
}
export -f process_rmdup
parallel -j30 process_rmdup ::: ./04_rmdup/*_rmdup.bam

function process_tags() {
    i="$1"
    base=$(basename $i "_CB_TB.bam")
    samtools view ./04_rmdup/${base}_rmdup.bam  -H | cat - ./04_rmdup/${base}_CB_TB.bam | samtools sort -t CB - -o ./04_rmdup/${base}_CB_TB_sorted_tags.bam
}
export -f process_tags
parallel -j30 process_tags ::: ./04_rmdup/*_CB_TB.bam


cd 04_rmdup/
rm *_CB_TB.bam
cd ..

mkdir 05_split







