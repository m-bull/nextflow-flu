#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


if( params.help ) {
    return help()
}
if( !nextflow.version.matches('0.25+') ) {
    return nextflow_version_error()
}
if( params.index ) { 
    index = Channel.fromPath(params.index).toSortedList() 
    if( index.size() == 0 ) return index_error(index)
}
if( params.host ) {     
    host = file(params.host)
    if( !host.exists() ) return host_error(host)
}
if( params.adapters ) {     
    adapters = file(params.adapters) 
    if( !adapters.exists() ) return adapter_error(adapters)
}
if( params.fqc_adapters ) {
    fqc_adapters = file(params.fqc_adapters)                             
    if( !fqc_adapters.exists() ) return fastqc_error(fqc_adapters)
}

threads = params.threads

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { return fastq_error(params.reads) }
    .into { read_pairs; fastqc_pairs; alignment_pairs }

process RunFastQC {
    tag { dataset_id }

    publishDir "${params.output}/RunFastQC", mode: 'copy'

    input:
        set dataset_id, file(forward), file(reverse) from fastqc_pairs

    output:
        set dataset_id, file("*_fastqc.zip") into (fastqc_logs)

    """
    mkdir output
    fastqc -f fastq ${forward} ${reverse} -t ${threads} -o output
    mv output/*.zip .
    """
}

process RunTrim {
    tag { dataset_id }

    publishDir "${params.output}/RunTrim", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (filename.indexOf(".fq.gz") > 0) "trimmed/$filename"
        }

    input:
    set dataset_id, file(forward), file(reverse) from read_pairs

    output:
    set dataset_id, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmed_reads

    script:
        """
        trim_galore --fastqc --paired --gzip $forward $reverse 
        """
}

if( !params.index ) {
    process BuildHostIndex {
        tag { host.baseName }

        publishDir "${params.output}/BuildHostIndex", mode: "copy"

        input:
            file(host)

        output:
            file '*' into index

        """
        bwa index ${host}
        """
    }
}

if( !params.index ) {
    process AlignReadsToHost {
        tag { host.baseName }
        
        publishDir "${params.output}/AlignReadsToHost", mode: "copy"
        
        input:
            set dataset_id, file(forward), file(reverse) from trimmed_reads
            file host
            
        output:
            set dataset_id, file("${dataset_id}.host.sam") into host_sam
            
        """ 
        bwa mem ${host} ${forward} ${reverse} -t ${threads} > ${dataset_id}.host.sam
        """
    }
}

if( params.index ) {
    process AlignReadsToHost {
        tag { dataset_id }

        publishDir "${params.output}/AlignReadsToHost", mode: "copy"

        input:
            set dataset_id, file(forward), file(reverse) from trimmed_reads
            file idx from index

        output:
            set dataset_id, file("${dataset_id}.host.sam") into host_sam

        """
        bwa mem ${idx.first()} ${forward} ${reverse} -t ${threads} > ${dataset_id}.host.sam
        """
    }
}

process RemoveHostDNA {
    tag { dataset_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy"

    input:
        set dataset_id, file(sam) from host_sam

    output:
        set dataset_id, file("${dataset_id}.host.sorted.removed.bam") into host_bam

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${dataset_id}.host.sorted.bam
    samtools view -h -f 4 -b ${dataset_id}.host.sorted.bam -o ${dataset_id}.host.sorted.removed.bam
    """
}

process BAMToFASTQ {
    tag { dataset_id }

    publishDir "${params.output}/BAMToFASTQ", mode: "copy"

    input:
        set dataset_id, file(bam) from host_bam

    output:
        set dataset_id, file("${dataset_id}.non.host.R1.fastq"), file("${dataset_id}.non.host.R2.fastq") into non_host_fastq_iva, non_host_fastq_variants

    """
    bedtools  \
       bamtofastq \
      -i ${bam} \
      -fq ${dataset_id}.non.host.R1.fastq \
      -fq2 ${dataset_id}.non.host.R2.fastq
    """
}

process RunIVA {
    tag { dataset_id }

    publishDir "${params.output}/RunIVA", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse) from non_host_fastq_iva

    output:
        set dataset_id, file("${dataset_id}.contigs.fa") into (iva_contigs)

    script:
    """
    iva \
      -t ${threads} \
      -f ${forward} \
      -r ${reverse} \
      output

    mv output/contigs.fasta .
    # Remove line breaks from FASTA files
    # Taken from Kent at https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file
    awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' contigs.fasta > ${dataset_id}.contigs.fa
    """
}

process RunBlast {
    tag { dataset_id }

    publishDir "${params.output}/RunBlast", mode: 'copy'

    input:
        set dataset_id, file(contigs) from iva_contigs

    output:
        file("${dataset_id}.contigs.annotated.fa") into (annotated_iva_contigs_variantref, annotated_iva_contigs_variantcall)

    script:
    """
    blastn -db InfluenzaDB -query ${contigs} -max_hsps 1 -max_target_seqs 1 -outfmt "10 stitle" -num_threads ${threads} -task megablast > ${dataset_id}.contig.description.tmp
    cat ${dataset_id}.contig.description.tmp | sed -e '/Influenza/s/^/>/' > ${dataset_id}.contig.description.txt
    sed -i 's/ /_/g' ${dataset_id}.contig.description.txt
    awk '/^>/ { getline <"${dataset_id}.contig.description.txt" } 1 ' ${contigs} > ${dataset_id}.contigs.annotated.fa
    """
}


process PreparePHEnixRef {
    tag { dataset_id }

    publishDir "${params.output}/PHEnixRef", mode: "copy"

    input:
        file(contigs) from annotated_iva_contigs_variantref

    output:
        file '*' into phenixref

    """
    phenix.py prepare_reference --reference ${contigs} --mapper bwa --variant gatk
    """
}


process PHEnixVariants {
    tag { dataset_id }

    publishDir "${params.output}/PHEnixVariants", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse) from non_host_fastq_variants
        file '*' from phenixref
        file(contigs) from annotated_iva_contigs_variantcall

    output:
        set dataset_id, file("${dataset_id}.filtered.vcf"), file("${dataset_id}.vcf.idx") into phenixvcf

    """
    phenix.py run_snp_pipeline -r1 ${forward} -r2 ${reverse} --reference ${contigs} --mapper bwa --variant gatk --sample-name ${dataset_id} --mapper-options "-t ${threads}" --filters "mq_score:30,min_depth:10,ad_ratio:0.2" --outdir .
    """
}

process PHEnixVariantsToFasta {
    tag { dataset_id }

    publishDir "${params.output}/PHEnixVariantsFasta", mode: "copy"

    input:
        set dataset_id, file(vcf), file(vcf_idx) from phenixvcf
        
    output:
        file("${dataset_id}.variants.fasta") into phenixvarfasta
        file("${dataset_id}.variants.fasta.stats.txt") into phenixfastastats

    """
    phenix.py vcf2fasta --input ${vcf} --out ${dataset_id}.variants.fasta --with-stats ${dataset_id}.variants.fasta.stats.txt
    """
}


multiQCReports = Channel.empty()
    .mix(
        fastqc_logs
    )
    .flatten().toList()

process RunMultiQC {
    publishDir "${params.output}/RunMultiQC", mode: 'copy'

    input:
        file('*') from multiQCReports

    output:
        file("*multiqc_report.html") into multiQCReport

    """
    multiqc -f -v .
    """
}

def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastqc_error(def input) {
    println ""
    println "[params.fqc_adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def host_error(def input) {
    println ""
    println "[params.host] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def index_error(def input) {
    println ""
    println "[params.index] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def help() {
    println ""
    println "Program: flu.nf"
    println "Version: $workflow.repository - $workflow.revision [$workflow.commitId]"
    println ""
    println "Usage:    nextflow run flu.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --adapters      STR      path to FASTA formatted adapter sequences"
    println "    --host          STR      path to FASTA formatted host genome"
    println "    --index         STR      path to BWA generated index files"
    println "    --output        STR      directory to write process outputs to"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of threads to use for each process"
    println "    --min-alt-count INT      requires this number of observations supporting an alternate allele"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
