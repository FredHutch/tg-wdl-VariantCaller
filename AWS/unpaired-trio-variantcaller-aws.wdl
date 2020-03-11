## Consensus variant calling workflow for human panel/PCR-based targeted DNA sequencing.
## Input requirements:
## - Pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile (a Picard tool)
## - - reads are provided in query-sorted order (not 100% sure if this is required as of 6/24/2019)
## - - all reads must have an RG tag
##
## Output Files:
## - recalibrated bam and it's index
## - GATK vcf
## - samtools/bcftools vcf
## - Annovar annotated vcfs and tabular variant list for each variant caller
## - Basic QC stats from bedtools for mean coverage over regions in panel
## 
workflow Panel_BWA_HC_Mutect_Samtools_AnnotatedVariants {
  File batchFile
  Array[Object] batchInfo = read_objects(batchFile)
  String ref_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
  # references such as b37 and hg19.
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File af_only_gnomad
  File af_only_gnomad_index
  # Note:  For Annovar, please reference: Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010
  File annovarTAR
  String annovar_protocols
  String annovar_operation

  # Docker containers this workflow has been designed for
  String GATKDocker = "broadinstitute/gatk:4.1.4.0"
  String bwaDocker = "fredhutch/bwa:0.7.17"
  String bedtoolsDocker = "fredhutch/bedtools:2.28.0" 
  String bcftoolsDocker = "fredhutch/bcftools:1.9"
  String perlDocker = "perl:5.28.0"
  String RDocker = "rocker/tidyverse:3.6.0"

  Int bwaThreads = 16

  ## Specific script GitHub info for consensus
  String githubRepoURL
  String scriptPath
  String githubTag

scatter (job in batchInfo){
  String sampleName = job.omics_sample_name
  String molecularID = job.molecular_id
  File bamLocation = job.bamLocation
  File bedLocation = job.bedLocation

  String base_file_name = sampleName + "_" + molecularID + "." + ref_name

  # Prepare bed file and check sorting
  call SortBed {
    input:
      unsorted_bed = bedLocation,
      ref_dict = ref_dict,
      docker = GATKDocker
  }

  # Convert unmapped bam to interleaved fastq
  call SamToFastq {
    input:
      input_bam = bamLocation,
      base_file_name = base_file_name,
      docker = GATKDocker
  }

  #  Map reads to reference
  call BwaMem {
    input:
      input_fastq = SamToFastq.output_fastq,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      threads = bwaThreads,
      docker = bwaDocker
  }

  # Merge original uBAM and BWA-aligned BAM
  call MergeBamAlignment {
    input:
      unmapped_bam = bamLocation,
      aligned_bam = BwaMem.output_bam,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      docker = GATKDocker
  }

  # Generate the recalibration model by interval and apply it
  call ApplyBaseRecalibrator {
    input:
      input_bam = MergeBamAlignment.output_bam,
      input_bam_index = MergeBamAlignment.output_bai,
      base_file_name = base_file_name,
      intervals = SortBed.intervals,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker = GATKDocker
    }

    call bedToolsQC {
      input: 
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        genome_sort_order = ApplyBaseRecalibrator.sortOrder,
        bed_file = SortBed.sorted_bed,
        base_file_name = base_file_name,
        docker = bedtoolsDocker
    }

    call CollectHsMetrics {
      input: 
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        intervals = SortBed.intervals,
        docker = GATKDocker
    }

    # Generate haplotype caller vcf
    call HaplotypeCaller {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        intervals = SortBed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf,
        docker = GATKDocker
    }

    # Generate mutect2 vcf
    call Mutect2TumorOnly {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        intervals = SortBed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        genomeReference = af_only_gnomad,
        genomeReferenceIndex = af_only_gnomad_index,
        docker = GATKDocker
    }

    # Annotate variants
    call annovar as annotateHaplotype {
      input:
        input_vcf = HaplotypeCaller.output_vcf,
        ref_name = ref_name,
        annovarTAR = annovarTAR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = perlDocker
    }

    # Generate bcftools vcf
    call bcftoolsMpileup {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        sorted_bed = SortBed.sorted_bed,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf,
        docker = bcftoolsDocker
    }

    # Annotate variants
    call annovar as annotateSAM {
      input:
        input_vcf = bcftoolsMpileup.output_vcf,
        ref_name = ref_name,
        annovarTAR = annovarTAR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = perlDocker
    }
    # Annotate variants
    call annovar as annotateMutect {
      input:
        input_vcf = Mutect2TumorOnly.output_vcf,
        ref_name = ref_name,
        annovarTAR = annovarTAR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = perlDocker
    }

  call consensusProcessingR {
    input:
      GATKVars = annotateHaplotype.output_annotated_table,
      MutectVars = annotateMutect.output_annotated_table,
      SAMVars = annotateSAM.output_annotated_table,
      base_file_name = base_file_name,
      githubRepoURL = githubRepoURL,
      githubTag = githubTag,
      scriptPath = scriptPath,
      docker = RDocker
  }

} # End scatter 
# Outputs that will be retained when execution is complete
output {
    Array[File] analysis_ready_bam = ApplyBaseRecalibrator.recalibrated_bam 
    Array[File] analysis_ready_bai = ApplyBaseRecalibrator.recalibrated_bai
    Array[File] GATK_vcf = HaplotypeCaller.output_vcf
    Array[File] SAM_vcf = bcftoolsMpileup.output_vcf
    Array[File] Mutect_Vcf = Mutect2TumorOnly.output_vcf
    Array[File] Mutect_VcfIndex = Mutect2TumorOnly.output_vcf_index
    Array[File] Mutect_AnnotatedVcf = annotateMutect.output_annotated_vcf
    Array[File] Mutect_AnnotatedTable = annotateMutect.output_annotated_table
    Array[File] GATK_annotated_vcf = annotateHaplotype.output_annotated_vcf
    Array[File] GATK_annotated = annotateHaplotype.output_annotated_table
    Array[File] SAM_annotated_vcf = annotateSAM.output_annotated_vcf
    Array[File] SAM_annotated = annotateSAM.output_annotated_table
    Array[File] panelQC = bedToolsQC.meanQC
    Array[File] PicardQC = CollectHsMetrics.picardMetrics
    Array[File] PicardQCpertarget = CollectHsMetrics.picardPerTarget
    Array[File] consensusVariants = consensusProcessingR.consensusTSV
  }
# End workflow
}

#### TASK DEFINITIONS
# annotate with annovar
task annovar {
  File input_vcf
  String ref_name
  File annovarTAR
  String annovar_protocols
  String annovar_operation
  String docker

  String base_vcf_name = basename(input_vcf, ".vcf.gz")
  command {
  set -eo pipefail
  
  tar -xzvf ${annovarTAR}
  
  perl annovar/table_annovar.pl ${input_vcf} annovar/humandb/ \
    -buildver ${ref_name} \
    -outfile ${base_vcf_name} \
    -remove \
    -protocol ${annovar_protocols} \
    -operation ${annovar_operation} \
    -nastring . -vcfinput
  }
  runtime {
    docker: "${docker}"
    memory: "4GB"
    cpu: "2"
  }
  output {
    File output_annotated_vcf = "${base_vcf_name}.${ref_name}_multianno.vcf"
    File output_annotated_table = "${base_vcf_name}.${ref_name}_multianno.txt"
  }
}
# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  File input_bam
  File intervals 
  File input_bam_index
  String base_file_name
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String docker
  command {
  set -e pipefail
  
  gatk --java-options "-Xms4g" \
    BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${base_file_name}.recal_data.csv \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      --intervals ${intervals} \
      --interval-padding 100 

  gatk --java-options "-Xms4g" \
    ApplyBQSR \
      -bqsr ${base_file_name}.recal_data.csv \
      -I ${input_bam} \
      -O ${base_file_name}.recal.bam \
      -R ${ref_fasta} \
      --intervals ${intervals} \
      --interval-padding 100 

  # finds the current sort order of this bam file
  samtools view -H ${base_file_name}.recal.bam|grep @SQ|sed 's/@SQ\tSN:\|LN://g' > ${base_file_name}.sortOrder.txt
  }
  runtime {
    docker: "${docker}"
    memory: "4GB"
    cpu: "4"
  }
  output {
    File recalibrated_bam = "${base_file_name}.recal.bam"
    File recalibrated_bai = "${base_file_name}.recal.bai"
    File sortOrder = "${base_file_name}.sortOrder.txt"
  }
}
# bcftools Mpileup variant calling
task bcftoolsMpileup {
  File input_bam
  File input_bam_index
  String base_file_name
  File sorted_bed
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf
  String docker
  command {
    set -eo pipefail
    bcftools mpileup \
      --max-depth 10000 \
      --max-idepth 10000 \
      --annotate "FORMAT/AD,FORMAT/DP" \
      --fasta-ref ${ref_fasta} \
      --regions-file ${sorted_bed} \
      --ignore-RG \
      --no-BAQ \
      ${input_bam} | bcftools call -Oz -mv \
          -o ${base_file_name}.SAM.vcf.gz
    }
  runtime {
    docker: "${docker}"
    memory: "8GB"
    cpu: "4"
  }
  output {
    File output_vcf = "${base_file_name}.SAM.vcf.gz"
  }
}
# use bedtools to find basic QC data
task bedToolsQC {
  File input_bam
  File bed_file
  File genome_sort_order
  String base_file_name
  String docker
  command {
  set -e pipefail
  bedtools sort -g ${genome_sort_order} -i ${bed_file} > correctly.sorted.bed
  bedtools coverage -mean -sorted -g ${genome_sort_order} -a correctly.sorted.bed \
      -b ${input_bam} > ${base_file_name}.bedtoolsQCMean.txt
  }
  runtime {
    docker: "${docker}"
    memory: "4GB"
    cpu: "4"
  }
  output {
    File meanQC = "${base_file_name}.bedtoolsQCMean.txt"
  }
}
# align to genome
task BwaMem {
  File input_fastq
  String base_file_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Int threads
  String docker
  command {
    set -e pipefail
    bwa mem \
      -p -v 3 -t ${threads} -M \
      ${ref_fasta} ${input_fastq} | samtools view -1bS -@ "${threads}-1" > ${base_file_name}.aligned.bam
  }
  runtime {
    docker: "${docker}"
    memory: "32GB"
    cpu: "${threads}"
  }
  output {
    File output_bam = "${base_file_name}.aligned.bam"
  }
}


# get hybrid capture based QC metrics via Picard
task CollectHsMetrics {
  File input_bam
  String base_file_name
  File ref_fasta
  File ref_fasta_index
  File intervals
  String docker
  command {
    set -e pipefail

    gatk --java-options "-Xmx4g" \
      CollectHsMetrics \
      --INPUT=${input_bam} \
      --OUTPUT=${base_file_name}.picard.metrics.txt \
      --REFERENCE_SEQUENCE=${ref_fasta} \
      --ALLELE_FRACTION=0.01 \
      --BAIT_INTERVALS=${intervals} \
      --TARGET_INTERVALS=${intervals} \
      --PER_TARGET_COVERAGE=${base_file_name}.picard.pertarget.txt 
  }
  runtime {
    docker: "${docker}"
    cpu: 2
    memory: "2GB"
  }
  output {
    File picardMetrics = "${base_file_name}.picard.metrics.txt"
    File picardPerTarget = "${base_file_name}.picard.pertarget.txt"
  }
}

task consensusProcessingR {
  File GATKVars
  File SAMVars
  File MutectVars
  String base_file_name
  String githubRepoURL
  String scriptPath
  String githubTag
  String docker
  command {
    set -eo pipefail
    git clone --branch ${githubTag} ${githubRepoURL}
    Rscript ${sciptPath} ${GATKVars} ${SAMVars} ${MutectVars} ${base_file_name} 
  }
  runtime {
    cpu: 1
    memory: "2GB"
    docker: "${docker}"
  }
  output {
    File consensusTSV = "${base_file_name}.consensus.tsv"
  }

# HaplotypeCaller per-sample
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String base_file_name
  File intervals
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf
  String docker
  command {
    set -e pipefail

    gatk --java-options "-Xmx4g" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${base_file_name}.GATK.vcf.gz \
      --intervals ${intervals} \
      --interval-padding 100 
  }
  runtime {
    docker: "${docker}"
    memory: "8GB"
    cpu: "4"
  }
  output {
    File output_vcf = "${base_file_name}.GATK.vcf.gz"
    File output_vcf_index = "${base_file_name}.GATK.vcf.gz.tbi"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  File aligned_bam
  String base_file_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String docker
  command {
    set -eo pipefail
    
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      MergeBamAlignment \
     --ALIGNED_BAM ${aligned_bam} \
     --UNMAPPED_BAM ${unmapped_bam} \
     --OUTPUT ${base_file_name}.merged.bam \
     --REFERENCE_SEQUENCE ${ref_fasta} \
     --PAIRED_RUN true \
     --SORT_ORDER coordinate \
     --CREATE_INDEX true \
     --CLIP_ADAPTERS true \
     --MAX_RECORDS_IN_RAM 2000000
  }
  runtime {
    docker: "${docker}"
    memory: "8GB"
    cpu: "4"
  }
  output {
    File output_bam = "${base_file_name}.merged.bam"
    File output_bai = "${base_file_name}.merged.bai"
  }
}
# Mutect 2 calling
task Mutect2TumorOnly {
  File input_bam
  File input_bam_index
  String base_file_name
  File intervals
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File genomeReference
  File genomeReferenceIndex
  String docker

  command {
    set -eo pipefail

    gatk --java-options "-Xmx10g" \
      Mutect2 \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O preliminary.vcf.gz \
        --intervals ${intervals} \
        --interval-padding 100 \
        --germline-resource ${genomeReference} 

    gatk --java-options "-Xmx10g" \
      FilterMutectCalls \
        -V preliminary.vcf.gz \
        -O ${base_file_name}.mutect2.vcf.gz \
        -R ${ref_fasta} 
    }
  runtime {
    docker: "${docker}"
    memory: "12 GB"
    cpu: "4"
  }
  output {
    File output_vcf = "${base_file_name}.mutect2.vcf.gz"
    File output_vcf_index = "${base_file_name}.mutect2.vcf.gz.tbi"
  }
}
# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  File input_bam
  String base_file_name
  String docker
  command {
    set -e pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      SamToFastq \
			--INPUT=${input_bam} \
			--FASTQ=${base_file_name}.fastq \
			--INTERLEAVE=true \
			--INCLUDE_NON_PF_READS=true 
  }
  runtime {
    docker: "${docker}"
    memory: "6GB"
    cpu: "2"
  }
  output {
    File output_fastq = "${base_file_name}.fastq"
  }
}

# Prepare bed file and check sorting
task SortBed {
  File unsorted_bed
  File ref_dict
  String docker
  command {
    set -e pipefail

    echo "Sort bed file"
    sort -k1,1V -k2,2n -k3,3n ${unsorted_bed} > sorted.bed

    echo "Transform bed file to intervals list with Picard----------------------------------------"
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      BedToIntervalList \
      -I=sorted.bed \
      -O=sorted.interval_list \
      -SD=${ref_dict}
  }
  runtime {
    docker: "${docker}"
    memory: "4GB"
    cpu: "2"
  }
  output {
    File intervals = "sorted.interval_list"
    File sorted_bed = "sorted.bed"
  }
}






