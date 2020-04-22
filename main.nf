import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process preflight {
    container "nkrumm/nextflow-demux:latest"
    input:
        file(samplesheet) from Channel.fromPath(params.samplesheet)
    output:
        file("${params.run_id}.samplesheet.csv") into samplesheet_ch
        file("${params.run_id}.config.json") into config_file_ch
        
    // this is high so that AWS provisions a big server immediately
    cpus 30
    memory '68 GB'

    script:
        umi_options = params.basemask != "" ? "--is-umi --basemask ${params.basemask}" : ""
        fwd_adapter = params.fwd_adapter ? "--fwd-adapter ${params.fwd_adapter}" : ""
        rev_adapter = params.rev_adapter ? "--rev-adapter ${params.rev_adapter}" : ""
        merge_lanes = params.merge_lanes ? "--merge-lanes" : ""
        """
        parse_samplesheet.py \
            --input ${samplesheet} \
            --output ${params.run_id} \
            --project-name ${params.project_name} \
            --library-type ${params.library_type} \
            ${merge_lanes} \
            ${umi_options} \
            ${fwd_adapter} ${rev_adapter}
        """
}

def jsonSlurper = new JsonSlurper()
demux_config = jsonSlurper.parseText(config_file_ch.first().text.getVal())

process demux {
    echo true
    cpus 30
    memory '68 GB'
    container "nkrumm/nextflow-demux:latest"
    publishDir params.output_path, pattern: 'output/Reports', mode: 'copy', saveAs: {f -> f.replaceFirst("output/", "")}, overwrite: true
    publishDir params.output_path, pattern: 'output/Stats', mode: 'copy', saveAs: {f -> f.replaceFirst("output/", "")}, overwrite: true
    publishDir params.output_path, pattern: 'output/*.fastq.gz', mode: 'copy', overwrite: true // these are the "Undetermined" fastq.gz files

    input:
        file(samplesheet) from samplesheet_ch
    output:
        file("output/**.fastq.gz") into demux_fastq_out_ch
        file("output/Reports")
        file("output/Stats")
        file("output/Stats/Stats.json") into stats_json_multiqc
        path("inputs/${params.run_id}/InterOp/*") into interop_input
        file("inputs/${params.run_id}/RunInfo.xml") into interop_input_xml

    script:
        rundir = "inputs/${params.run_id}"
        basemask = demux_config.basemask ? "--use-bases-mask " + demux_config.basemask : ""
        merge_lanes = params.merge_lanes ? "--no-lane-splitting" : ""
        """
        mkdir -p ${rundir}
        aws s3 sync --only-show-errors ${params.run_folder} ${rundir}

        if [ -f ${rundir}/Data.tar ]; then
         tar xf ${rundir}/Data.tar -C ${rundir}/
         rm ${rundir}/Data.tar
        fi
        
        bcl2fastq \
            --runfolder-dir ${rundir} \
            --output-dir output/ \
            --sample-sheet ${samplesheet} \
            --min-log-level WARNING \
            --ignore-missing-positions \
            --ignore-missing-bcls \
            --ignore-missing-filter \
            --barcode-mismatches 0,0 \
            --create-fastq-for-index-reads \
            ${basemask} \
            ${merge_lanes} \
            --mask-short-adapter-reads 0
        """
}

demux_fastq_out_ch.flatMap()
    // process output directories from bcl2fastq and re-associate them with config entries
    // we need a key of (lane, project_name, sample_id) to map back to the nested config file.
    .view()
    .filter { path -> !("${path.getName()}" =~ /^Undetermined_S0_/) } // first filter ignore Undetermined read files
    .filter { path -> !("${path.getName()}" =~ /I\d_001.fastq.gz$/) } // filter out indexing reads
    .view()
    .map { path -> 
          def (filename, project_name, rest) = path.toString().tokenize('/').reverse() // tokenize path
          if (params.merge_lanes){
            // if bcl2fastq was run with --no-lane-splitting, no lane token will be availabe in filename
            // e.g., "294R09-A02-MONCv1-NA12878_S1_I1_001.fastq.gz"
            (extension, read_num, sample_num, sample_id) = filename.tokenize("_").reverse()
            lane = "all"
          } else {
            // otherwise, if lanes are split token will be availabe in filename
            // e.g., "294R09-A02-MONCv1-NA12878_S1_L001_I1_001.fastq.gz"
            (extension, read_num, lane, sample_num, sample_id) = filename.tokenize("_").reverse()
            lane = lane.replaceAll("L00", "") // just use lane number
          }
          def key = tuple(lane, project_name, sample_id)
          return tuple(key, path)
    }
    .view()
    .groupTuple() // group FASTQ files by key
    .map { key, files -> 
          // attach config information
          def c = demux_config.lanes[key[0]][key[1]][key[2]]
          if (c.is_umi){
            // if umi, read 2 is the UMI read, and read 3 is the reverse read, so re-order appropriately
            [key, [files[0], files[2], files[1]], c] 
          } else {
            // for non-umi read 1 and read 2 are ordered sequentially
            [key, [files[0], files[1]], c] 
          }
         } 
    .view{ JsonOutput.prettyPrint(JsonOutput.toJson(it[2])) } // diagnostic print'
    .branch { 
        // send files to postprocess_ch.umi_true if is_umi is set.
        umi_true: it[2].is_umi
        umi_false: it[2].is_umi==false
    }
    .set{postprocess_ch}


process postprocess_umi {
    container "nkrumm/nextflow-demux:latest"
    echo true
    cpus 4
    memory '7 GB'
    input:
        set key, file(fastqs), config from postprocess_ch.umi_true
    output:
        set key, file("processed/*.fastq.gz"), config into umi_out_ch
    script:
        fastq1 = fastqs[0]
        fastq2 = fastqs[1]
        umi = fastqs[2]
        """
        mkdir -p processed/
        paste \
            <(zcat ${fastq1}) \
            <(zcat ${umi}) \
            | awk 'NR%4==1{readname=\$1}
                   NR%4==2{seq=\$1; umi=\$2}
                   NR%4==0 {print readname " RX:Z:"umi"\\n"seq"\\n+\\n"\$1;}' \
            | gzip > processed/1.fastq.gz ;

        paste \
            <(zcat ${fastq2}) \
            <(zcat ${umi}) \
            | awk 'NR%4==1{readname=\$1}
                   NR%4==2{seq=\$1; umi=\$2}
                   NR%4==0 {print readname " RX:Z:"umi"\\n"seq"\\n+\\n"\$1;}' \
            | gzip > processed/2.fastq.gz ;
        """
}


// merge umi_out_ch with postprocess_ch.umi_false
umi_out_ch.mix(postprocess_ch.umi_false)
    // route samples to trim process if needed
    .branch {
        trim_true: it[2].fwd_adapter || it[2].rev_adapter
        trim_false: true // default, if no fwd or rev adapters are specified
    }
    .set { trim_in_ch }


process trim {
    cpus 4
    memory '7 GB'
    container 'nkrumm/atropos-paper:latest'
    input:
        set key, file(fastqs), config from trim_in_ch.trim_true
        //file(params.adapters)
    output:
        set key, file("trimmed/*.fastq.gz"), config into trim_out_ch
    
    script:
        // TODO: consider using --stats and --report-file here to get report of read trimming
        def trim_options = config.additional_trim_options ? config.additional_trim_options : ""
        if (fastqs.size() == 2)
            """
            mkdir -p trimmed/
            atropos \
              -T ${task.cpus} \
              -a ${config.fwd_adapter} \
              -A ${config.rev_adapter} \
              ${trim_options} \
              -pe1 ${fastqs[0]} \
              -pe2 ${fastqs[1]} \
              -o trimmed/1.fastq.gz \
              -p trimmed/2.fastq.gz \
              > report.txt
            """
        else
            """
            mkdir -p trimmed/
            atropos \
              -T ${task.cpus} \
              -a ${config.fwd_adapter} \
              ${trim_options} \
              -se ${fastqs[0]} \
              -o trimmed/1.fastq.gz \
              > report.txt
            """
}

// merge together trimmed and non-trimmed readgroups
// if params.qc_merge_lanes is set, we will output merged lanes from FASTQC
// note that we *never* combine lanes in the library (fastq.gz) output.

trim_out_ch.mix(trim_in_ch.trim_false)
    .into { qc_in_ch; finalize_libraries_in_ch }

if (params.merge_lanes || params.qc_merge_lanes){
    // merge together lanes in QC step
    // group by sample project + sample ID (omits lane from the key)
    qc_in_ch.map{ key, files, config -> [["all", key[1], key[2]], files, config]}
            .groupTuple(size: demux_config.number_of_lanes, remainder: true)
            .map{ key, files, config -> [key, files.flatten(), config[0]] }
            .set{ fastqc_in_ch }
} else {
    // consider lanes separately when running FASTQC
    qc_in_ch.set{ fastqc_in_ch }
}

process fastqc {
    cpus 2
    memory '4 GB'
    container 'quay.io/biocontainers/fastqc:0.11.8--1'

    publishDir params.output_path, pattern: "*.html", mode: "copy", overwrite: true
    
    input:
        set key, path("fastq??.fq.gz"), config from fastqc_in_ch
    output:
        path "fastqc/*", type:"dir" into fastqc_report_ch
   
    script:
        lane = key[0] // could be a number or "all" if qc_merge_lanes is true (see above) 
        readgroup = "${params.fcid}.${lane}.${config.index}-${config.index2}"
        sample_name = config.Sample_Name //"${config.Sample_Name}:${config.library_type}:${readgroup}:${lane}"
        fastqc_path = "fastqc/${sample_name}/"
        
        """
        mkdir -p ${fastqc_path}
        zcat fastq*.fq.gz | fastqc --quiet -o ${fastqc_path} stdin:${sample_name}
        """
}   

process finalize_libraries {
    container "nkrumm/nextflow-demux:latest"
    cpus 2
    memory '4 GB'
    // consider moving the code to put final libraries here
    // would handle external s3 destinations (+ multiple destinations)
    // would not confound fastqc step above
    // would be an extra copy/move of the data (as an extra process)
    // would handle final s3 tagging to specify lifecyle policies
    publishDir params.output_path, pattern: "**.fastq.gz", mode: "copy", overwrite: true
    input:
        set key, file(fastqs), config from finalize_libraries_in_ch
    output:
        path("libraries/**.fastq.gz")
    script: 
        lane = key[0] // note this is either the lane number or "all" if params.merge_lanes == true
        readgroup = "${params.fcid}.${lane}.${config.index}-${config.index2}"
        library_path = "libraries/${config.Sample_Name}/${config.library_type}/${readgroup}"
        if (fastqs.size() == 2)
            """
            mkdir -p ${library_path}
            mv ${fastqs[0]} ${library_path}/1.fastq.gz
            mv ${fastqs[1]} ${library_path}/2.fastq.gz
            """
        else
            """
            mkdir -p ${library_path}
            mv ${fastqs[0]} ${library_path}/1.fastq.gz
            """
}


process interop {
    container 'quay.io/biocontainers/illumina-interop:1.1.8--hfc679d8_0'
    memory '4 GB'
    cpus 1
    input:
        path("${params.run_id}/InterOp/*") from interop_input
        file("${params.run_id}/RunInfo.xml") from interop_input_xml
    output:
        path("*.csv") into interop_output

    script:
        """
        interop_summary --csv=1 ${params.run_id}/ > interop_summary.csv
        interop_index-summary --csv=1 ${params.run_id}/ > interop_index-summary.csv
        """
}

process multiqc {
    echo true
    container 'quay.io/biocontainers/multiqc:1.7--py_4'
    memory '4 GB'
    cpus 4
    input:
        path('*') from fastqc_report_ch.flatMap().collect()
        //file("Stats.json") from stats_json_multiqc
        path("interop/*") from interop_output
        //file('*') from trim_reads_report_ch.collect()
    output:
       file "multiqc_report.${params.fcid}.html"
    publishDir params.output_path, saveAs: {f -> "multiqc/${f}"}, mode: "copy", overwrite: true
    script:
        """
        multiqc -v --filename "multiqc_report.${params.fcid}.html" .
        """
}