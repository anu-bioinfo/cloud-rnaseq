task s3cp {
  String source
  String destination
  String flags

  command {
    aws s3 cp ${source} ${destination} ${flags}
  }

  output {
      String destination_path = destination
  }
}

task bcl2fastq {
    String bcl_path
    String fastq_path
    String sample_sheet
    command {
        mkdir -p ${fastq_path}; bcl2fastq --sample-sheet ${sample_sheet} -R ${bcl_path} -o ${fastq_path}
    }
    output {
       String output_fastq_path = fastq_path
    }
}

task find_wc_path {
    String wildcard_string
    command {
        ls -d ${wildcard_string}
    }
    output {
        String output_dir = read_string(stdout())
    }
}

task filter_for_fastq_dirs {
    String output_dir 
    command {
      ls ${output_dir} |grep -v '.gz' |grep -v Reports|grep -v Stats

    }
    output {
        Array[String] fastq_dirs = read_lines(stdout())
    }
}

workflow wf {
  String s3_dir
  String experiment_id
  String sample_sheet_file_name 

  call s3cp {
      input: 
          source=s3_dir+"/"+experiment_id+"/bcl/", 
          destination="/mnt/data/" + experiment_id + "/bcl",
          flags=" --recursive"
  }
  call s3cp as s3cpconfig {
      input: 
          source=s3_dir+"/"+experiment_id+"/config/"+sample_sheet_file_name, 
          destination="/mnt/data/" + experiment_id + "/"+sample_sheet_file_name,
          flags=""
  }
  call bcl2fastq {
      input:
          bcl_path=s3cp.destination_path,
          fastq_path="/mnt/data/" + experiment_id + "/fastqs", 
          sample_sheet = s3cpconfig.destination_path
  }
  call find_wc_path {
      input: wildcard_string = bcl2fastq.output_fastq_path + "/Reports/html/*/all/all/all/"
  }

  call s3cp as s3cpreports {
      input: 
          source=find_wc_path.output_dir,
          destination=s3_dir + "/" + experiment_id + "/reports/bcl2fastq/",
          flags=" --recursive "
  }

  call filter_for_fastq_dirs {
      input: output_dir=bcl2fastq.output_fastq_path
  } 

  scatter (fastq_dir in filter_for_fastq_dirs.fastq_dirs) {
     call s3cp as s3cpdata {
         input: 
             source=bcl2fastq.output_fastq_path+"/"+fastq_dir,
             destination=s3_dir + "/" + experiment_id + "/rawdata/",
             flags=" --recursive "
     }
  }
}
