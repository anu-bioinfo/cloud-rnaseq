cwlVersion: v1.0
class: CommandLineTool
baseCommand: python
stdout: log.txt
stderr: error.txt
inputs: 
  script: 
    type: File
    inputBinding: 
      position: 1
    default:
      class: File
      location: bcl2fastq-cwl.py
  s3_path:
    type: string
    inputBinding:
      position: 2
      prefix: -d
  sample_sheet_s3_path: 
    type: string
    inputBinding:
      position: 3
      prefix: -s
outputs:
  output:
    type: stdout
  error:
    type: stderr
