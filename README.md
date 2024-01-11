# split_spliced_reads
Split spliced reads in a similar way compared to the GATK tool  
  
Usage:  
  
split_spliced_reads [options]  
  
-i --input-bam-file <text>     Single input bam file (required)  
-o --output-bam-file <text>    Single output bam file (optional default _split.bam)  
-s --max-alignment-score <int> Maximum alignment score (optional default 255)  
-c --hard-clipping <void>      When introduced use hard clipping  
-h --help <void>               This help  
  
-s sets a maximum alignment score on the reads  
-c introduces hard clipping  
