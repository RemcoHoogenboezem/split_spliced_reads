# split_spliced_reads
Split spliced reads in a similar way compared to the GATK tool\n
\n
Usage:\n
\n
split_spliced_reads [options]\n
\n
-i --input-bam-file <text>     Single input bam file (required)\n
-o --output-bam-file <text>    Single output bam file (optional default _split.bam)\n
-s --max-alignment-score <int> Maximum alignment score (optional default 255)\n
-c --hard-clipping <void>      When introduced use hard clipping\n
-h --help <void>               This help\n
\n
-s sets a maximum alignment score on the reads\n
-c introduces hard clipping \n
