#### MIUMS version:V1.0 
MIcrobial Identification Using Marker Sequence (MIUMS): A state-of-the-art taxonomy prediction tool for (meta)genomic sequences that can identify novel sequences for which there are no closely related sequences available in NCBI RefSeq database.

#### MIUMS Description:
Identification and classification of metagenomic sequences by sequence homology search against publicly available sequence databases (e.g. GenBank, RefSeq etc.) is a common practice. Advances in NGS techniques complemented with growing popularity of metagenomics and lower sequencing cost has been constantly populating these databases with new sequences. While the abundance of new sequences has benefited the identification/classification of metagenomic sequences, it has also significantly increased the amount of computation and data storage requirements. Although it is possible to construct the reference database from an entire database such as RefSeq, the storage and memory requirements to run the prediction becomes impossible in a non-server-like computing environment with limited resources. To address these issues, MIUMS (version 1.0) utilizes a small database of protein sequence fragments that are highly specific to their source organism, The protein sequence fragments (also referred as the protein marker sequences in MIUMS) are collected from full length protein sequences of 224 archaeal, 2810 bacterial and 3958 viral species (minimum sequence length of 5000 nucleotides; published before 1st of March 2017; Refseq release version 79). In addition, an eukaryotic protein sequence database was constructed from 16179736 eukaryotic proteins to be used as negative control. 

#### MIUMS features:
MIUMS (V1.0) is a complete package comes in a docker container. It supports:

* Downloading raw reads from SRA and ENA databases just specifying relative accessions of the sequence library
* Supports both fasta raw reads or FASTQ raw reads (single, paired, or interleaved) files as input
* Read error coorections (using BBmap suite of tools) 
* Assembly of raw reads (using MegaHit or Spades assemblers) or provide your own assembled contigs file
* Taxonomic classification (archaea, bacteria and viruses) of reads and contigs
* Binning of classified sequences at a given taxa level
* Can identify novel sequences for which there are no closely related sequences available in RefSeq database
* It also offers deep learning techniques to identify more closely related sequences in iterative manner (recommended for metaviromes)
* Generate output in formats similar to other widely used tools (e.g. kraken, kaiju, CLARK)   

#### MIUMS output:
By default MIUMS generates a main taxonomy output table and a secondary taxonomy prediction table. The main output taxa table lists inputted sequence IDs followed by predicted taxon ID (determined by LCA), organism (e.g. archaea, bacteria or viruses), genetic compartment (e.g. chromosome, plasmid, Bacteriophage, prophage etc.) and six level of taxonomy (i.e. superkingdom, phylum, class, order, family, genus and species) in tabular format. The second output table shows sequences that are potentially archaea, bacteria or viruses but contains one or more conflicting (i.e. non-superkingdom specific) reference marker protein fragments in their protein sequences. MIUMS also generates a summary file from the main output table, which shows the total number of predicted archaea, bacteria, viruses as well as the total number of bacteriophage, prophage and plasmids. It also gives the error corrected reads fasta file, assembled contigs and a file containing predicted proteins from the contigs.

#### MIUMS installation & displaying help:

     docker run --rm -it --name miums ambarishbiswas/miums:v1 -h

#### MIUMS examples:

     docker run -v $(pwd):/tmp --rm -it ambarishbiswas/miums:v1 -test fast

     docker run -v $(pwd):/tmp --rm -it ambarishbiswas/miums:v1 -test full

     docker run -v $(pwd):/tmp --rm -it ambarishbiswas/miums:v1 -sample_id ABCD -f test.fastq -o EFGH

     docker run -v /tmp:/tmp --rm -it ambarishbiswas/miums:v1 -sample_id ABCD -1 test_R1.fastq -2 test_R2.fastq -o EFGH

     docker run -v /tmp:/tmp --rm -it ambarishbiswas/miums:v1 -sample_id SRR390728 -sra SRR390728 -o SRR390728_output

     docker run -v /tmp/ABCD:/tmp --rm -it ambarishbiswas/miums:v1 -sample_id XYZ -ena ERR2105523 -o output_of_ERR2105523


#### MIUMS commandline parmeters:

<h5>Options must be provided:</h5>
  
	-sample_id                            XYZ1                      [Every MIUMS run must be supplied with a sample_id.]  
 	-o/-out_dir                           a_folder_name             [A folder with the provided name will be created in the current directory; No '/' or '\' is accepted]

<h5>Options for input sequence(s) in FASTA format:</h5>
  
 	-i/-f                                 input_fasta_sequence_file [either gzip compressed (i.e. with extension .gz) or uncompressed FASTA sequence file. Supported extensions are .fa, .fna, .fasta ]



<h5>Options for input sequence(s) in FASTQ format:</h5>
  
 	-s/-r                                 input_fastq_sequence_file [either gzip compressed (i.e. with extension .gz) or uncompressed FASTQ file. Supported extensions are .fq, .fastq ]
 	-1                                    forward_fastq_file        [Specify the forward reads FASTQ file of paired-end libraries. Supported extensions are .fq, .fastq with/without .gz]
 	-2                                    reverse_fastq_file        [Specify the reverse reads FASTQ file of paired-end libraries. Supported extensions are .fq, .fastq with/without .gz]
 	--12                                  interleaved_fastq_file    [Specify a FASTQ file with forward and reverse reads interleaved. Supported extensions are .fq, .fastq with/without .gz]


<h5>Options for publicly available (meta)genomic NGS reads:</h5>
  
 	-sra                                  SRA_Accession             [Accession of (meta)genomic NGS reads from NCBI SRA database]
 	-ena                                  ENA_Accession             [Accession of (meta)genomic NGS reads from ENA database]


<h5>Options for reads specific operations:</h5>
  
 	-read_error_correction                0/1                       [Default is set to 1; MIUMS uses BBmap suite of tools for error correction] 
 	-min_read_length                      20                        [Default value set to 20]
 	-subsample_reads_to                   N                         [A positive integer 'N' can be provided to create a subsample_reads file using randomly selected reads; This is useful where multiple metagenomic libraries with different reads depth are being analyzed; ]

 	-subsample_reads_min_length           N                         [Default value set to 75; Any positive integer is supported; Remaining reads after read_error_correction Reads will be checked for length >= N before including them in the subsample_reads file ]



<h5>Options for reads assembly:</h5>
  
 	-use_reference_sequence_from_assembly 0/1                       [Default is set to 1, which means assembly will be done on the inputted reads/sequences in the inputted prior constructing a sample specific reference library of marker sequence, which will then be used for classifying archaea, bacteria and viral reads ]

 	-assembler                            Spades/Megahit            [Default set to SPAdes; to use Megahit assembler specify '-assembler megahit']
 	-assembled_contigs_file               FASTA_seq_file            [specify an existing (multi)FASTA sequence file. Taxonomy prediction of the contigs will be carried out prior constructing a sample specific reference library of marker sequence. Either gzip compressed (i.e. with extension .gz) or uncompressed FASTA sequence file is expected ]



<h5>Options for Taxonomy search:</h5>
  
 	-taxa_assignment_of_contigs           0/1                       [Default is set to 1, i.e. the taxonomy will be predicted on assembled contigs;]
 	-taxa_assignment_of_reads             0/1                       [Default is set to 1, i.e. the taxonomy will be predicted on the inputted reads (or reads from the subsampled reads file if opted); ]



<h5>Options for Iterative classification:</h5>
  
 	-iterative_search                     0/1                       [Default is set to 1, i.e. the taxonomy will be predicted on the contigs in an iterative way, where in each run 'Reference Library of Marker Sequence' created in the previous run(s) will be used for classifying more 'related' contigs ]

 	-is_sensitivity                       HIGH/MEDIUM/LOW           [Default is set to HIGH; The options (i.e. high, medium or low) controls the level of effort given while constructing the 'Reference Library of Marker Sequence (RLMS)'; LOW= only proteins classified in the previous round; MEDIUM= (proteins classified in the previous round + MIUMS RLMS); HIGH= (proteins classified in the previous round + MIUMS RLMS + Negative_control); ]

 	-is_min_id                            35                        [The minimum identitiy on the overlapping regions between reference and query sequences; Default is set to 35, but any positive integer is supported; ]

 	-is_min_bitscore                      51                        [The bitscore cutoff of the overlapping regions between reference and query sequences; Default is set to 51, but any positive integer is supported; ]

 	-is_min_coverage                      51                        [The minimum coverage cutoff of the overlapping regions between reference and query sequences; Default is set to 51, but any positive integer is supported; ]

 	-is_increment                         5                         [Default is set to 5, i.e. in each iteration the is_min_id, is_min_bitscore and is_min_coverage will be incremented by 5 (until they reach 90); this option ensures that only the more and more sequences with homologous proteins are classified; any positive/negative real number can be used; ]



<h5>Options for additional reporting:</h5>
  
 	-bin_sequences_at                     X                         [Supported values for X are: superkingdom, phylum, class, order, family, genus or species; By default sequence binning will not be done ]

 	-report_dual_predictions              0/1                       [Default is set to 1, MIUMS can identify sequences that contains proteins from organisms belongs to inter-superkingdom; typical example is viral proteins in bacterial sequence ]

 	-kraken_style_output                  0/1                       [Default is set to 0; Use '1' to create a output file similar to what kraken, CLARKS and kaiju provides]


<h5>Other options:</h5>
  
 	-T/-threads                           N                         [By default the program uses 4 threads; Any positive integer is accepted]
 	-q/-quiet                             0/1                       [Default is set to 1, which shows program step-by-step logs; Use '0' to turn of the logging; Note, a file log.log will still be created in the output_folder ]

 	-h/--h/-help/--help                                             [Shows this help]
 	-v/-version                                                     [Shows MIUMS Version ]
 	-rm                                                             [removes the output folder specified by -o or -out_dir (if exists);]
 	-clean                                                          [removes all intermediate temporary files once a process is finished;]
 	-test                                                           [supported options are: full, fast or custom; If -test is provided, MIUMS uses pre-supplied fastq files for testing the pipeline; the option 'full' sets parameters for most-comprehensive MIUMS run; The option 'fast' does the least-comprehensive MIUMS run where reads are classified directly using MIUMS RLMS; The 'custom' option takes all other user parameters except input fasta/fastq related parameters; ]


#### General information:

  For all queries, bugs and updates (and even suggestions how to make MIUMS better) please email Ambarish Biswas [ambarishbiswas@gmail.com]
  
  Latest version and other informationcan can be found at: https://github.com/ambarishbiswas/miums_v1.0

  The MIUMS V1.0 is developed as a proof-of-concept for the upcoming/associated paper, and is free for academic use. If you use MIUMS for your reasearch, please cite the associated paper.  
