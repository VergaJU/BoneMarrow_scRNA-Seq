#! /bin/bash


source /usr/local/bin/functions.sh # source the files with functions

# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
exec > >(tee -ia main.LOG)

# source conda commands
source /opt/conda/etc/profile.d/conda.sh


# path for the reference index to obtain the counts (for kb)
ref="/lib/ref"

# reference index (for kb)
index=$ref/index.idx

# reference transcripts to gene (for kb)
t2g=$ref/t2g.txt

# first argument for the csv files containing the samples, format:
# patient_ID,sample_ID,entry(url or sra entry),tech(10xv2 or 10xv3)
csv=${1}

# run the pipeline
while IFS=, read -r patient name entry tech; do
    cd /var/
    check_filetype ${entry} # check if url for bam or sra entry
    printf "Processing sample $name\n"
    create_folder patient ${patient}
    create_folder sample ${name}
    if [[ "$(find .)" =~ "counts_filtered" ]] # Check if final output already exist and skip 1 loop
    then
        printf "Sample already pre-processed!\n"
        cd ../../
    elif [[ "$(find .)" =~ "fastq" ]] # check if fastq files are already present
    then
        if [[ "${filetype}" == "BAM" ]]
        then
            check_bam # if entry is bam file check if converted correctly
        elif [[ "${filetype}" == "SRA" ]]
        then
            check_sra # if entry is sra file check if converted correctly
        fi

        #determine_tech
        get_counts $tech # get counts using kb
        cd ../../
    else
        if  [[ "${filetype}" == "BAM" ]] # check if entry is bam or sra
        then
            download_bam ${entry} # download the bam file (10x genomics)
            cellranger ${name} # convert file in fastq files using cellranger
            check_bam # check if converted correctly
        else
            download_sra ${entry} # if entry is sra, automatically download and convert file in fastq using fasterq-dump
            check_sra # check if converted correctly
        fi

        #determine_tech
        get_counts $tech # get counts using kb
        cd ../../
    fi
done < "$csv" 