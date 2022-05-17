#! /usr/bin/bash


source ~/Documents/PROJECT/scripts_env/scripts/functions.sh # source the files with functions

# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
exec > >(tee -ia main.LOG)

# source conda commands
source /home/jacopo/miniconda3/etc/profile.d/conda.sh


ref="/home/jacopo/Documents/PROJECT/MM/datasets/SRA_data/ref"

index=$ref/index.idx
t2g=$ref/t2g.txt


csv=${1}
while IFS=, read -r patient name entry tech; do
    check_filetype ${entry}
    printf "Processing sample $name\n"
    create_folder patient ${patient}
    create_folder sample ${name}
    if [[ "$(find .)" =~ "counts_filtered" ]] # Check if final output already exist and skip 1 loop
    then
        printf "Sample already pre-processed!\n"
        cd ../../
    elif [[ "$(find .)" =~ "fastq" ]]
    then
        if [[ "${filetype}" == "BAM" ]]
        then
            check_bam
        elif [[ "${filetype}" == "SRA" ]]
        then
            check_sra
        fi

        #determine_tech
        get_counts $tech
        cd ../../
    else
        if  [[ "${filename}" == "BAM" ]]
        then
            download_bam ${entry}
            cellranger ${name}
            check_bam
        else
            download_sra ${entry}
            check_sra
        fi

        #determine_tech
        get_counts $tech
        cd ../../
    fi
done < "$csv" 