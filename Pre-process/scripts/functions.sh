#! /usr/bin/bash



# Check if the folder already exist, if not create it
# First argument: str (patient or sample)
# Second argument: variable with patient or sample ID
# e.g. create_folder "patient" 51449
# e.g. create_folder "sample" SAMN17249432

create_folder () {
    folders=($(ls))
    if [[ "${folders[@]}" =~ "${2}" ]]
    then
        printf "Entering the ${1} ${2} folder\n"
        cd ${2}
    else
        printf "Creating and entering the ${1} ${2} folder\n"
        mkdir ${2}
        cd ${2}
    fi
}


# check if the csv file contains links to download bam files or SRA entries
check_filetype () {
    if [[ "${1}" =~ "https" ]]
    then
        filetype="BAM"
    elif [[ "${1}" =~ "SRR" ]]
    then
        filetype="SRA"
    fi
}


# DOWNLOAD THE BAM FILE
# Argument is the url (${url})

download_bam () {
    printf "Downloading bam file...\n"
    wget -c ${1}
    filename=$(ls) # update filename
}

# Download the SRA file and gunzip it in order to be correctly processed in next steps
download_sra () {
    printf "Downloading SRA file..\n";
    fasterq-dump SRR12506863 -o SRR12506863 -O fastq/ -S --include-technical |& tee -a fasterq-dump.LOG
    cache-mgr -c ~/.sratoolkit.2.11.3-ubuntu64/cache/

    #gzip ./fastq/${1}*
    for file in ./fastq/* ; do mv $file ${file//_/_R} ; done

}

# Run cellranger to convert the bam file to fastq files, as input it takes the filename 
# First argument: sample name (${name})
# Then control if the conversion was done correctly, rename the fastq files with the sample ID
# Move the files in the fastqs/ folder
# Remove the other folders and bam files

cellranger () {
    printf "running CellRanger bamtofastq to obtain the fastq files\n"
    ## Use cellranger to get the fastq files
    ~/.cellranger-6.1.2/lib/bin/bamtofastq \
        --nthreads=6 \
        --reads-per-fastq=100000000000 \
        --traceback \
        ${filename} \
        ./fastq/ \
        |& tee -a bamtofastq.LOG
}

check_bam () {
    # Control the log file and reformat the output for forwarding analyses
    if [[ "$(tail -n1 bamtofastq.LOG)" =~ "Writing finished" ]]
    then
        printf "bamtofastq conversion finished correctly \nRemoving bam files and moving fastq files\n"
        fastqfiles=($(find ./fastq/ | grep bamto)) # store fastqfiles names and path
        for file in ${fastqfiles[@]};do
            suffix=$(echo ${file} | awk -F "/" '{print  $NF}' | cut -d "_" -f 2-) # save the suffix name
            mv ${file} ./fastq/$name"_"$suffix # rename the files with the name given and move in fastq directory
        done
        rm -fr ./fastq/*/* # remove empty folders
        rm -f *.bam* # remove bam file
    else
        printf "ERROR: bamtofastq conversion failed, EXITING.\n"
        exit 1
    fi
}

check_sra () {
    # control fasterq-dump output and reformat the output for fowrarding analyses
    if [[ "$(tail -n1 fasterq-dump.LOG)" =~ "reads written" ]]
    then
        printf "SRA file converted correctly\n"
    else
        print "ERROR: File not downloaded correctly, EXITING\n"
        exit 1
    fi
}

# Determine the technology used to exert the scRNA-Seq experiment. Now it recognize only 10x genomics experiments
# No arguments needed
determine_tech () {
    printf "Checking technology\n"

    fastqs=($(find . | grep _R[1-2] | sort))

    echo ${fastqs[@]}

    # detect chemistry version by checking the length of the first 20 R1s and taking the floor of the average
    r1_len="$(zcat fastq/*R1* | head -80 | sed -n 2~4p | awk '{ print length }' | awk -F : '{sum+=$1} END {print sum/NR}' | cut -f1 -d".")"
    if [[ $r1_len =~ 28 ]]
    then
        tech="10xv3"
    elif [[ $r1_len =~ 26 ]]
    then
        tech="10xv2"
    else
    echo 'Unable to determine version\n'
    fi
}


# Get the counts using the kb wrapper for kallisto|bustools.
# Firstly activate the conda evironment, the environment need to bre created using the pre_process.yml file 
get_counts () {
    #activate conda env
    conda activate pre_process

    fastqs=($(find . | grep _R[1-2] | sort))

    # run kallisto|bustools
    printf "running kb to obtain counts\n"
    kb count \
        --verbose \
        -t 6 \
        --cellranger \
        -i $index \
        -g $t2g \
        -x 10xv3 \
        -o ./kb_out \
        ${fastqs[@]} \
        --overwrite \
        |& tee -a kallisto.LOG 

    if [[ "$(find kb_out) " =~ "cellranger" ]]
    then
        # filter empty droplets
        printf "filtering empty droplets with FDR threshold of 0.1\n"
        ~/Documents/PROJECT/scripts_env/scripts/filter_empty_v2.R ./kb_out/counts_unfiltered/cellranger/ 0.1 \
        |& tee -a filter_empty.LOG
    else
        printf "ERROR: kallisto|bustool failed obtaining the counts. EXITING"
        exit 1
    fi

    if [[ "$(find kb_out) " =~ "counts_filtered" ]]
    then
        # remove fastq and bus files
        printf "Files preprocessed correctly\nRemove fastq and bus files to save hard drive space\nReturn in root folder\n"
        find . -name "*.fastq" -delete
        find . -name "*.bus" -delete
    else
        printf "ERROR: Cells not correctly filtered. EXITING"
        exit 1
    fi
}