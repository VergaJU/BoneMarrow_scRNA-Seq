#! /bin/bash



# Check if the folder already exist, if not create it
# First argument: str (patient or sample)
# Second argument: variable with patient or sample ID
# e.g. create_folder "patient" 51449
# e.g. create_folder "sample" SAMN17249432

# create folder for the patient and samples (from first 2 columns of the csv input file)
## input:
# 1- string: patient or sample
# 2- string: entry (first or second column input csv)

## output
# working directory for the patient/sample
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
## input: 3rd column of the csv file, SRA entry or url to download bam file
## output: variable "filetype" for next steps
check_filetype () {
    if [[ "${1}" == https* ]]
    then
        filetype="BAM"
    elif [[ "${1}" == SRR* ]]
    then
        filetype="SRA"
    fi
}


# DOWNLOAD THE BAM FILE
## input: url (${url})
## output: bam file
download_bam () {
    printf "Downloading bam file...\n"
    wget -c ${1}
    filename=$(ls) # update filename
}

# Download and convert SRA file with fasterq-dump
## input: SRR entry
## output: 3 fastq files (R1, R2, R3)
download_sra () {
    printf "Downloading SRA file..\n";
    fasterq-dump ${1} -o ${1} -O fastq/ -S --include-technical |& tee -a fasterq-dump.LOG
    cache-mgr -c ~/.sratoolkit.2.11.3-ubuntu64/cache/

}

# Run cellranger to convert the bam file to fastq files, as input it takes the filename 
## Input: string, filename, bam file previously downloaded
## Output: 3 fastq files

cellranger_bam () {
    printf "running CellRanger bamtofastq to obtain the fastq files\n"
    ## Use cellranger to get the fastq files
    bamtofastq --reads-per-fastq=100000000000 --traceback ${filename} ./fastq/ |& tee -a bamtofastq.LOG
}

# Check if cellranger converted correctly the files
# Move the files in the fastqs/ folder
# Remove the other folders and bam files

check_bam () {
    # Control the log file and reformat the output for forwarding analyses
    if [[ "$(tail -n1 bamtofastq.LOG)" =~ "Writing finished" ]]
    then
        printf "bamtofastq conversion finished correctly \nRemoving bam files and moving fastq files\n"
        fastqfiles=($(find ./fastq/ | grep bamto)) # store fastqfiles names and path
        for file in ${fastqfiles[@]};do
            suffix=$(echo ${file} | awk -F "/" '{print  $NF}' | cut -d "_" -f 2-) # save the suffix name
            mv ${file} ./fastq/$name"_"$suffix # rename the files with the name given in the csv file 
                                               # and move in fastq directory
        done
        rm -fr ./fastq/*/* # remove empty folders
        rm -f *.bam* # remove bam file
    else
        printf "ERROR: bamtofastq conversion failed, EXITING.\n"
        exit 1
    fi
}

check_sra () {
    # control fasterq-dump output and reformat the output for forwarding analyses
    if [[ "$(tail -n1 fasterq-dump.LOG)" =~ "reads written" ]]
    then
        printf "SRA file converted correctly\n"
        # rename fastq files in R1, R2 anr I1
        for file in ./fastq/* ; do
        len="$(cat ${file} | head -80 | sed -n 2~4p | awk '{ print length }' | awk -F : '{sum+=$1} END {print sum/NR}' | cut -f1 -d".")";
            if [[ $len -lt 26 ]]; then
            mv $file ${file%_*.fastq}_S1_L001_I1_001.fastq
            elif [[ $len -eq 26 || $len -eq 28 ]]; then
            mv $file ${file%_*.fastq}_S1_L001_R1_001.fastq
            elif [[ $len -gt 26 && $len -gt 28 ]]; then
            mv $file ${file%_*.fastq}_S1_L001_R2_001.fastq           
            fi
        done

    else
        printf "ERROR: File not downloaded correctly, EXITING\n"
        exit 1
    fi
}

# Determine the technology used to exert the scRNA-Seq experiment. Now it recognize only 10x genomics experiments
# No arguments needed
# Not used now, technology inputted in the csv file
# TODO: FIX function for different technologies
determine_tech () {
    printf "Checking technology\n"

    fastqs=($(find . | grep _R[1-2] | sort)) # find fastq files of interest

    #echo ${fastqs[@]}

    # detect chemistry version by checking the length of the first 20 R1s and taking the floor of the average
    r1_len="$(cat fastq/*1* | head -80 | sed -n 2~4p | awk '{ print length }' | awk -F : '{sum+=$1} END {print sum/NR}' | cut -f1 -d".")"
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
# inputs: fastq files, index, t2g, technology
# output: count matrices, unfiltered/filtered
get_counts () {
    #activate conda env
    conda activate pre_process

    fastqs=($(find . | grep _R[1-2] | sort))  # find fastq files of interest

    # run kallisto|bustools
    # -t number of threads
    #  --cellranger cellranger output
    # -i index
    # -g t2g
    # -x technology
    # -o output
    printf "running cellranger to obtain counts\n"
    #kb count --verbose -t 6 --cellranger -i $index -g $t2g -x $tech -o ./kb_out ${fastqs[@]} --overwrite |& tee -a kallisto.LOG 
    cellranger count \
    --id cr \
    --transcriptome /home/refdata-gex-GRCh38-2020-A \
    --fastqs ./fastq \
    --sample ${1} \
    --nosecondary \
    --include-introns false \
    --no-bam |& tee -a cellranger.LOG

    mkdir out

    if [[ "$(find cr) " =~ "outs" ]] # check if counts obtained correctly
    then
        # filter empty droplets using EmptyDroplets, FDR<0.1
        # filter doublets using scDblFinder
        # output: folder with filtered cells
        printf "filtering empty droplets with FDR threshold of 0.1\n"
        filter_empty_v2.R ./cr/outs/raw_feature_bc_matrix/ 0.1 |& tee -a filter_empty.LOG
    else
        printf "ERROR: cellranger failed obtaining the counts. EXITING"
        exit 1
    fi

    if [[ "$(find out) " =~ "counts_filtered" ]]
    then
        # remove fastq and bus files
        printf "Files preprocessed correctly\nRemove fastq and bus files to save hard drive space\nReturn in root folder\n"
        mv $(pwd)/cr/outs/raw_feature_bc_matrix/ out/
        mv $(pwd)/cr/outs/web_summary.html out/
        mv $(pwd)/cr/outs/molecule_info.h5 out/
        rm -fr cr
        rm -fr fastq/
        find . -name "*.fastq" -delete
        find . -name "*.bam" -delete
        find . -name "*.fastq.gz" -delete
    else
        printf "ERROR: Cells not correctly filtered. EXITING"
        exit 1
    fi
}