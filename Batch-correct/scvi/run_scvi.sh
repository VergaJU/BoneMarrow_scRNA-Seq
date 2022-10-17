#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Run parameter optimization for scVI/scANVI."
   echo
   echo "Syntax: scriptTemplate [-csv|b|c|h|]"
   echo "options:"
   echo "c     csv file containing file paths, hvg and latent space dimension for each run"
   echo "        e.g."
   echo "             path/to/file,number hvg,number latent space"
   echo "b     Batch key"
   echo "f     first key (cell type)"
   echo "s     second key (optional, another label key)"
   echo "h     Print this Help."
   echo
}

############################################################
############################################################

# Get the options
while getopts ":hf:c:b:f:s:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      c) # Enter csv file
         csv=$OPTARG;;
      b) # Batch key
         batch=$OPTARG;;
      f) # Enter a name
         first_key=$OPTARG;;
      s) # Enter a name
         second_key=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option."
         echo "run the script with -h to see usage."
         exit;;
   esac
done

if [ $OPTIND -eq 1 ]; then
	echo "Error: No options were passed."
	echo "run the script with -h to see usage."
	exit
fi

# argument:
# csv file containing file paths, latent space and high variable genes to consider:
#  filepath,batch_key,number hvg, number latent space dimensions
#

filename=$(head -n1 params | cut -d"," -f1)
hvg=$(head -n1 params | cut -d"," -f2)

./bin/pre_process.py\
 --input_object ${filename}\
 --output_object ${filename%.h5ad}tmp.h5ad\
 --batch_key ${batch}\
 --hvg ${hvg}\


input=${filename%.h5ad}tmp.h5ad

mkdir latent # create output folder for scVI

# run scVI

while IFS=, read -r filepath hvg latent; do
	file=$(echo ${filepath} | awk -F "/" '{print $NF}')
	output=$(echo latent/${file%.h5ad}_${hvg}_${latent})
	./bin/run_scvi.py\
        --input_object ${filepath}\
        --output_prefix ${output}\
        --batch_key ${batch}\
        --latent ${latent}
done < "$csv"

# run umaps


mkdir umap #create umap output folder
mkdir umap/figures #png figures
mkdir umap/csv # csv outputs

dims=(5 15 25 35 45 50)

for f in ./latent/*;do
	file=$(echo ${f} | awk -F "/" '{print $NF}')
	echo ${file}
	for k in ${dims[@]};do
        ./bin/run_umap.py\
            --input_object ${f}\
            --output_prefix ${file%latent.h5ad}${k}\
            --knn $k\
            --batch_key ${batch}\
            --celltype_key ${first_key}
	done
done

mkdir eval #eval output directory

for f in ./umap/csv/*;do
	file=$(echo ${f} | awk -F "/" '{print $NF}')
	output=$(echo eval/${file%.csv})
	./bin/eval.R \
      --input_object ${f} \
      --output_prefix ${output} \
      --batch_key ${batch} \
      --first_key ${first_key} \
      --second_key ${second_key}
done

rm -f ${filename%.h5ad}tmp.h5ad