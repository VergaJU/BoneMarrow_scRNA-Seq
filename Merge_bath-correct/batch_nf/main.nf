!/usr/bin/nextflow

ch_input = Channel.fromPath("data/*.h5ad").map{ it -> [it.baseName, it]}

// stage, split parameters for each
hvg_stage = params.hvg? params.hvg.split(',').collect{it.trim()} : []
knn_stage = params.knn? params.knn.split(',').collect{it.trim()} : []
latent_stage = params.latent ? params.latent.split(',').collect{it.trim()} : []

// place in channel
hvg_ch = Channel.value(hvg_stage).flatten()
knn_ch = Channel.value(knn_stage).flatten()
latent_ch = Channel.value(latent_stage).flatten()

// Cartesian product
scvi_params = hvg_ch.combine(latent_ch)
scanvi_params = hvg_ch.combine(latent_ch)

process RUN_SCVI{
        publishDir "./data/scVI/latent"
        tag "scVI $database $hvg hvgs $latent latent dims"
        
        input:
        tuple val(base), file(database), val(hvg), val(latent) from ch_input.combine(scvi_params)

        output:
        set val(base), val(hvg), val(latent), file('*_scVI_latent.h5ad') into UMAPvi

        script:
        """
        run_scvi.py\
                --input_object ${database}\
                        --output_prefix ${base}_${hvg}_${latent}\
                        --batch_key ${params.batch_key}\
                --hvg ${hvg}\
                --latent ${latent}
        """
}

process umap_SCVI{
        publishDir "./data/scVI/umap"
        tag "umap $datain"
        
        input:
        set val(knn), val(base), val(hvg), val(latent), file(datain) from knn_ch.combine(UMAPvi)

        output:
        set val(base), val(hvg), val(latent), val(knn), file('*.csv') into EVALvi

        script:
        """
        run_umap.py\
                --input_object ${datain}\
                --output_prefix ${base}_${hvg}_${latent}_${knn}\
                --knn $knn\
                        --batch_key ${params.batch_key}\
                        --celltype_key ${params.celltype_key}
        """
}

process eval_SCVI{
        publishDir "./data/scVI/eval"
        tag "eval $datain"
        
        input:
        set val(base), val(hvg), val(latent), val(knn), file(datain) from EVALvi

        output:
        stdout to out
        file('*_scVI_eval.csv')

        script:
        """
        eval.R\
                --input_object ${datain}\
                --output_prefix ${base}_${hvg}_${latent}_${knn}\
                        --batch_key ${params.batch_key}\
                        --celltype_key ${params.celltype_key}
        """
}

process RUN_SCANVI{
        publishDir "./data/scANVI/latent"
        tag "scVI $database $hvg hvgs $latent latent dims"
        
        input:
        tuple val(base), file(database),i val(hvg), val(latent) from ch_anvi_input.combine(scanvi_params)

        output:
        set val(base), val(hvg), val(latent), file('*_scVI_latent.h5ad') into UMAPanvi

        script:
        """
        run_scvi.py\
                --input_object ${database}\
                        --output_prefix ${base}_${hvg}_${latent}\
                        --batch_key ${params.batch_key}\
                --celltype_key ${params.celltype_key}\
                --hvg ${hvg}\
                --latent ${latent}
        """
}

process umap_SCANVI{
        publishDir "./data/scANVI/umap"
        tag "umap $datain"
        
        input:
        set val(knn), val(base), val(hvg), val(latent), file(datain) from knn_anvi_ch.combine(UMAPanvi)

        output:
        set val(base), val(hvg), val(latent), val(knn), file('*.csv') into EVALanvi

        script:
        """
        run_umap.py\
                --input_object ${datain}\
                --output_prefix ${base}_${hvg}_${latent}_${knn}\
                --knn $knn\
                        --batch_key ${params.batch_key}\
                        --celltype_key ${params.celltype_key}
        """
}

process eval_SCANVI{
        publishDir "./data/scVI/eval"
        tag "eval $datain"
        
        input:
        set val(base), val(hvg), val(latent), val(knn), file(datain) from EVALanvi

        output:
        stdout to out
        file('*_scANVI_eval.csv')

        script:
        """
        eval.R\
                --input_object ${datain}\
                --output_prefix ${base}_${hvg}_${latent}_${knn}\
                        --batch_key ${params.batch_key}\
                        --celltype_key ${params.celltype_key}
        """
}