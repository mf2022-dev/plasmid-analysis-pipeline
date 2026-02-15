process MASH_DIST {
    tag "clustering"
    publishDir "${params.outdir}/07_comparative/mash", mode: 'copy'

    input:
    path plasmids

    output:
    path "mash_distances.tsv",    emit: distances
    path "mash_clusters.tsv",     emit: clusters
    path "plasmid_network.html",  emit: network, optional: true

    script:
    """
    # Sketch all plasmids
    mash sketch -o plasmid_sketches ${plasmids}

    # Calculate pairwise distances
    mash dist plasmid_sketches.msh plasmid_sketches.msh -t > mash_distances.tsv

    # Cluster plasmids
    cluster_plasmids.py \
        --distances mash_distances.tsv \
        --threshold ${params.mash_threshold} \
        --output mash_clusters.tsv \
        --network plasmid_network.html
    """
}
