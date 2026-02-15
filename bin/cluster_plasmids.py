#!/usr/bin/env python3
"""
Cluster plasmids based on Mash distance matrix and generate network visualization.
Part of the PlasmidScope pipeline.
"""

import argparse
import sys
import csv
from collections import defaultdict


def parse_mash_distances(distance_file):
    """Parse Mash distance matrix (tab-separated)."""
    distances = {}
    samples = []

    with open(distance_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        samples = header[1:]  # First column is row label

        for row in reader:
            query = row[0]
            for i, dist in enumerate(row[1:]):
                try:
                    d = float(dist)
                    ref = samples[i]
                    if query != ref:
                        pair = tuple(sorted([query, ref]))
                        distances[pair] = d
                except (ValueError, IndexError):
                    continue

    return samples, distances


def cluster_plasmids(samples, distances, threshold):
    """Simple single-linkage clustering based on distance threshold."""
    # Initialize each sample as its own cluster
    cluster_map = {s: i for i, s in enumerate(samples)}
    next_cluster = len(samples)

    # Merge clusters for pairs below threshold
    for (s1, s2), dist in sorted(distances.items(), key=lambda x: x[1]):
        if dist <= threshold:
            c1 = cluster_map[s1]
            c2 = cluster_map[s2]
            if c1 != c2:
                # Merge: assign all members of c2 to c1
                for s in samples:
                    if cluster_map[s] == c2:
                        cluster_map[s] = c1

    # Renumber clusters sequentially
    unique_clusters = sorted(set(cluster_map.values()))
    remap = {old: new for new, old in enumerate(unique_clusters, 1)}
    cluster_map = {s: remap[c] for s, c in cluster_map.items()}

    return cluster_map


def generate_network_html(samples, distances, clusters, threshold, output_file):
    """Generate an interactive HTML network visualization using vis.js."""

    # Color palette for clusters
    colors = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
        '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9'
    ]

    # Build nodes
    nodes = []
    for s in samples:
        c = clusters.get(s, 0)
        color = colors[c % len(colors)]
        label = s.split('/')[-1].replace('.fasta', '')
        nodes.append(f'{{id: "{s}", label: "{label}", color: "{color}", '
                     f'title: "Cluster {c}", group: {c}}}')

    # Build edges (only for pairs below threshold)
    edges = []
    for (s1, s2), dist in distances.items():
        if dist <= threshold:
            width = max(1, int((1 - dist / threshold) * 5))
            edges.append(f'{{from: "{s1}", to: "{s2}", value: {1 - dist:.4f}, '
                         f'title: "Distance: {dist:.6f}", width: {width}}}')

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>PlasmidScope - Plasmid Network</title>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        h1 {{ color: #333; text-align: center; }}
        #network {{ width: 100%; height: 700px; border: 1px solid #ccc; background: white; border-radius: 8px; }}
        .info {{ text-align: center; color: #666; margin: 10px 0; }}
    </style>
</head>
<body>
    <h1>PlasmidScope - Plasmid Similarity Network</h1>
    <p class="info">Mash distance threshold: {threshold} | Nodes: {len(samples)} | Edges: {len(edges)} | Clusters: {len(set(clusters.values()))}</p>
    <div id="network"></div>
    <script>
        var nodes = new vis.DataSet([{', '.join(nodes)}]);
        var edges = new vis.DataSet([{', '.join(edges)}]);
        var container = document.getElementById('network');
        var data = {{ nodes: nodes, edges: edges }};
        var options = {{
            physics: {{ stabilization: {{ iterations: 200 }}, barnesHut: {{ gravitationalConstant: -3000 }} }},
            interaction: {{ hover: true, tooltipDelay: 100 }},
            nodes: {{ shape: 'dot', size: 16, font: {{ size: 12 }} }},
            edges: {{ smooth: {{ type: 'continuous' }} }}
        }};
        var network = new vis.Network(container, data, options);
    </script>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(description='Cluster plasmids by Mash distance')
    parser.add_argument('--distances', required=True, help='Mash distance matrix TSV')
    parser.add_argument('--threshold', type=float, default=0.04, help='Distance threshold')
    parser.add_argument('--output', required=True, help='Output clusters TSV')
    parser.add_argument('--network', help='Output network HTML file')
    args = parser.parse_args()

    samples, distances = parse_mash_distances(args.distances)
    clusters = cluster_plasmids(samples, distances, args.threshold)

    # Write clusters
    with open(args.output, 'w') as f:
        f.write("sample\tcluster\n")
        for s in sorted(samples):
            f.write(f"{s}\t{clusters[s]}\n")

    print(f"Clustered {len(samples)} plasmids into {len(set(clusters.values()))} clusters "
          f"at threshold {args.threshold}")

    # Generate network
    if args.network:
        generate_network_html(samples, distances, clusters, args.threshold, args.network)
        print(f"Network visualization saved to {args.network}")


if __name__ == '__main__':
    main()
