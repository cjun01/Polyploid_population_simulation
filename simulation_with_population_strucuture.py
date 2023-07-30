import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
import os
import seaborn as sns
import random
from tqdm import tqdm
from numba import jit, prange
script_directory = os.path.dirname(os.path.abspath(__file__))
def simulate_and_visualize_population(num_accessions, num_markers, ploidy_level, num_clusters, min_block_size=20, max_block_size=80):
    accessions_per_cluster = num_accessions // num_clusters

    # Create names for accessions and markers
    accession_names = ["Accession_" + str(i + 1) for i in range(num_accessions)]
    marker_names = ["Marker_" + str(i + 1) for i in range(num_markers)]

    # Initialize the data array as int32 type
    data = np.random.randint(0, ploidy_level + 1, (num_accessions, num_markers))

    for i in range(num_clusters):
        chaos = random.uniform(0.1, 2)
        print(chaos)
        bias_mean = np.random.uniform(1, 1)
        bias = np.random.normal(bias_mean/chaos, 1/chaos, (1, num_markers))

        start_idx = i * accessions_per_cluster
        end_idx = (i + 1) * accessions_per_cluster

        biased_data = data[start_idx:end_idx] + np.round(bias)
        data[start_idx:end_idx] = np.clip(biased_data, 0, ploidy_level)

        # Introduce LD blocks within the cluster
        current_marker = 0
        cluster_data = data[start_idx:end_idx]
        while current_marker < num_markers:
            block_size = np.random.randint(min_block_size, max_block_size + 1)

            start_marker_idx = current_marker
            end_marker_idx = start_marker_idx + block_size
            if end_marker_idx > num_markers:
                end_marker_idx = num_markers

            key_snp = cluster_data[:, start_marker_idx]

            for j in range(start_marker_idx + 1, end_marker_idx):
                noise = np.random.normal(0, 0.5, accessions_per_cluster).astype(int)
                cluster_data[:, j] = key_snp + noise
                cluster_data[:, j] = np.clip(cluster_data[:, j], 0, ploidy_level)

            current_marker += block_size

        data[start_idx:end_idx] = cluster_data
    # Intermixing 5% of the data between each cluster
    mix_percentage = 0.005
    mix_size = int(accessions_per_cluster * mix_percentage)

    for i in range(num_clusters):
        next_cluster = (i + 1) % num_clusters
        source_indices = np.random.choice(np.arange(i * accessions_per_cluster, (i + 1) * accessions_per_cluster), mix_size, replace=False)
        target_indices = np.random.choice(np.arange(next_cluster * accessions_per_cluster, (next_cluster + 1) * accessions_per_cluster), mix_size, replace=False)
        data[source_indices], data[target_indices] = data[target_indices].copy(), data[source_indices].copy()

    # Convert the data into a pandas DataFrame
    df = pd.DataFrame(data, index=accession_names, columns=marker_names)

    # Save the dataframe to a CSV file in the script's directory

    csv_path = os.path.join(script_directory, 'hexaploid_population_data.csv')
    df.to_csv(csv_path)

    # Performing PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data)

    # Visualizing the PCA results
    fig, ax = plt.subplots()
    colors = ['r', 'g', 'b', 'c', 'm']

    for i, color in enumerate(colors):
        ax.scatter(principal_components[i * accessions_per_cluster:(i + 1) * accessions_per_cluster, 0],
                   principal_components[i * accessions_per_cluster:(i + 1) * accessions_per_cluster, 1],
                   c=color, label=f'Cluster {i + 1}', s=5)

    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.legend()
    csv_path_image=csv_path = os.path.join(script_directory, 'simulated_PCA_plot')
    plt.title('PCA of Ambiguously Simulated Hexaploid Population with Integer Genotypes')
    plt.savefig(csv_path_image, dpi=600,
                bbox_inches='tight')
    plt.show()

    return df,data


@jit(nopython=True, parallel=True)
def compute_r2_for_marker(i, matrix, num_markers):
    r2_row = np.zeros(num_markers)
    for j in prange(num_markers):
        x = matrix[:, i] - np.mean(matrix[:, i])
        y = matrix[:, j] - np.mean(matrix[:, j])
        r2_row[j] = (np.sum(x * y) / (np.sqrt(np.sum(x ** 2) * np.sum(y ** 2)))) ** 2
    return r2_row


def calculate_r2(matrix):
    num_markers = matrix.shape[1]
    r2_values = np.zeros((num_markers, num_markers))

    for i in tqdm(range(num_markers), desc="Calculating LD"):
        r2_values[i] = compute_r2_for_marker(i, matrix, num_markers)

    # Show the heatmap
    plt.figure(figsize=(10, 10))
    sns.heatmap(r2_values, cmap="YlGnBu", xticklabels=500, yticklabels=500)
    plt.title("Linkage Disequilibrium (r^2 values) Heatmap")
    csv_path = os.path.join(script_directory, 'simulated_LD_heatmap')
    plt.savefig(csv_path, dpi=600,
                bbox_inches='tight')
    plt.show()

    return r2_values
df,data = simulate_and_visualize_population(2000, 15000, 6, 5)
ld_r2 = calculate_r2(data)