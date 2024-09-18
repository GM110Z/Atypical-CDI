import pandas as pd
import sys

# Step 1: Edit MMseqs clustering table
def load_clustering_table(file_path):
    return pd.read_csv(file_path, sep="\t", header=None, names=["ClusterRep", "SingleMember"])

def assign_numerical_ids(df):
    unique_clusters = df['ClusterRep'].unique()
    cluster_to_id = {cluster: idx + 1 for idx, cluster in enumerate(unique_clusters)}
    df['ClusterRep'] = df['ClusterRep'].map(cluster_to_id)
    return df, cluster_to_id

def save_updated_table(df, file_path):
    df.to_csv(file_path, sep="\t", index=False, header=False)

# Step 2: Calculate protein overlap and generate similarity report
def load_updated_clustering_table(file_path):
    return pd.read_csv(file_path, sep="\t", header=None, names=["ClusterRep", "SingleMember"])

def extract_genome_proteins(df):
    genome_proteins = {}
    for _, row in df.iterrows():
        protein = row['SingleMember']
        genome = extract_genome(protein)
        cluster = row['ClusterRep']
        if genome not in genome_proteins:
            genome_proteins[genome] = set()
        genome_proteins[genome].add(cluster)
    return genome_proteins

def extract_genome(protein_id):
    parts = protein_id.split("_")
    return "_".join(parts[:2])

def calculate_protein_overlap(genome1_proteins, genome2_proteins):
    overlap = len(genome1_proteins.intersection(genome2_proteins))
    total_proteins = len(genome1_proteins.union(genome2_proteins))
    return overlap / total_proteins * 100

def generate_similarity_report(genome_proteins, threshold=70):
    report = []
    genome_list = list(genome_proteins.keys())
    for i in range(len(genome_list)):
        for j in range(i + 1, len(genome_list)):
            genome1 = genome_list[i]
            genome2 = genome_list[j]
            overlap = calculate_protein_overlap(genome_proteins[genome1], genome_proteins[genome2])
            if overlap >= threshold:
                report.append({
                    'Genome1': genome1,
                    'Genome2': genome2,
                    'Overlap (%)': overlap
                })
    return pd.DataFrame(report)

# Step 3: Assign genome groups based on similarity
def load_similarity_report(file_path):
    return pd.read_csv(file_path)

def assign_genome_groups(similarity_df, threshold=70):
    genome_groups = {}
    group_id = 1
    visited = set()
    for _, row in similarity_df.iterrows():
        genome1 = row['Genome1']
        genome2 = row['Genome2']
        overlap = row['Overlap (%)']
        if overlap >= threshold:
            if genome1 not in visited and genome2 not in visited:
                genome_groups[genome1] = group_id
                genome_groups[genome2] = group_id
                visited.add(genome1)
                visited.add(genome2)
                group_id += 1
            elif genome1 in visited:
                genome_groups[genome2] = genome_groups[genome1]
                visited.add(genome2)
            elif genome2 in visited:
                genome_groups[genome1] = genome_groups[genome2]
                visited.add(genome1)
    all_genomes = set(similarity_df['Genome1']).union(set(similarity_df['Genome2']))
    for genome in all_genomes:
        if genome not in genome_groups:
            genome_groups[genome] = group_id
            group_id += 1
    return genome_groups

def generate_grouping_table(genome_groups):
    df = pd.DataFrame(list(genome_groups.items()), columns=['Genome', 'Group'])
    return df

# Main workflow
if __name__ == "__main__":
    # Paths to files
    original_file_path = '/Users/giusym/Desktop/Test/clusterRes_cluster.tsv'
    updated_file_path = '/Users/giusym/Desktop/Test/MMseqs-clusternumb.tsv'
    similarity_report_path = '/Users/giusym/Desktop/Test/genome_similarity_report.csv'
    grouping_table_path = '/Users/giusym/Desktop/Test/genome_grouping_table.csv'
    
    # Step 1: Edit MMseqs clustering table
    clustering_df = load_clustering_table(original_file_path)
    updated_df, _ = assign_numerical_ids(clustering_df)
    save_updated_table(updated_df, updated_file_path)
    
    # Step 2: Calculate protein overlap and generate similarity report
    clustering_df = load_updated_clustering_table(updated_file_path)
    genome_proteins = extract_genome_proteins(clustering_df)
    similarity_report = generate_similarity_report(genome_proteins, threshold=10)  # Use threshold from sys.argv[1] if needed
    similarity_report.to_csv(similarity_report_path, index=False)
    
    # Step 3: Assign genome groups based on similarity
    similarity_df = load_similarity_report(similarity_report_path)
    genome_groups = assign_genome_groups(similarity_df, threshold=70)
    grouping_table = generate_grouping_table(genome_groups)
    grouping_table.to_csv(grouping_table_path, index=False)
    
    # Print confirmation
    print(f"Updated clustering table saved to {updated_file_path}.")
    print(f"Similarity report saved to {similarity_report_path}.")
    print(f"Genome grouping table saved to {grouping_table_path}.")
