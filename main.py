import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, normalize, LabelEncoder, MinMaxScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

if __name__ == '__main__':
    data = pd.ExcelFile('PROPPR_documents/PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx')
    df = pd.read_excel(data, 'PROPPR_longitudinal_dataset_edm')
    df = df[df['TIMEPOINT'] == 0]
    selected_cols = list(df.iloc[:, 7:50])
    data = df[selected_cols]
    print(data)

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp = imp.fit(data)
    data = imp.transform(data)
    # data = data.fillna(0)
    print(data)


    scalar = StandardScaler()
    std = scalar.fit_transform(data)
    pca = PCA()
    pca.fit(std)

    # Determine how many compnents.
    plt.figure(figsize=(10, 8))
    plt.plot(range(1, 44), pca.explained_variance_ratio_.cumsum(), marker='o', linestyle='--')
    plt.title("Explained Variance by Components")
    plt.xlabel("Number of Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.show()

    pca = PCA(n_components=6)
    pca.fit(std)
    scores_pca = pca.transform(std)

    # Determine how many clusters.
    wcss = []
    for i in range(1, 21):
        kmeans_pca = KMeans(n_clusters=i, init='k-means++', random_state=42)
        kmeans_pca.fit(scores_pca)
        wcss.append(kmeans_pca.inertia_)

    plt.figure(figsize=(10, 8))
    plt.plot(range(1, 21), wcss, marker='o', linestyle='--')
    plt.xlabel("Number of clusters")
    plt.ylabel("WCSS")
    plt.title("K-means with PCA Clustering")
    plt.show()

    kmeans_pca = KMeans(n_clusters=2, init='k-means++', random_state=42)
    kmeans_pca.fit(scores_pca)

    df_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
    df_pca_kmeans.columns.values[-6:] = ['Component 1', 'Component 2', 'Component 3', 'Component 4', 'Component 5',
                                         'Component 6']
    df_pca_kmeans['K-means PCA'] = kmeans_pca.labels_

    df_pca_kmeans['Clusters'] = df_pca_kmeans['K-means PCA'].map({0: 'first',
                                                                 1: 'second'})
    # Plot data by PCA components
    x_axis = df_pca_kmeans['Component 3']
    y_axis = df_pca_kmeans['Component 1']
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=x_axis, y=y_axis, hue=df_pca_kmeans['Clusters'])
    plt.title("Clusters by PCA Components")
    plt.show()
