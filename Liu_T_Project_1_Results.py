#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import silhouette_score
from sklearn.model_selection import train_test_split

data = pd.ExcelFile('PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx')
df = pd.read_excel(data, 'PROPPR_longitudinal_dataset_edm')
df = df[df['TIMEPOINT'] == 0]
# selected_cols = list(df.iloc[:, 7:50])
selected_cols = list(df.iloc[:, 50:93])
data = df[selected_cols]
data = data.drop(['log2_Hu_IL_2__38', 'log2_Hu_IL_15__73', 'log2_Hu_IL_12_p70__75', 'log2_Hu_IL_17__76',
                 'log2_Hu_FGF_basic__44', 'log2_Hu_GM_CSF__34', 'log2_Hu_VEGF__45'], axis=1)
percent_missing = data.isna().sum() * 100 / len(data)
percent_missing


# In[4]:


# Impute for missing data entries using mean.
imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp = imp.fit(data)
data = imp.transform(data)
# data = data.fillna(0)
# print(data)

# Standardize data.
scalar = StandardScaler()
std = scalar.fit_transform(data)
pca = PCA()
pca.fit(std)

# Determine how many compnents.
plt.figure(figsize=(10, 8))
plt.plot(range(1, 37), pca.explained_variance_ratio_.cumsum(), marker='o', linestyle='--')
plt.title("Explained Variance by Components")
plt.xlabel("Number of Components")
plt.ylabel("Cumulative Explained Variance")
plt.grid()
plt.show()


# In[5]:


pca = PCA(n_components=15)
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
plt.grid()
plt.show()


# In[160]:


range_clusters = range(2, 21)
for k in range_clusters:
    clusterer = KMeans(n_clusters=k)
    preds = clusterer.fit_predict(data)
    centers = clusterer.cluster_centers_
    score = silhouette_score(data, preds)
    print("For n_clusters = {}, silhouette score is {})".format(k, score))


# In[6]:


# Do k-means
kmeans_pca = KMeans(n_clusters=2, init='k-means++', random_state=42)
kmeans_pca.fit_predict(scores_pca)

score = silhouette_score(data, kmeans_pca.labels_, metric='euclidean')
print('Silhouetter Score: %.3f' % score)

df_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
df_pca_kmeans.columns.values[-15:] = ['Component 1', 'Component 2', 'Component 3', 'Component 4', 'Component 5',
                                     'Component 6', 'Component 7', 'Component 8', 'Component 9', 'Component 10', 
                                     'Component 11', 'Component 12', 'Component 13', 'Component 14', 'Component 15']
df_pca_kmeans['K-means PCA'] = kmeans_pca.labels_

df_pca_kmeans['Clusters'] = df_pca_kmeans['K-means PCA'].map({0: 'first', 1: 'second'}) #, 2: 'third', 3: 'fourth',
                                                             # 4: 'fifth'})
df_pca_kmeans.head()


# In[7]:


# Plot data by PCA components
x_axis = df_pca_kmeans['Component 1']
y_axis = df_pca_kmeans['Component 2']
plt.figure(figsize=(10, 8))
sns.scatterplot(x=x_axis, y=y_axis, hue=df_pca_kmeans['INJ_MECH'])
plt.title("Clusters by PCA Components")
plt.show()


# Permanova, LDA, Shannon Index, 
