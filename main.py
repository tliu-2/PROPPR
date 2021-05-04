import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import silhouette_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import export_graphviz
import pydot

if __name__ == "__main__":
    data = pd.ExcelFile('PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx')
    df = pd.read_excel(data, 'timepoint_0')
    inj_col = df['INJ_MECH']
    # target_col = df['_30DAYST_SURV'] + 1
    target_col = df['ARDS'] + 1
    # selected_cols = list(df.iloc[:, 7:50])
    selected_cols = list(df.iloc[:, 50:93])
    print(selected_cols)
    data = df[selected_cols]
    data = data.drop(['log2_Hu_IL_2__38', 'log2_Hu_IL_15__73', 'log2_Hu_IL_12_p70__75', 'log2_Hu_IL_17__76',
                      'log2_Hu_FGF_basic__44', 'log2_Hu_GM_CSF__34', 'log2_Hu_VEGF__45'], axis=1)
    selected_cols.remove('log2_Hu_IL_2__38')
    selected_cols.remove('log2_Hu_IL_15__73')
    selected_cols.remove('log2_Hu_IL_12_p70__75')
    selected_cols.remove('log2_Hu_IL_17__76')
    selected_cols.remove('log2_Hu_FGF_basic__44')
    selected_cols.remove('log2_Hu_GM_CSF__34')
    selected_cols.remove('log2_Hu_VEGF__45')
    print(selected_cols)
    percent_missing = data.isna().sum() * 100 / len(data)

    # Random Forest Analysis
    features = data
    print(data)
    imp2 = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp2 = imp2.fit(features)
    features = imp2.transform(features)
    features = pd.DataFrame(features)
    features.columns = selected_cols
    features = features.assign(INJ_MECH=inj_col)
    features = pd.get_dummies(features)
    print(features.head())

    labels = np.array(target_col)
    feature_list = list(features.columns)
    features = np.array(features)

    train_features, test_features, train_labels, test_labels = train_test_split(features,
                                                                                labels, test_size=0.25,
                                                                                random_state=42)
    print('Training Features Shape: ', train_features.shape)
    print('Training Labels Shape: ', train_labels.shape)
    print('Testing Features Shape: ', test_features.shape)
    print('Testing Labels Shape: ', test_labels.shape)

    rf = RandomForestRegressor(n_estimators=1000, random_state=42)
    rf.fit(train_features, train_labels)

    predictions = rf.predict(test_features)
    errors = abs(predictions - test_labels)
    print('Mean Absolute Error: ', np.mean(errors))

    mape = 100 * (errors / test_labels)

    acc = 100 - np.mean(mape)
    print('Accuracy: ', round(acc, 2), "%")
    tree = rf.estimators_[5]
    export_graphviz(tree, out_file='tree.dot', feature_names=feature_list, rounded=True, precision=1)
    (graph, ) = pydot.graph_from_dot_file('tree.dot')
    graph.write_png('tree.png')

    rf_small = RandomForestRegressor(n_estimators=10, max_depth=3)
    rf_small.fit(train_features, train_labels)
    tree_small = rf_small.estimators_[5]
    export_graphviz(tree_small, out_file='small_tree.dot', feature_names=feature_list, rounded=True, precision=1)
    (graph, ) = pydot.graph_from_dot_file('small_tree.dot')
    graph.write_png('small_tree.png')

    importances = list(rf.feature_importances_)

    # feature_importances = [(features, importances) for feature, importance in zip(feature_list, importances)]
    # feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
    # [print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances]
    for feature, importance in zip(feature_list, importances):
        print(f'Var: {feature} Importance: {importance}')

    x_val = list(range(len(importances)))
    plt.bar(x_val, importances, orientation='vertical')
    plt.xticks(x_val, feature_list, rotation='vertical')
    plt.show()

    runKmeans = False
    if runKmeans:
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

        range_clusters = range(2, 21)
        for k in range_clusters:
            clusterer = KMeans(n_clusters=k)
            preds = clusterer.fit_predict(data)
            centers = clusterer.cluster_centers_
            score = silhouette_score(data, preds)
            print("For n_clusters = {}, silhouette score is {})".format(k, score))

        # Do k-means
        kmeans_pca = KMeans(n_clusters=2, init='k-means++', random_state=42)
        kmeans_pca.fit(scores_pca)

        df_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
        df_pca_kmeans.columns.values[-15:] = ['Component 1', 'Component 2', 'Component 3', 'Component 4', 'Component 5',
                                              'Component 6', 'Component 7', 'Component 8', 'Component 9',
                                              'Component 10',
                                              'Component 11', 'Component 12', 'Component 13', 'Component 14',
                                              'Component 15']
        df_pca_kmeans['K-means PCA'] = kmeans_pca.labels_

        df_pca_kmeans['Clusters'] = df_pca_kmeans['K-means PCA'].map(
            {0: 'first', 1: 'second'})  # , 2: 'third', 3: 'fourth',
        # 4: 'fifth'})
        df_pca_kmeans.head()

        # Plot data by PCA components
        x_axis = df_pca_kmeans['Component 1']
        y_axis = df_pca_kmeans['Component 2']
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=x_axis, y=y_axis, hue=df_pca_kmeans['INJ_MECH'])
        plt.title("Clusters by PCA Components")
        plt.show()



