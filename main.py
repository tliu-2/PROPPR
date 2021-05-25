import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import silhouette_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import export_graphviz
import pydot
from scipy import stats as st
from yellowbrick.cluster import SilhouetteVisualizer


def rand_forest(df):
    inj_col = df['INJ_MECH']
    # target_col = df['_30DAYST_SURV'] + 1
    target_col = df['BLUNT_INJ'] + 1
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
    # features = features.assign(INJ_MECH=inj_col)
    # features = pd.get_dummies(features)
    print(features.head())

    labels = np.array(target_col)
    feature_list = list(features.columns)
    features = np.array(features)

    train_features, test_features, train_labels, test_labels = train_test_split(features,
                                                                                labels, test_size=0.25,
                                                                                random_state=42)
    sc = StandardScaler()
    train_features = sc.fit_transform(train_features)
    test_features = sc.transform(test_features)
    print('Training Features Shape: ', train_features.shape)
    print('Training Labels Shape: ', train_labels.shape)
    print('Testing Features Shape: ', test_features.shape)
    print('Testing Labels Shape: ', test_labels.shape)

    rf = RandomForestRegressor(n_estimators=200, random_state=42)
    rf.fit(train_features, train_labels)

    predictions = rf.predict(test_features)
    errors = abs(predictions - test_labels)
    # print("errors = ", errors)
    print('Mean Absolute Error: ', np.mean(errors))

    mape = 100 * (errors / test_labels)

    acc = 100 - np.mean(mape)
    print('Accuracy: ', round(acc, 2), "%")
    tree = rf.estimators_[5]
    export_graphviz(tree, out_file='tree.dot', feature_names=feature_list, rounded=True, precision=1)
    (graph,) = pydot.graph_from_dot_file('tree.dot')
    graph.write_png('tree.png')

    rf_small = RandomForestRegressor(n_estimators=10, max_depth=3)
    rf_small.fit(train_features, train_labels)
    tree_small = rf_small.estimators_[5]
    export_graphviz(tree_small, out_file='small_tree.dot', feature_names=feature_list, rounded=True, precision=1)
    (graph,) = pydot.graph_from_dot_file('small_tree.dot')
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
    plt.ylabel("Var Importance")
    plt.show()


def pca_kmeans(df, t):
    selected_cols = list(df.iloc[:, 50:93])
    # print(selected_cols)
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

    fig1_name = "Num Components " + t + ".png"
    # Determine how many compnents.
    plt.figure(figsize=(10, 8))
    plt.plot(range(1, 37), pca.explained_variance_ratio_.cumsum(), marker='o', linestyle='--')
    plt.title("Explained Variance by Components")
    plt.xlabel("Number of Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.grid(True)
    plt.savefig(fig1_name)
    plt.close()

    pca = PCA(n_components=15)
    pca.fit(std)
    scores_pca = pca.transform(std)

    features = range(pca.n_components_)
    plt.figure(figsize=(10, 8))
    plt.bar(features,pca.explained_variance_ratio_)
    plt.xlabel('PCA Components')
    plt.ylabel('Variance %')
    plt.xticks(features)
    plt.savefig('PCA Components Variance % ' + t + ".png")
    plt.close()

    fig2_name = "Number of Clusters " + t + ".png"
    # Determine how many clusters.
    wcss = []
    for i in range(1, 21):
        kmeans_pca = KMeans(n_clusters=i, init='k-means++', random_state=42)
        kmeans_pca.fit(scores_pca)
        wcss.append(kmeans_pca.inertia_)

    plt.figure(figsize=(10, 8))
    plt.grid(True)
    plt.plot(range(1, 21), wcss, marker='o', linestyle='--')
    plt.xlabel("Number of clusters")
    plt.ylabel("WCSS")
    plt.title("K-means with PCA Clustering")
    plt.xticks(range(1, 21))
    plt.savefig(fig2_name)
    plt.close()

    # range_clusters = range(2, 21)
    # for k in range_clusters:
    #     clusterer = KMeans(n_clusters=k)
    #     preds = clusterer.fit_predict(data)
    #     centers = clusterer.cluster_centers_
    #     score = silhouette_score(data, preds)
    #     print("For n_clusters = {}, silhouette score is {})".format(k, score))

    fig3_name = "Silhouette Scores " + t + ".png"
    # Silhouette Scores:
    fig, ax = plt.subplots(2, 2, figsize=(15, 8))
    for i in [2, 3, 4, 5]:
        km = KMeans(n_clusters=i, init='k-means++', n_init=10, max_iter=100, random_state=42)
        q, mod = divmod(i, 2)
        visualizer = SilhouetteVisualizer(km, colors='yellowbrick', ax=ax[q-1][mod])
        visualizer.fit(std)
    plt.savefig(fig3_name)
    plt.close()

    # Do k-means
    kmeans_pca = KMeans(n_clusters=2, init='k-means++', random_state=42)
    kmeans_pca.fit(scores_pca)

    df_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
    df_pca_kmeans.columns.values[-15:] = ['Component 0', 'Component 1', 'Component 2', 'Component 3', 'Component 4',
                                          'Component 5', 'Component 6', 'Component 7', 'Component 8',
                                          'Component 9',
                                          'Component 10', 'Component 11', 'Component 12', 'Component 13',
                                          'Component 14']
    df_pca_kmeans['K-means PCA'] = kmeans_pca.labels_

    df_pca_kmeans['Clusters'] = df_pca_kmeans['K-means PCA'].map(
        {0: 'first', 1: 'second'})  # , 2: 'third', 3: 'fourth',
    # 4: 'fifth'})
    df_pca_kmeans.head()

    fig4_name = "K-Means Clustering with PCA " + t + ".png"
    # Plot data by PCA components
    x_axis = df_pca_kmeans['Component 0']
    y_axis = df_pca_kmeans['Component 1']
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=x_axis, y=y_axis, hue=df_pca_kmeans['Clusters'], palette='tab10')
    plt.title("Clusters by PCA Components " + t)
    plt.savefig(fig4_name)
    plt.close()


def compare_blunt_pen(df, t):
    f = open(t, 'w')
    df_blunt = df.loc[df['INJ_MECH'] == 'Blunt Injury Only']
    df_pen = df.loc[df['INJ_MECH'] == 'Penetrating Injury Only']
    selected_cols = list(df.iloc[:, 50:93])
    # selected_cols = list(df.iloc[:, 7:50])
    # df_blunt_il6 = df_blunt['log2_Hu_IL_6__19']
    # df_pen_il6 = df_pen['log2_Hu_IL_6__19']

    # df_blunt_all = df_blunt.iloc[:, 50:93]
    # df_pen_all = df_pen.iloc[:, 50:93]
    #
    # print("T-test IL-6:")
    # print(st.ttest_ind(df_blunt_il6, df_pen_il6, nan_policy='omit'))

    df_blunt = df_blunt[selected_cols]
    df_pen = df_pen[selected_cols]

    df_blunt = df_blunt.drop(['log2_Hu_IL_2__38', 'log2_Hu_IL_15__73', 'log2_Hu_IL_12_p70__75', 'log2_Hu_IL_17__76',
                                 'log2_Hu_FGF_basic__44', 'log2_Hu_GM_CSF__34', 'log2_Hu_VEGF__45'], axis=1)
    df_pen = df_pen.drop(['log2_Hu_IL_2__38', 'log2_Hu_IL_15__73', 'log2_Hu_IL_12_p70__75', 'log2_Hu_IL_17__76',
                             'log2_Hu_FGF_basic__44', 'log2_Hu_GM_CSF__34', 'log2_Hu_VEGF__45'], axis=1)
    selected_cols.remove('log2_Hu_IL_2__38')
    selected_cols.remove('log2_Hu_IL_15__73')
    selected_cols.remove('log2_Hu_IL_12_p70__75')
    selected_cols.remove('log2_Hu_IL_17__76')
    selected_cols.remove('log2_Hu_FGF_basic__44')
    selected_cols.remove('log2_Hu_GM_CSF__34')
    selected_cols.remove('log2_Hu_VEGF__45')

    # Test if vars are in a normal distribution.
    for var in selected_cols:
        temp_b = df_blunt[var]
        temp_p = df_pen[var]
        test_stat_b, p_val_b = st.normaltest(temp_b, nan_policy='omit')
        test_stat_p, p_val_p = st.normaltest(temp_p, nan_policy='omit')
        if p_val_b < 0.05:  # Null Hypothesis is that var comes from a normal distribution
            line = var + " blunt does not come from a normal distribution. \n"
            f.write(line)
        if p_val_p < 0.05:  # Null Hypothesis is that var comes from a normal distribution
            line = var + " pen does not come from a normal distribution. \n"
            f.write(line)

    f.write("\n")
    # Levene Test for equal variances.
    eq_variance_vars = []
    for var in selected_cols:
        temp_b = df_blunt[var].dropna()
        temp_p = df_pen[var].dropna()
        stat, p = st.levene(temp_b, temp_p, center='median')
        if p < 0.05:
            line = var + " has equal variance \n"
            f.write(line)
            eq_variance_vars.append(var)

    f.write("\n")
    for var in selected_cols:
        temp_b = df_blunt[var]
        temp_p = df_pen[var]
        line = "T-test on: " + var + "\n"
        f.write(line)
        if var in eq_variance_vars:
            # print("Doing equal variance t-test")
            t_stat, pval = st.ttest_ind(temp_b, temp_p, equal_var=True, nan_policy='omit')
        else:
            t_stat, pval = st.ttest_ind(temp_b, temp_p, equal_var=False, nan_policy='omit')

        s_p_val = str(pval)
        s_t_stat = str(t_stat)

        t_stat_line = "T statistic: " + s_t_stat + "\n"
        p_val_line = "p-value: " + s_p_val + "\n"
        f.write(t_stat_line)
        f.write(p_val_line)
        if pval < 0.05 / 36:  # Bonferroni correction - divide the threshold by number of tests - in this case 36.
            f.write("Meets threshold for statistical significance. \n\n")
        else:
            f.write("\n\n")

    # df_blunt_st = df_blunt[selected_cols]
    # df_pen_st = df_pen[selected_cols]
    # df_blunt_st.describe().to_csv('blunt_description.csv')
    # df_pen_st.describe().to_csv('pen_description.csv')
    f.close()


if __name__ == "__main__":
    org_data = pd.ExcelFile('PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx')
    timepoints = ['timepoint_0', 'timepoint_2', 'timepoint_4', 'timepoint_6', 'timepoint_12', 'timepoint_24',
                  'timepoint_48', 'timepoint_72']

    for t in timepoints:
        df = pd.read_excel(org_data, t)
        compare_blunt_pen(df, t)
        pca_kmeans(df, t)
        print("Done with: ", t)

