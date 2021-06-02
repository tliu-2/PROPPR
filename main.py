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


def remove_cols(df):
    """
    Calculates percent missingness of each column in the dataframe and removes a column if its missing more than 20%
    :param df: dataframe whose columns are going to be examined for missingness
    :return: a dataframe with any column with more than 20% missingness removed.
    """
    removed = list()
    # df.rename(columns=lambda x: x.strip(), inplace=True)
    for col in df.columns:
        num_na = df[col].isna().sum()
        if num_na / len(df) > 0.2:
            df = df.drop(col, axis=1)
            removed.append(col)
    return df, removed


def rand_forest(df, t):
    """
    Conducts random forest regression. Creates trees and saves them to the computer. Identifies variable importance
    in predictions and saves graphs of variable importance to the computer.
    :param df: A dataframe of biomarker data.
    :param t: The timepoint being examined.
    :return: no return
    """
    f = open('rand_forest ' + t + '.txt', 'w')
    inj_col = df['INJ_MECH']
    # target_col = df['_30DAYST_SURV'] + 1
    target_col = df['BLUNT_INJ'] + 1
    # biomarkers = range(7, 50)
    # selected_cols = list(df.iloc[:, 7:50])
    biomarkers_log = range(50, 93)
    selected_cols = list(df.iloc[:, biomarkers_log])  # log adjusted values.

    data = df[selected_cols]
    percent_missing = data.isna().sum() * 100 / len(data)
    # These columns all have a % missing > 20%
    data, to_remove = remove_cols(data)
    if not data.empty:
        for col in to_remove:
            selected_cols.remove(col)
        # selected_cols.remove('log2_Hu_IL_2__38')
        # selected_cols.remove('log2_Hu_IL_15__73')
        # selected_cols.remove('log2_Hu_IL_12_p70__75')
        # selected_cols.remove('log2_Hu_IL_17__76')
        # selected_cols.remove('log2_Hu_FGF_basic__44')
        # selected_cols.remove('log2_Hu_GM_CSF__34')
        # selected_cols.remove('log2_Hu_VEGF__45')

        # Random Forest Analysis
        features = data
        # print(data)
        imp2 = SimpleImputer(missing_values=np.nan, strategy='mean')
        imp2 = imp2.fit(features)
        features = imp2.transform(features)
        features = pd.DataFrame(features)
        features.columns = selected_cols

        labels = np.array(target_col)
        feature_list = list(features.columns)
        features = np.array(features)

        train_features, test_features, train_labels, test_labels = train_test_split(features,
                                                                                    labels, test_size=0.25,
                                                                                    random_state=42)
        sc = StandardScaler()
        train_features = sc.fit_transform(train_features)
        test_features = sc.transform(test_features)
        rf = RandomForestRegressor(n_estimators=200, random_state=42)
        rf.fit(train_features, train_labels)

        predictions = rf.predict(test_features)
        errors = abs(predictions - test_labels)
        line = 'Mean Absolute Error: ' + str(np.mean(errors))
        f.write(line)
        f.write("\n")

        mape = 100 * (errors / test_labels)

        acc = 100 - np.mean(mape)
        line = 'Accuracy: ' + str(round(acc, 2)) + '%'
        f.write(line)
        f.write("\n")

        tree = rf.estimators_[5]
        treename = 'tree_' + t + '.dot'
        pngname = 'tree_' + t + '.png'
        export_graphviz(tree, out_file='rand_forest_graphs/' + treename, feature_names=feature_list, rounded=True, precision=1)
        (graph,) = pydot.graph_from_dot_file('rand_forest_graphs/' + treename)
        graph.write_png('rand_forest_graphs/' + pngname)

        small_treename = 'snall_tree_' + t + '.dot'
        small_pngname = 'small_tree_' + t + '.png'
        rf_small = RandomForestRegressor(n_estimators=10, max_depth=3)
        rf_small.fit(train_features, train_labels)
        tree_small = rf_small.estimators_[5]
        export_graphviz(tree_small, out_file='rand_forest_graphs/' + small_treename, feature_names=feature_list, rounded=True, precision=1)
        (graph,) = pydot.graph_from_dot_file('rand_forest_graphs/' + small_treename)
        graph.write_png('rand_forest_graphs/' + small_pngname)

        importances = list(rf.feature_importances_)
        feature_importances = list()
        for feature, importance in zip(feature_list, importances):
            feature_importances.append([feature, importance])

        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
        for feature, importance in feature_importances:
            f.write(f'Var: {feature} Importance: {importance} \n')

        x_features = list()
        y_importance = list()
        for feature, importance in enumerate(feature_importances):
            x_features.append(importance[0])
            y_importance.append(importance[1])

        f.write("\n\n")
        x_val = list(range(len(y_importance)))
        plt.bar(x_val, y_importance, orientation='vertical')
        plt.xticks(x_val, x_features, rotation='vertical')
        plt.ylabel("Var Importance")
        plt.title("Variable Importance in BLUNT_INJ " + t)
        plt.tight_layout()
        plt.savefig('var_importance/' + "Variable Importance BLUNT_INJ" + t + ".png")
        plt.close()

        f.close()


def pca_kmeans(df, t):
    """
    Does principal coordinates analysis for dimension reduction and applies k-means clustering onto a dataframe of
    biomarker data. Saves any associated graphs to the machine.
    :param df: Dataframe with biomarker data.
    :param t: timepoint being examined.
    :return: no return
    """
    biomarkers_log = range(50, 93)
    selected_cols = list(df.iloc[:, biomarkers_log])
    data = df[selected_cols]

    # These columns are dropped because more than 20% of their entries are empty.
    data, to_remove = remove_cols(data)
    for col in to_remove:
        selected_cols.remove(col)

    if not data.empty:
        # selected_cols.remove('log2_Hu_IL_2__38')
        # selected_cols.remove('log2_Hu_IL_15__73')
        # selected_cols.remove('log2_Hu_IL_12_p70__75')
        # selected_cols.remove('log2_Hu_IL_17__76')
        # selected_cols.remove('log2_Hu_FGF_basic__44')
        # selected_cols.remove('log2_Hu_GM_CSF__34')
        # selected_cols.remove('log2_Hu_VEGF__45')

        # Impute for missing data entries using mean.
        range_data = range(1, len(selected_cols) + 1)
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        imp = imp.fit(data)
        data = imp.transform(data)

        # Standardize data.
        scalar = StandardScaler()
        std = scalar.fit_transform(data)
        # Fit the standardized data to PCA.
        pca = PCA()
        pca.fit(std)

        fig1_name = "Num Components " + t + ".png"
        # Determine how many compnents.
        plt.figure(figsize=(10, 8))
        plt.plot(range_data, pca.explained_variance_ratio_.cumsum(), marker='o', linestyle='--')
        plt.title("Explained Variance by Components " + t)
        plt.xlabel("Number of Components")
        plt.ylabel("Cumulative Explained Variance")
        plt.grid(True)
        plt.savefig('k-means_graphs/' + fig1_name)
        plt.close()

        if t == 'timepoint_24':
            pca = PCA(n_components=2)
        else:
            pca = PCA(n_components=15)
        pca.fit(std)
        scores_pca = pca.transform(std)

        features = range(pca.n_components_)
        plt.figure(figsize=(10, 8))
        plt.bar(features, pca.explained_variance_ratio_)
        plt.xlabel('PCA Components')
        plt.ylabel('Variance %')
        plt.xticks(features)
        plt.title('PCA Components and Variance ' + t)
        plt.savefig('k-means_graphs/' + 'PCA Components Variance % ' + t + ".png")
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
        plt.title("K-means with PCA Clustering " + t)
        plt.xticks(range(1, 21))
        plt.savefig('k-means_graphs/' + fig2_name)
        plt.close()

        fig3_name = "Silhouette Scores " + t + ".png"
        # Silhouette Scores:
        fig, ax = plt.subplots(2, 2, figsize=(15, 8))
        for i in [2, 3, 4, 5]:
            km = KMeans(n_clusters=i, init='k-means++', n_init=10, max_iter=100, random_state=42)
            q, mod = divmod(i, 2)
            visualizer = SilhouetteVisualizer(km, colors='yellowbrick', ax=ax[q - 1][mod])
            visualizer.fit(std)
        plt.savefig('k-means_graphs/' + fig3_name)
        plt.close()

        # Do k-means
        kmeans_pca = KMeans(n_clusters=2, init='k-means++', random_state=42)
        kmeans_pca.fit(scores_pca)
        if t == 'timepoint_24':
            df_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
            df_pca_kmeans.columns.values[-2:] = ['Component 0', 'Component 1']

        else:
            df_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
            df_pca_kmeans.columns.values[-15:] = ['Component 0', 'Component 1', 'Component 2', 'Component 3',
                                                  'Component 4',
                                                  'Component 5', 'Component 6', 'Component 7', 'Component 8',
                                                  'Component 9',
                                                  'Component 10', 'Component 11', 'Component 12', 'Component 13',
                                                  'Component 14']


        df_pca_kmeans['K-means PCA'] = kmeans_pca.labels_
        df_pca_kmeans['Clusters'] = df_pca_kmeans['K-means PCA'].map(
            {0: 'first', 1: 'second'})  # , 2: 'third', 3: 'fourth',
        # 4: 'fifth'})
        df_pca_kmeans = df_pca_kmeans.sort_values('Clusters')

        fig4_name = "K-Means Clustering with PCA " + t + ".png"
        # Plot data by PCA components
        x_axis = df_pca_kmeans['Component 0']
        y_axis = df_pca_kmeans['Component 1']
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=x_axis, y=y_axis, hue=df_pca_kmeans['INJ_MECH'], style=df_pca_kmeans['Clusters'], palette='tab10')
        plt.title("Clusters by PCA Components " + t)
        # handles, labels = plt.gca().get_legend_handles_labels()
        # order = [0, 1, 2, 3, 4, 5, 6]
        # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])
        plt.savefig('k-means_graphs/' + fig4_name)
        plt.close()


def compare_blunt_pen(df, t):
    """
    Compares biomarker data between blunt injury patients and penetrating injury patients. Tests for statistical
    significant differences through t-tests.
    :param df: dataframe of biomarker data.
    :param t: timepoint being examined.
    :return: no return
    """
    f = open('t-tests ' + t + '.txt', 'w')
    df_blunt = df.loc[df['INJ_MECH'] == 'Blunt Injury Only']
    df_pen = df.loc[df['INJ_MECH'] == 'Penetrating Injury Only']
    biomarkers = range(7, 50)
    biomarkers_log = range(50, 93)
    selected_cols = list(df.iloc[:, biomarkers_log])

    df_blunt = df_blunt[selected_cols]
    df_pen = df_pen[selected_cols]

    df_blunt, to_remove = remove_cols(df_blunt)
    df_pen, to_remove2 = remove_cols(df_pen)
    for col in to_remove:
        selected_cols.remove(col)
    # Test if vars are in a normal distribution.
    for var in selected_cols:
        if var in df_blunt.columns and var in df_pen.columns:
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
        if var in df_blunt.columns and var in df_pen.columns:
            temp_b = df_blunt[var].dropna()
            temp_p = df_pen[var].dropna()
            stat, p = st.levene(temp_b, temp_p, center='median')
            if p < 0.05:
                line = var + " has equal variance \n"
                f.write(line)
                eq_variance_vars.append(var)

    f.write("\n")
    # T-Test
    for var in selected_cols:
        if var in df_blunt.columns and var in df_pen.columns:
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

                graph_name = var + '_' + t + '_t-test'
                temp_bar = {var + '_blunt': np.mean(df_blunt[var]), var + '_pen': np.mean(df_pen[var])}
                names = list(temp_bar.keys())
                values = list(temp_bar.values())
                blunt_std = np.std(temp_b)
                pen_std = np.std(temp_p)
                error = [blunt_std, pen_std]

                fig, ax = plt.subplots()
                ax.bar(names, values, yerr=error, capsize=10, error_kw={'markeredgewidth': 2})
                ax.set_ylabel('log concentration')
                ax.set_title("Comparison of " + var + '_' + t)
                plt.savefig('t-test_bar_graphs/' + graph_name + '.png')
                plt.close()
            else:
                f.write("\n\n")

    # df_blunt_st = df_blunt[selected_cols]
    # df_pen_st = df_pen[selected_cols]
    # df_blunt_st.describe().to_csv('blunt_description.csv')
    # df_pen_st.describe().to_csv('pen_description.csv')
    f.close()


def get_patient_make_up(df, t):
    """
    Get basic information on the patient population such as % male, female and % type of injury.
    :param df: A dataframe that holds information on the patients in the set.
    :param t: timepoint being examined.
    :return: no return.
    """
    total_patients = len(df['biomarker_key_PROPPRID'])
    blunt_patients = len(df.loc[df['INJ_MECH'] == 'Blunt Injury Only'])
    pen_patients = len(df.loc[df['INJ_MECH'] == 'Penetrating Injury Only'])
    perc_blunt = blunt_patients / total_patients
    perc_pen = pen_patients / total_patients

    male_patients = len(df.loc[df['female1'] == 0])
    female_patients = len(df.loc[df['female1'] == 1])
    perc_male = male_patients / total_patients
    perc_female = female_patients / total_patients

    print("For ", t)
    print("Total patients in set: ", total_patients)
    print("Percent blunt: ", perc_blunt)
    print("Percent pen: ", perc_pen)
    print("Percent male: ", perc_male)
    print("Percent female: ", perc_female)


if __name__ == "__main__":
    org_data = pd.ExcelFile('PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx')
    timepoints = ['timepoint_0', 'timepoint_2', 'timepoint_4', 'timepoint_6', 'timepoint_12', 'timepoint_24',
                  'timepoint_48', 'timepoint_72']

    for t in timepoints:
        df_org = pd.read_excel(org_data, t)
        get_patient_make_up(df_org, t)
        compare_blunt_pen(df_org, t)
        pca_kmeans(df_org, t)
        rand_forest(df_org, t)
        print("Done with: ", t)
