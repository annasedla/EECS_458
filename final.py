# importing the libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
import seaborn as sns
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot

pd.set_option('display.max_columns', 100)  # or 1000


def main():

    # CREATE INSTANCES
    labelencoder_Y = LabelEncoder()
    sc = StandardScaler()

    # IMPORT DATA SET
    dataset = pd.read_csv("breast-cancer-wisconsin.data")
    # print(dataset.head())
    print("Cancer data set dimensions : {}".format(dataset.shape))

    # CLEANUP
    dataset = dataset.drop(["bare_nucleoli"], axis = 1)

    for label in ['clump_thickness', 'size_uniformity', 'shape_uniformity', 'marginal_adhesion']:
        dataset[label] = dataset[label].fillna(method='ffill')

    # FEATURE CORRELATION
    # graph = dataset.drop(["id"], axis = 1)
    # graph_benign = graph[graph["class"] != 4]
    # graph_benign = graph_benign.drop(["class"], axis = 1)
    # graph_malignant = graph[graph["class"] != 2]
    # graph_malignant = graph_malignant.drop(["class"], axis = 1)
    # sns.set(style="ticks", color_codes = True)
    # plt.figure(figsize=(14, 12))
    # benign = sns.heatmap(graph_benign.astype(float).corr(), linewidths=0.1, square=True, linecolor="white", annot = True, cmap="YlGnBu")
    # benign.set_xticklabels(benign.get_xticklabels(), rotation=45, fontsize=8)
    # # plt.show()
    #
    # plt.figure(figsize=(14, 12))
    # mal = sns.heatmap(graph_malignant.astype(float).corr(), linewidths=0.1, square=True, linecolor="white", annot = True)
    # mal.set_xticklabels(mal.get_xticklabels(), rotation=45, fontsize=8)
    # # plt.show()

    fig = plt.figure()
    ax = sns.countplot(x="shape_uniformity", hue ="class", data=dataset, palette="Set3")
    ax.set(xlabel="Shape Uniformity", ylabel="Number of Cases")
    fig.suptitle("Shape Uniformity vs Class", y=0.96)
    plt.show()

    # TRAIN VARIOUS CLASSIFIERS
    X = dataset.iloc[:, 1:8].values
    Y = dataset.iloc[:, 9].values
    Y = labelencoder_Y.fit_transform(Y)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25, random_state=0)
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)

    #LOG REG
    classifier = LogisticRegression(random_state=0, solver="lbfgs")
    scores = cross_val_score(classifier, X, Y, cv=5, scoring='f1')
    print(scores)
    print("logreg accruacy: ", np.mean(scores))

    # SVC
    svm_classifier = SVC(kernel='linear', random_state=0, gamma="scale", C =2, probability=True)
    # classifier.fit(X_train, Y_train)
    scores = cross_val_score(svm_classifier, X, Y, cv=5, scoring='f1')
    print(scores)
    print("SVM: ", np.mean(scores))

    # DECISION TREES
    classifier = DecisionTreeClassifier(criterion='entropy', random_state=0)
    scores = cross_val_score(classifier, X, Y, cv=5, scoring='f1')
    print(scores)
    print("Decision tree: ", np.mean(scores))

    # K-NN
    classifier = KNeighborsClassifier(n_neighbors=5, metric='minkowski', p=2)
    scores = cross_val_score(classifier, X, Y, cv=5, scoring='f1')
    print(scores)
    print("KNN: ", np.mean(scores))

    # NB
    classifier = GaussianNB()
    scores = cross_val_score(classifier, X, Y, cv=5, scoring='f1')
    print(scores)
    print("Naive Bayes: ", np.mean(scores))

    # Random Forest
    classifier = RandomForestClassifier(n_estimators=10, criterion='entropy', random_state=0)
    scores = cross_val_score(classifier, X, Y, cv=5, scoring='f1')
    print(scores)
    print("Random Forest: ", np.mean(scores))

    # CONFUCION METRICS
    # Y_pred = classifier.predict(X_test)
    # cm = confusion_matrix(Y_test, Y_pred)

    # ROC curve
    # generate a no skill prediction (majority class)
    ns_probs = [0 for _ in range(len(Y_test))]
    # fit a model
    svm_classifier.fit(X_train, Y_train)
    # predict probabilities
    lr_probs = svm_classifier.predict_proba(X_test)
    # keep probabilities for the positive outcome only
    lr_probs = lr_probs[:, 1]
    # calculate scores
    ns_auc = roc_auc_score(Y_test, ns_probs)
    lr_auc = roc_auc_score(Y_test, lr_probs)
    # summarize scores
    print('No Skill: ROC AUC=%.3f' % (ns_auc))
    print('SVM: ROC AUC=%.3f' % (lr_auc))
    # calculate roc curves
    ns_fpr, ns_tpr, _ = roc_curve(Y_test, ns_probs)
    lr_fpr, lr_tpr, _ = roc_curve(Y_test, lr_probs)
    # plot the roc curve for the model
    pyplot.plot(ns_fpr, ns_tpr, linestyle='--', label='No Skill')
    pyplot.plot(lr_fpr, lr_tpr, marker='.', label='SVM')
    # axis labels
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    # show the legend
    pyplot.legend()
    # show the plot
    # pyplot.show()



if __name__ == "__main__":
    main()
