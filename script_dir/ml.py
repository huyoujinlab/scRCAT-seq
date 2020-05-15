from sklearn.utils import shuffle
import sys
import warnings
from plotnine import *
import re 
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import os
import glob
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import cross_val_score
cv = ShuffleSplit(n_splits=5, test_size=0.3, random_state=0)
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split
import argparse
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.externals import joblib

parser = argparse.ArgumentParser()
parser.add_argument('--train', '-t', help='train data path')
parser.add_argument('--test', '-e', help='test data path')
parser.add_argument('--model', '-m', help="existed traning model")
parser.add_argument('--threshold', '-n', default = 0, help="threshold value")
args = parser.parse_args()

train = None
input_file = args.test
print(input_file)
if ("3tail" in input_file) | ("tes" in input_file):
    tss = "tes"
    label = "isinPAS"
elif ("5cap" in input_file) | ("5cap" in input_file):
    tss = "tss"
    label = "isindsc"
output_file_gene = input_file.split(".csv")[0] + "_gene.csv"
output_file_intergenic = input_file.split(".csv")[0] + "_intergenic.csv"
outputdir = os.path.join("output", re.search("(tc_)([a-zA-Z0-9]+)",input_file).group())

tmp_dfs = pd.read_csv(input_file)
cols = tmp_dfs.filter(regex="^dominant_t\w+_RPM").columns[0]
tmp_dfs = tmp_dfs.loc[tmp_dfs[cols] >=int(args.threshold)]
#remove miRNA 
tmp_df = tmp_dfs[tmp_dfs["gene_type"] != "miRNA"]
gene = tmp_df.loc[(tmp_df["gene"] != "intergenic")]
intergenic = tmp_df.loc[(tmp_df["gene"] == "intergenic")]
gene.to_csv(output_file_gene, index=False)
intergenic.to_csv(output_file_intergenic, index=False)

tmp = tmp_df.columns.values
if tss == "tes":
    features_basic = tmp[[5,7,9]]
    features_motif = tmp[13:-15]
    features_internal = tmp[-15:-1]
elif tss == "tss":
    features_basic = tmp[[5,7,9]]
    features_motif = tmp[13:-7]
    features_internal = tmp[-7:-1]

print("first peak0:\t", len(tmp_dfs.loc[tmp_dfs[label] == 0]))
print("first peak1:\t", len(tmp_dfs.loc[tmp_dfs[label] == 1]))
print("remove miRNA peak0:\t", len(tmp_df.loc[tmp_df[label] == 0]))
print("remove miRNA peak1:\t", len(tmp_df.loc[tmp_df[label] == 1]))
print("remove miRNA all peak:\t", len(tmp_df))
print("gene peak0:\t", len(gene.loc[gene[label] == 0]))
print("gene peak1:\t", len(gene.loc[gene[label] == 1]))
print("gene number:\t", len(gene.loc[gene.gene_type == "protein_coding"].gene.unique()))

class MajorityVoteClassifier(BaseEstimator):
    def __init__(self, C = 1, method = "vote"):
        self.C = C
        self.method = method
        return
    def predict(self, X):
        if self.method == "vote":
            y_knn_pred  = self.knn.predict(X)
            y_rf_pred  = self.rf.predict(X)
            y_svc_pred = self.svc.predict(X)
            y_majority_vote = (y_knn_pred & y_rf_pred) | (y_knn_pred & y_svc_pred) | (y_rf_pred & y_svc_pred)
            return y_majority_vote
        elif self.method == "lr":
            y_lr_pred  = self.lr.predict(X)
            return y_lr_pred
        elif self.method == "rf":
            y_rf_pred  = self.rf.predict(X)
            return y_rf_pred
        elif self.method == "svm":
            y_svm_pred = self.svm.predict(X)
            return y_svm_pred
        elif self.method == "svc":
            y_svc_pred = self.svc.predict(X)
            return y_svc_pred
        elif self.method == "knn":  #add knn model for prediction
            y_knn_pred = self.knn.predict(X)
            return y_knn_pred
    def fit(self, X, y):
        if self.method == "vote":
            self.knn = KNeighborsClassifier(n_neighbors=25,n_jobs = -1, leaf_size = 50)
            self.svc =  svm.SVC(C=self.C, gamma="scale")
            self.rf = RandomForestClassifier(n_estimators=10, random_state=1, max_depth = 10, n_jobs = -1, max_features = None)
            self.knn.fit(X, y)
            self.svc.fit(X, y)
            self.rf.fit(X, y)
        elif self.method == "lr":
            self.lr = LogisticRegression(C=self.C, random_state=1, solver="lbfgs", max_iter=1000000, n_jobs = -1)
            self.lr.fit(X, y)
        elif self.method == "rf":
            self.rf = RandomForestClassifier(n_estimators=100, random_state=1, max_samples = 0.7, max_depth = 10, n_jobs = -1, max_features = None)
            self.rf.fit(X, y)
        elif self.method == "svm":
            self.svm = svm.LinearSVC(C=self.C, tol=1e-4, max_iter=1000000)
            self.svm.fit(X, y)
        elif self.method == "svc":
            self.svc = svm.SVC(C=self.C, gamma="scale")
            self.svc.fit(X, y)
        elif self.method == "knn":
            self.knn = KNeighborsClassifier(n_neighbors=25,n_jobs = -1, leaf_size = 50)
            self.knn.fit(X, y)

class ML():
    def __init__(self, features, test, method="load_model", train = None, OutputDir = "output"):
        self.features = features.copy()
        self.train = train
        self.test = test
        self.method = method
        self.train_data = None
        self.load_default_data()
        self.OutputDir = OutputDir
        if not os.path.exists(self.OutputDir):
            os.makedirs(self.OutputDir)
        return

    def process_data(self, data):
        if type(self.features) == list:
            self.return_features = [x + "_normalized" for x in self.features]
            for feature in self.features:
                return_feature = feature + "_normalized"
                if feature == "peak_length":
                    peak_length_min = data[feature].min()
                    peak_length_max = data[feature].max()
                    data[return_feature] = (data[[feature]].values - peak_length_min) / (peak_length_max - peak_length_min)
                elif "percentage" in feature:
                    data[return_feature] = data[[feature]].values
                else:
                    data[return_feature] = preprocessing.QuantileTransformer().fit_transform(data[[feature]].values)
        else:
            print("feature must be list type")
            q()
        return data

    def load_default_data(self):
        self.test_data = self.load_data(file = self.test)
        if self.train:
            self.train_data = self.load_data(file = self.train)


    def load_data(self, file):
        data = pd.read_csv(file)
        if not data.filter(regex="^isin").columns.empty:
            colname = data.filter(regex="^isin").columns[0]
            data["label"] = data[colname]
#             data = data.loc[data[colname] == data["annotated_peak"]]
        data = self.process_data(data = data)
        return data


    def best_c_value(self, train_data, train_label, CList=[0.1, 1, 10, 100, 1000, 10000]):
        value_list = []
        if len(CList) == 1:
            self.c_value =  CList[0]
        else:
            if train_data.shape[1] == 1:
                train_data = train_data.values.reshape(-1,1)
            for c in CList:
                tmp_value = np.mean(cross_val_score(MajorityVoteClassifier(c, method = self.method), train_data, train_label, cv=cv, scoring='f1'))
                value_list.append(tmp_value)
            if np.array(value_list).argmax() != 6:
                C = CList[np.array(value_list).argmax()]
            else:
                first = np.array(value_list[:-1]) 
                last = np.array(value_list[1:])
                arr = last - first
                tmp = np.where(arr<0.01)[0][0]
                C = CList[tmp]
            self.c_value = C
    
    def training(self, CList = [0.1, 1, 10, 100, 1000, 10000]):
        if not isinstance(self.train_data, pd.core.frame.DataFrame):
            print("use 70% data for traning and 30% data for prediction")
            self.train_data, self.test_data = train_test_split(self.test_data, test_size = 0.3, random_state=42)
        train_data = self.train_data[self.return_features]
        train_label = self.train_data["label"]
        self.best_c_value(train_data = train_data, train_label = train_label, CList = CList)
        self.model = MajorityVoteClassifier(self.c_value, method=self.method)
        self.model.fit(train_data, train_label)
        if self.train != None:
            model_name = os.path.join(self.OutputDir, self.train.split("/")[-1].split(".")[0] + "_" + self.method + "_model.m")
        else:
            model_name = os.path.join(self.OutputDir, self.test.split("/")[-1].split(".")[0] + "_" + self.method + "_model.m")
        joblib.dump(self.model, model_name)
        
    def prediction(self, test=None):
        train_data = self.train_data
        train_data_process = self.train_data[self.return_features]
        train_label = self.train_data["label"]
        if test != None:
            test_data = self.load_data(file = test)
            test_data_process = test_data[self.return_features]
            test_label = test_data["label"]
        else:
            test_data = self.test_data
            test_data_process = self.test_data[self.return_features]
            test_label = self.test_data["label"]
            test = self.test
        train_prediction = self.model.predict(train_data_process)
        test_prediction = self.model.predict(test_data_process)
        train_data["model prediction"] = train_prediction
        test_data["model prediction"] = test_prediction
        train_f1_score = metrics.f1_score(train_label, train_prediction)
        test_f1_score = metrics.f1_score(test_label, test_prediction)
        train_accuracy = metrics.accuracy_score(train_label, train_prediction)
        test_accuracy = metrics.accuracy_score(test_label, test_prediction)
        train_recall = metrics.recall_score(train_label, train_prediction)
        test_recall = metrics.recall_score(test_label, test_prediction)
        train_precision = metrics.precision_score(train_label, train_prediction)
        test_precision = metrics.precision_score(test_label, test_prediction)
        train_roc = metrics.roc_auc_score(train_label, train_prediction)
        test_roc = metrics.roc_auc_score(test_label, test_prediction)
        if self.train == None and test == None:
            output = pd.concat([train_data, test_data])
        else:
            output = test_data
        output_name = os.path.join(self.OutputDir, test.split("/")[-1].split(".")[0] + "_" + self.method + "_prediction.csv")
        output.to_csv(output_name, index = False)
        return [train_f1_score, test_f1_score, train_accuracy, test_accuracy, train_recall, test_recall, train_precision, test_precision, train_roc, test_roc]
    
    def gridsearchcv(self, param_grid, **kw):
        if self.method == "rf":
            model = RandomForestClassifier()
        grid_search = GridSearchCV(model, param_grid=param_grid, cv=5, n_jobs=-1, scoring = "accuracy")
        if not isinstance(self.train_data, pd.core.frame.DataFrame):
            data = self.test_data[self.return_features]
            label = self.test_data["label"]
        else:
            data = self.train_data[self.return_features]
            label = self.test_data["label"]
        grid_search.fit(data, label)
        return grid_search

    def model(self, model=None, test=None):
        if model != None:
            try:
                model = joblib.load(model)
            except ModelError:
                print("no model hava found")
        if test != None:
            test_data = self.load_data(file = test)
            test_data_process = test_data[self.return_features]
            test_label = test_data["label"]
        else:
            test_data = self.test_data
            test_data_process = self.test_data[self.return_features]
            test_label = self.test_data["label"]
            test = self.test
        test_prediction = model.predict(test_data_process)
        test_data["model prediction"] = test_prediction
        test_f1_score = metrics.f1_score(test_label, test_prediction)
        test_accuracy = metrics.accuracy_score(test_label, test_prediction)
        test_recall = metrics.recall_score(test_label, test_prediction)
        test_precision = metrics.precision_score(test_label, test_prediction)
        test_roc = metrics.roc_auc_score(test_label, test_prediction)
        output = test_data
        output_name = os.path.join(self.OutputDir, test.split("/")[-1].split(".")[0] + "_" + self.method + "_prediction.csv")
        output.to_csv(output_name, index = False)
        return [test_f1_score, test_accuracy, test_recall, test_precision, test_roc]
    

l=0
dicts = {}
method = "model_load_define"
features = np.append(features_basic, np.append(features_motif, features_internal)).tolist()
ml = ML(features=features, train =train, test =output_file_gene, OutputDir=outputdir)
results = ml.model(model=args.model)
results_intergenic = ml.model(model=args.model, test=output_file_intergenic)
dicts[l] = results + [method, tss] 

df = pd.DataFrame.from_dict(dicts, orient="index", columns=["test_f1_score", "test_accuracy", "test_recall", "test_precision", "test_roc", "method",   "tss"])
df = pd.melt(df, id_vars=["test_accuracy", "test_recall", "test_precision", "test_roc", "method", "tss"])
p = ggplot(df, aes(x="tss", y="value", fill = "method")) + geom_bar(stat="identity", position="dodge") + facet_grid([".","variable"]) + theme_light() + coord_cartesian(ylim=[0,1])
ggsave(p, filename=os.path.join(outputdir, "all.png"), dpi=300)
df = pd.DataFrame.from_dict(dicts, orient="index", columns=["test_f1_score", "test_accuracy", "test_recall", "test_precision", "test_roc", "method",   "tss"])
df.sort_values(["method"])

if not os.path.exists(outputdir):
    os.makedirs(outputdir)
df = pd.DataFrame.from_dict(dicts, orient="index", columns=["test_f1_score", "test_accuracy", "test_recall", "test_precision", "test_roc", "method",   "tss"])
df.to_csv(os.path.join(outputdir, "f1_score.csv"), index=False)
df = pd.melt(df, id_vars=["test_accuracy", "test_recall", "test_precision", "test_roc", "method", "tss"])
p=ggplot(df, aes(x="type", y="value", fill = "method")) + geom_bar(stat="identity", position="dodge") + facet_grid([".","variable"]) + theme_light() + coord_cartesian(ylim=[0.7,1])
ggsave(p, filename=os.path.join(outputdir, "f1_score.png"), dpi=300)
df = pd.DataFrame.from_dict(dicts, orient="index", columns=["test_f1_score", "test_accuracy", "test_recall", "test_precision", "test_roc", "method",   "tss"])
df = pd.melt(df, id_vars=["test_f1_score", "test_recall", "test_precision", "test_roc", "method", "tss"])
p=ggplot(df, aes(x="type", y="value", fill = "method")) + geom_bar(stat="identity", position="dodge") + facet_grid([".","variable"]) + theme_light() + coord_cartesian(ylim=[0.7,1])
ggsave(p, filename=os.path.join(outputdir, "accuracy_score.png"), dpi=300)


