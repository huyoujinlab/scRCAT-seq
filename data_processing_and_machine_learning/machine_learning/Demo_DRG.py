import six
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = [15, 15]
import seaborn as sns
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.tree import export_graphviz
import scikitplot as skplt
import os
from mpl_toolkits.mplot3d import Axes3D
from sklearn import datasets
from sklearn.decomposition import PCA, KernelPCA
from sklearn.feature_selection import f_regression, mutual_info_regression
import glob
import math

feature_other_study = ["peak_RPM_normalized", "dominant_RPM_normalized", "peak_length_normalized", "percentage_normalized"]
pd_tss_train        = pd.read_csv('./data/ERCC_tss_peak_new.csv')
pd_tss_test          = pd.read_csv('./data/20190711DRGtes_peak_new.csv')

pd_tes_train     = pd.read_csv('./data/ERCC_tes_peak_new.csv')
pd_tes_test       = pd.read_csv('./data/20190711DRGtss_peak_new_new.csv')

# Get the min and max for each feature in tss

tss_gene_length_min = min(pd_tss_train["gene_length"].min(), pd_tss_test["gene_length"].min())
tss_gene_length_max = max(pd_tss_train["gene_length"].max(), pd_tss_test["gene_length"].max())

tss_peak_length_min = min(pd_tss_train["peak_length"].min(), pd_tss_test["peak_length"].min())
tss_peak_length_max = max(pd_tss_train["peak_length"].max(), pd_tss_test["peak_length"].max())

# Get the min and max for each feature in tes

tes_gene_length_min = min(pd_tes_train["gene_length"].min(), pd_tes_test["gene_length"].min())
tes_gene_length_max = max(pd_tes_train["gene_length"].max(), pd_tes_test["gene_length"].max())

tes_peak_length_min = min(pd_tes_train["peak_length"].min(), pd_tes_test["peak_length"].min())
tes_peak_length_max = max(pd_tes_train["peak_length"].max(), pd_tes_test["peak_length"].max())

def load_other_study(url):
  raw_data = pd.read_csv(url)
  
  raw_data['peak_RPM_normalized']              = preprocessing.QuantileTransformer().fit_transform(raw_data[['peak_RPM']].values)
  
  if "tss_" in url or "cap_" in url: 
    raw_data['dominant_RPM_normalized']      = preprocessing.QuantileTransformer().fit_transform(raw_data[['dominant_tss_RPM']].values)
  else:
    raw_data['dominant_RPM_normalized']      = preprocessing.QuantileTransformer().fit_transform(raw_data[['dominant_tes_RPM']].values)

  if "tss_" in url or "cap_" in url:
    raw_data['gene_length_normalized']           = (raw_data[['gene_length']].values - tss_gene_length_min) / (tss_gene_length_max - tss_gene_length_min) 
  else:
    raw_data['gene_length_normalized']           = (raw_data[['gene_length']].values - tes_gene_length_min) / (tes_gene_length_max - tes_gene_length_min)
  
  if "tss_" in url or "cap_" in url:
    raw_data['peak_length_normalized']        = (raw_data[['peak_length']].values - tss_peak_length_min) / 1620.0
  else:
    raw_data['peak_length_normalized']        = (raw_data[['peak_length']].values - tes_peak_length_min) / 1620.0

  raw_data['percentage_normalized'] = raw_data[['percentage']].values
  
  return raw_data[feature_other_study], raw_data["annotated_peak"]

  ## Use only 4 features

train_tss, train_tss_label       = load_other_study('./data/ERCC_tss_peak_new.csv')
test_tss, test_tss_label           = load_other_study('./data/ERCC_tss_from_DRGdeep_peak_new.csv')

train_tes, train_tes_label       = load_other_study('./data/ERCC_tes_peak_new.csv')
test_tes, test_tes_label           = load_other_study('./data/ERCC_tes_from_DRGdeep_peak_new.csv')

from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import cross_val_score
cv = ShuffleSplit(n_splits=5, test_size=0.3, random_state=0)
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier

from sklearn.base import BaseEstimator

# cross validation to select hyperparameter: damping parameter "C".

class MajorityVoteClassifier(BaseEstimator):
  def __init__(self, C = 1):
    self.C = C
    return

  def predict(self, X):
    y_lr_pred  = self.lr.predict(X)
    y_rf_pred  = self.rf.predict(X)
    y_svm_pred = self.svm.predict(X)
    y_majority_vote = (y_lr_pred & y_rf_pred) | (y_lr_pred & y_svm_pred) | (y_rf_pred & y_svm_pred)
    return y_majority_vote

  def fit(self, X, y, **kwargs):
    self.lr = LogisticRegression(C=self.C, random_state=1, solver="lbfgs", max_iter=1000000)
    self.svm = svm.LinearSVC(C=self.C, tol=1e-4, max_iter=1000000)
    self.rf = RandomForestClassifier(n_estimators=100, random_state=0)
    
    self.lr.fit(X, y)
    self.svm.fit(X, y)
    self.rf.fit(X, y)

print(np.mean(cross_val_score(MajorityVoteClassifier(1000), train_tss, train_tss_label, cv=cv, scoring='f1')))


def result(X_train, y_train, X_test, y_test, test_csv, data="tss"):
  print("-------------models trained on " + data + " dataset-----------")
  
  clc = MajorityVoteClassifier(1000)
  clc.fit(X_train, y_train)
  y_train_prediction = clc.predict(X_train)  

  print("Model performance using f1-score on training data")
  print(metrics.f1_score(y_train, y_train_prediction))

  y_test_prediction = clc.predict(X_test)

  print("Model performance using f1-score on test data")
  print(metrics.f1_score(y_test, y_test_prediction))

  # save predictions to file
  output = test_csv
  output["model prediction"] = y_test_prediction
  output.to_csv("./output/output_DRG_" + data + ".csv")

  print("-------------------------------------------------------------------------------------------------------------------")

if os.path.exists("output"):
    pass
else:
    os.mkdir("output")
result(train_tss, train_tss_label, test_tss, test_tss_label, pd_tss_test, data="tss")
result(train_tes, train_tes_label, test_tes, test_tes_label, pd_tes_test, data="tes")
