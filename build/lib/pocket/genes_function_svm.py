import random
import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import RobustScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
class gene_function_svm():
    def __init__(self, feature_df, good_genes, repeat_time, n_jobs):
        self.feature_df = feature_df
        self.good_genes = good_genes
        self.repeat_time = repeat_time
        self.n_jobs = n_jobs

    def gene_svm(self, test_gene_list, test_chrom=None, negative_genes=None):
        good_genes = set(self.good_genes) & set(self.feature_df.index)
        feature_df = self.feature_df
        if 'Gene' in feature_df.columns:
            feature_df = feature_df.drop('Gene',axis=1)
        if test_chrom is None:
            good_genes = set(good_genes) - set(test_gene_list)
            if negative_genes is None:
                negative_genes = random.sample(set(feature_df.index) - set(good_genes) - set(test_gene_list), len(good_genes)*2)
        else:
            good_genes = set(good_genes) - set([x for x in good_genes if test_chrom in x])
            if negative_genes is None:
                negative_genes = random.sample(set(feature_df.index) - set(good_genes) - set([x for x in feature_df.index if test_chrom in x]), len(good_genes)*2)
        used_genes = list(good_genes) + negative_genes
        X = feature_df.loc[used_genes,:]
        Y=[0]*len(good_genes) +[1]*len(negative_genes)
        clf = Pipeline([
            ('feature_selection', SelectFromModel(ExtraTreesClassifier(n_estimators=500,class_weight='balanced', max_features='auto',n_jobs=4))),
            ('classification', SVC(C=1.8, kernel='rbf', gamma='auto',probability=True, shrinking=True,class_weight='balanced',cache_size=10000))
        ])
        clf.fit(X, Y)
        scores = cross_val_score(clf, X, Y, cv=10)
        mean_score = np.mean(scores)
        score = clf.score(X,Y)
        print('Mean accuracy: %0.3f' %mean_score, 'Predict accuracy: %0.3f' %score)
        res = clf.predict(feature_df.loc[test_gene_list,:])
        prob = clf.predict_proba(feature_df.loc[test_gene_list,:])
        prob = pd.DataFrame(prob, index= test_gene_list, columns=['Y','N'])
        return [res, prob, mean_score]

    def mutiple_svm(self, test_gene_list):
        c_list = []
        Y_prob_list = []
        r2_list = []
        from joblib import Parallel, delayed
        print('Start SVM, %s jobs running.' %self.n_jobs)
        res_list = Parallel(n_jobs= self.n_jobs)(delayed(self.gene_svm)(test_gene_list) for i in range(self.repeat_time))
        for c,prob,score in res_list:
            c_list.append(c)
            Y_prob_list.append(prob.Y)
            r2_list.append(score)
        cat_df = pd.DataFrame(c_list).T
        prob_df = pd.DataFrame(Y_prob_list)
        cat_res = [x.value_counts().index[0] for i,x in cat_df.iterrows()]
        prob_res = prob_df.mean(axis=0)
        self.cat_res = cat_res
        self.prob_res = prob_res
        self.r2_list = np.mean(r2_list)


