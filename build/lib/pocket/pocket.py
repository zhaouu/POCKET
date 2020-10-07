import pandas as pd
import numpy as np
import os
import re
import ipdb
from pandas_plink import read_plink
#from genes_function_svm import gene_function_svm
from imp import reload
from .genes_function_svm import gene_function_svm
from .eqtl import limix_gwas
from os.path import expanduser
class pocket():
    """
    """
    def __init__(self, leadSNP, chrom, region_start, region_end, gwas_res, pheno, bed_f, gene_annotation_df, vars_annotation_df, kinship, positive_sets_dict, feature_df, expression_df_dict=None, twas_df_dict=None, repeat_time=10, n_jobs=10):
        self.leadSNP = leadSNP
        self.chrom = chrom
        self.region_start = region_start
        self.region_end = region_end
        gwas_res = gwas_res.query("(Chr == '%s') & (ChrPos <= %s) & (ChrPos >= %s)" %(chrom, region_end, region_start))
        gwas_res = gwas_res.sort_values(['PValue'])
        #gwas_res.index = gwas_res.SNP
        self.gwas_res = gwas_res
        anno_df = gene_annotation_df.query("chrom == '%s' & start >= %s & end <= %s"  %(chrom, region_start, region_end))
        self.anno_df = anno_df
        self.region_genes = list(anno_df.index)
        self.vars_anno = vars_annotation_df.query("Var_chrom == '%s' & Var_pos >= %s & Var_pos <= %s"  %(chrom, region_start, region_end))
        self.pheno = pheno
        self.bed_f = bed_f
        self.kinship = kinship
        self.positive_sets_dict = positive_sets_dict
        self.feature_df = feature_df
        self.expression_df_dict = expression_df_dict
        self.repeat_time = repeat_time
        self.n_jobs = n_jobs
        self.twas_df_dict = twas_df_dict
        self.hap_res = None
        self.expression_effect =None
        self.region_geno_mat = None
        self.region_snp_info_df = None
        self.snp_info_df = None
        self.ld_res = None
        self.gene_effect_res = None
        self.gf_score = None



    def make_geno_df_region(self, chrom, start, end):
        from limix.io import plink
        (bim, fam, bed) = plink.read(self.bed_f, verbose=False)
        bim.columns = ['chrom', 'rsid', 'cm', 'position', 'a0', 'a1', 'i']
        d = bim.query(" chrom == '{0}' & position >={1}  & position <= {2}".format(chrom, start, end))
        geno_mat = bed[d.i.values, :].compute()
        geno_mat = pd.DataFrame(geno_mat, index = d.rsid, columns = fam.index )
        self.region_geno_mat = geno_mat
        self.region_snp_info_df = d

    def make_geno_df_vars(self, var_ids):
        from limix.io import plink
        (bim, fam, bed) = plink.read(self.bed_f, verbose=False)
        bim.columns = ['chrom', 'rsid', 'cm', 'position', 'a0', 'a1', 'i']
        d = bim.query("rsid in @var_ids")
        geno_mat = bed[d.i.values, :].compute()
        geno_mat = pd.DataFrame(geno_mat, index = d.rsid, columns = fam.index )
        self.geno_mat = geno_mat
        self.snp_info_df = d

    def ld_caculation(self, plink_path=None, temp_path=None, interval =1000):
        """
        """
        import subprocess
        import random
        import glob
        import os
        SNP =self.leadSNP
        if temp_path is None:
            temp_path = expanduser("~")+'/temp'
        if not os.path.isdir(temp_path):
            os.mkdir(temp_path)
        rand_n = str(random.random()*1e7)[0:6]
        if plink_path is None:
            subprocess.call("plink --bfile %s --r2 --ld-snp %s --allow-extra-chr --ld-window-kb %d --ld-window 100000 --ld-window-r2 0  --out %s/%s_ld_res" %(self.bed_f, SNP, interval,temp_path, rand_n), shell=True)
        else:
            subprocess.call("%s --bfile %s --r2 --ld-snp %s --allow-extra-chr --ld-window-kb %d --ld-window 100000 --ld-window-r2 0  --out %s/%s_ld_res" %(plink_path, self.bed_f, SNP, interval,temp_path, rand_n), shell=True)
        ld_res = pd.read_table('%s/%s_ld_res.ld' %(temp_path,rand_n),sep='\s+')
        ld_res.index = ld_res.loc[:,'SNP_B']
        ld_res = ld_res[~ld_res.index.duplicated(keep='first')]
        files = glob.glob('%s/%s_ld_res*' %(temp_path,rand_n))
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass
        self.ld_res = ld_res

    def ev_caculation(self, ld_res=None, plink_path=None, temp_path=None, region_geno=None, ld_thresh= 0.4, effect_w_dict={'high':5,'moderate':3,'modifier':1,'low':1}):
        lmm_df = self.gwas_res
        gene_df = self.anno_df
        effect_anno = self.vars_anno
        i = (effect_anno.Annotation == 'upstream_gene_variant')|(effect_anno.Annotation == 'intragenic_variant')| (effect_anno.Annotation == 'intergenic_region')|(effect_anno.Annotation == 'downstream_gene_variant')|(effect_anno.Annotation == 'intron_variant')|(effect_anno.Annotation == 'N') ##remove non cds variant
        effect_anno = effect_anno[~i]
        if region_geno is None:
            self.make_geno_df_region(self.chrom, self.region_start, self.region_end)
            region_geno = self.region_snp_info_df
        geno_used = region_geno
        geno_used.index = region_geno.rsid
        effect_anno =  effect_anno[effect_anno.index.isin(geno_used.index)]
        snps = [set(g.loc[['Ref_geno','Alt_geno']]) == set(geno_used.loc[i, ['a0','a1']]) for i,g in effect_anno.iterrows()]
        effect_anno = effect_anno[snps]
        if ld_res is None:
            self.ld_caculation()
            ld_res = self.ld_res
        lmm_df.index = lmm_df.SNP
        lmm_df.index.name = None
        lmm_df.loc[:,'r2'] = ld_res.R2
        #print(lmm_df.head())
        #print(effect_anno.head())
        #ipdb.sset_trace()
        anno_res = pd.merge(lmm_df, effect_anno, on='SNP',how='inner')
        anno_res = anno_res.query("r2 > {}".format(ld_thresh))
        if anno_res.shape[0] == 0:
            gene_effect_res = pd.Series(0,index = gene_df.gene)
        else:
            anno_res.loc[:,'scaled_p'] = (-np.log10(anno_res.PValue) -min(-np.log10(anno_res.PValue)))/(max(-np.log10(anno_res.PValue)) - min(-np.log10(anno_res.PValue)))
            score_list = []
            g_list = []
            for g, res in anno_res.groupby('Gene'):
                score = max([x.scaled_p*effect_w_dict[x.Imapct.lower()] for i,x in res.iterrows()])
                score_list.append(score)
                g_list.append(g)
            gene_effect_res = pd.Series(score_list,index=g_list)
            gene_effect_res = (gene_effect_res - min(gene_effect_res))/(max(gene_effect_res) - min(gene_effect_res))
        self.gene_effect_res = gene_effect_res

    def gf_caculation(self):
        prob_list = []
        repeat_time = self.repeat_time
        n_jobs = self.n_jobs
        feature_df = self.feature_df
        positive_sets_dict = self.positive_sets_dict
        for k, posi_genes in positive_sets_dict.items():
            feature_df_used = feature_df.loc[:,~feature_df.columns.str.lower().str.contains(k)]
            svm_res = gene_function_svm(feature_df_used, posi_genes, repeat_time, n_jobs)
            svm_res.mutiple_svm(self.region_genes)
            c, prob, score = svm_res.cat_res, svm_res.prob_res, svm_res.r2_list
            prob_list.append(prob)
        prob_df = pd.DataFrame(prob_list, index = positive_sets_dict.keys()).T
        mean_score = prob_df.mean(axis=1)
        self.gf_score = mean_score

    def topN_eQTL_asso(self, geno_mat, gene_express, snp_info_df):
        res = limix_gwas(geno_mat, gene_express, snp_info_df, kinship=self.kinship)
        #ipdb.sset_trace()
        res.data_check()
        res.maf_filter()
        res.mean_impute()
        g = res.do_gwas()
        g = g.sort_values(['pv20'])
        return g.iloc[0,:]


    def ee_caculation(self,  ld_res =None, p_thresh=1e-5, ld_thresh=0.4, combine='mean'):
        from joblib import Parallel, delayed
        lmm_df = self.gwas_res
        twas_res_dict = self.twas_df_dict
        lmm_df.index = lmm_df.SNP
        n_job = self.n_jobs
        if ld_res is None:
            ld_res = self.ld_res
        lmm_df.loc[:,'r2'] = ld_res.R2
        lmm_df_eqtl = lmm_df.query("r2 > {}".format(ld_thresh))
        i = lmm_df_eqtl.PValue < p_thresh
        if sum(i) <=10:
            top_SNPs = lmm_df_eqtl.SNP[i]
        else:
            top_SNPs = lmm_df_eqtl.SNP[:10]
        if len(top_SNPs) < 2:
            top_SNPs = lmm_df_eqtl.SNP[:3]
        self.make_geno_df_vars(top_SNPs)
        #ipdb.sset_trace()
        twas_eqtl_res_df = pd.DataFrame()
        twas_eqtl_res_df.loc[:,'Gene'] = self.region_genes
        twas_eqtl_res_df.index = self.region_genes
        for k,express_df in self.expression_df_dict.items():
            express_df = express_df[express_df.index.isin(self.region_genes)]
            #for i,g in express_df.iterrows():
            #    self.topN_eQTL_asso(self.geno_mat, g, self.snp_info_df)
            eqtl_res = Parallel(n_jobs= n_job)(delayed(self.topN_eQTL_asso)(self.geno_mat, g, self.snp_info_df) for i,g in express_df.iterrows())
            eqtl_res = pd.DataFrame(eqtl_res)
            eqtl_res.index = express_df.index
            eqtl_res.loc[:,'Gene'] = express_df.index
            if twas_res_dict is not None:
                #eqtl_res =  pd.merge(eqtl_res, self.anno_df, left_on='Gene',right_on='gene',how = 'inner')
                twas_eqtl_res = pd.merge(twas_res_dict[k], eqtl_res, on='Gene', how='inner')
                twas_eqtl_res.loc[:,'twas_p'] = -np.log10(twas_eqtl_res.lmm_p)
                twas_eqtl_res.loc[:,'eqtl_p'] = -np.log10(twas_eqtl_res.pv20)
                region_genes = twas_eqtl_res.Gene
                twas_eqtl_res = twas_eqtl_res.loc[:,'twas_p']*twas_eqtl_res.loc[:,'eqtl_p']
                twas_eqtl_res.index = region_genes
                twas_eqtl_res_df.loc[:,k] = twas_eqtl_res
                twas_eqtl_res.index = region_genes
            else:
                twas_eqtl_res = -np.log10(eqtl_res.pv20)
            twas_eqtl_res_df.loc[:,k] = twas_eqtl_res
        twas_eqtl_res_df = twas_eqtl_res_df.fillna(0)
        if combine == 'max':
            expression_effect = twas_eqtl_res_df.loc[:, list(self.expression_df_dict.keys())]
            expression_effect = expression_effect.max(axis=1)
        else:
            expression_effect = twas_eqtl_res_df.loc[:, list(self.expression_df_dict.keys())]
            expression_effect = expression_effect.mean(axis=1)
        self.expression_effect = (expression_effect - min(expression_effect))/(max(expression_effect) - min(expression_effect))

    def make_gene_hap_pos(self, promoter_len=2000):
        gene_df = self.anno_df
        position_list = []
        for i,g in gene_df.iterrows():
            if g.strand == '-':
                position_list.append([g.start -100, g.end + promoter_len])
            else:
                position_list.append([g.start- promoter_len, g.end + 100])
        return position_list

    def make_gene_hap_df(self, geno_df, start, end, return_min=True, n_clusters=[2,3,4,5,6], min_num=10):
        from scipy import cluster
        from limix.qc import mean_impute
        d = self.region_snp_info_df.query("position >= @start & position <= @end")
        snps = d.rsid
        geno_df = geno_df[geno_df.index.isin(snps)]
        geno = mean_impute(geno_df.values)
        geno_df.index.name =None
        if geno.shape[0] < 5:
            return 1
        Z = cluster.hierarchy.linkage(geno.T)
        cutree = pd.DataFrame(cluster.hierarchy.cut_tree(Z, n_clusters=n_clusters),index= geno_df.columns, columns=n_clusters)
        hap_count = cutree.apply(pd.value_counts).fillna(0)
        hap_list =[]
        i_list = []
        for i in cutree.columns:
            used_h = np.where(hap_count.loc[:,i] >= min_num)[0]
            #ipdb.sset_trace()
            if len(used_h) >1:
                h = np.array(cutree.loc[:,i]).astype('float')
                #ipdb.sset_trace()
                h[~np.isin(h, used_h)[0]] = np.nan
                hap_list.append(h)
                i_list.append(i)
            else:
                continue
        hap_df = pd.DataFrame(hap_list,index=i_list, columns=geno_df.columns)
        if len(hap_list) < 1:
            return 1
        else:
            if len(hap_df.shape) ==1:
                hap_df = pd.DataFrame([hap_df,hap_df], index=['h1','h2'])
            elif hap_df.shape[0] ==1:
                hap_df = hap_df.append(hap_df)
                hap_df.index = ['h1','h2']
            snp_info_df = pd.DataFrame()
            snp_info_df.loc[:,'rsid'] = range(hap_df.shape[0])
            snp_info_df.loc[:,'chrom'] = 1
            snp_info_df.loc[:,'position'] = 1
            #print(hap_df)
            res = limix_gwas(hap_df, self.pheno, snp_info_df, kinship=self.kinship)
            res.data_check()
            res.mean_impute()
            g = res.do_gwas()
            g = g.sort_values(['pv20'])

        if return_min:
            return min(g.pv20)
        else:
            return res

    def eh_caculation(self, region_geno=None, promoter_len=2000, n_clusters=[2,3,4,5,6]):
        from joblib import Parallel, delayed
        n_job= self.n_jobs
        if region_geno is None:
            self.make_geno_df_region(self.chrom, self.region_start, self.region_end)
            geno_df = self.region_geno_mat
        else:
            geno_df = region_geno
        position_list = self.make_gene_hap_pos(promoter_len = promoter_len)
        #pos = position_list[0]
        #p_list = []
        #for gene, pos in zip(self.anno_df.gene, position_list):
        #    print(gene)
        #    p_list.append(self.make_gene_hap_df(geno_df, pos[0], pos[1], n_clusters= n_clusters))

        #xx = self.make_gene_hap_df(geno_df, pos[0], pos[1], n_clusters= n_clusters)
        #ipdb.sset_trace()
        p_list = Parallel(n_jobs= n_job)(delayed(self.make_gene_hap_df)(geno_df, pos[0], pos[1], n_clusters= n_clusters) for pos in position_list)
        hap_p = pd.Series(p_list, index = self.anno_df.gene)
        self.hap_res = (-np.log10(hap_p) - min(-np.log10(hap_p)))/(max(-np.log10(hap_p))- min(-np.log10(hap_p)))


    def pockt_summary(self, save_path = None):
        res_df = pd.DataFrame()
        res_df.loc[:,"Gene"] = self.anno_df.gene
        res_df.loc[:,'symbol'] = self.anno_df.loc[:,'symbol']
        if self.gene_effect_res is None:
            self.ev_caculation()
        if self.gf_score is None:
            self.gf_caculation()
        if self.expression_df_dict is not None:
            if self.expression_effect is None:
                self.ee_caculation()
        if self.hap_res is None:
            self.eh_caculation()
        res_df.loc[:,'variation_effect'] = self.gene_effect_res
        res_df.loc[:,'expression_effect'] = self.expression_effect
        res_df.loc[:,'haplotype_effect'] = self.hap_res
        res_df.loc[:,'gene_function'] = self.gf_score
        if res_df.expression_effect is None:
            res = res_df.gene_function*(res_df.variation_effect  + 2*res_df.haplotype_effect)
        else:
            res = res_df.gene_function*(res_df.variation_effect + res_df.expression_effect + 2*res_df.haplotype_effect)
        res_df.loc[:,'summary_score'] = res
        res_df = res_df.sort_values(['summary_score'],ascending=False)
        res_df.loc[:,'Description'] = self.anno_df.loc[:,'description']
        res_df.loc[:,'Description'] = ['' if pd.isnull(x) else x.split('[')[0] for x in res_df.loc[:,'Description']]
        self.summary_score = res_df
        if save_path is not None:
            res_df.to_csv('{0}/{1}_{2}_{3}_{4}_gene_prori_res.csv'.format(save_path, self.leadSNP, self.chrom, self.region_start, self.region_end))




