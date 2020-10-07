import os
import numpy as np
import pandas as pd
class limix_gwas():
    def __init__(self, geno_matrix, pheno_list, snp_info_df, kinship=None):
        self.geno_matrix = geno_matrix
        self.geno_matrix.index = snp_info_df.rsid
        self.pheno_list = pheno_list.dropna()
        self.kinship = kinship
        self.SNPinfo = snp_info_df

    def data_check(self):
        if self.kinship is not None:
            culs = set(self.geno_matrix.columns) & set(self.pheno_list.index) & set(self.kinship.index)
            self.kinship = self.kinship.loc[culs,culs]
        else:
            culs = set(self.geno_matrix.columns) & set(self.pheno_list.index)
        self.geno_matrix = self.geno_matrix.loc[:,culs]
        self.pheno_list = self.pheno_list.loc[culs]

    def maf_filter(self, thresh =0.05):
        from limix.qc import compute_maf
        print("Start SNP filtering....(MAF > %s)" %thresh)
        maf_list = compute_maf(self.geno_matrix.T.values)
        self.geno_matrix = self.geno_matrix[np.array(maf_list) > thresh]
        self.SNPinfo = self.SNPinfo[np.array(maf_list) > thresh]
        print('%s variations used to perform GWAS.' %self.SNPinfo.shape[0])

    def mean_impute(self):
        from limix.qc import mean_impute
        print("Genotype imputation with mean value")
        geno_matrix = mean_impute(self.geno_matrix.T.values)
        self.geno_matrix = geno_matrix.T
    def do_gwas(self, geno_mat =None):
        from limix.qtl import scan
        print("Start to perform GWAS......")
        if geno_mat is None:
            geno_mat = self.geno_matrix.T
        else:
            geno_mat = geno_mat
        if self.kinship is not None:
            res = scan(geno_mat, self.pheno_list.values, "normal", K= self.kinship.values ,verbose=False)
        else:
            res = scan(self.geno_matrix.T.values, self.pheno_list.values, "normal", K= None ,verbose=False)
        res_p = res.stats
        res_p.index = self.SNPinfo.rsid
        res_p.loc[:,'rsid'] = self.SNPinfo.rsid
        res_p.loc[:,'chrom'] = self.SNPinfo.chrom
        res_p.loc[:,'position'] = self.SNPinfo.position
        betas = np.array(res.effsizes['h2'].effsize[res.effsizes['h2'].effect_type=='candidate'])
        se = np.array(res.effsizes['h2'].effsize_se[res.effsizes['h2'].effect_type=='candidate'])
        res_p.loc[:,'beta'] = betas
        res_p.loc[:,'se'] = se
        res_p.loc[:,'z_score'] = betas/se
        self.res_p = res_p
        return res_p
    def chrom_gwas_paiallel(self, njobs=10):
        from joblib import Parallel, delayed
        geno_list = []
        for g,d in self.SNPinfo.groupby('chrom'):
            snps = d.rsid
            geno_mat =  self.geno_matrix.loc[snps,:].T.values
            geno_list.append(geno_mat)
        res_list = Parallel(n_jobs= n_jobs)(delayed(do_gwas)(mat) for mat in geno_list)
        res_p = pd.DataFrame()
        for res in res_list:
            res_p = res_p.append(res, ignore_index=True)
        self.res_p = res_p
        return res_p
    def lead_SNP_caculation(self, plink_bed_f, temp_path, r2_thresh=0.2, interval=5000, r2_N_filter=20, r2_N_threshold=0.4, r2_max_filter= 0.8, p_thresh=1e-5, topN = 1000, plink_path=None):
        print("Identify leadSNP ......")
        self.res_p = self.res_p.sort_values(['pv20'])
        if os.path.isdir(temp_path):
            pass
        else:
            os.mkdir(temp_path)
        sig_list = []
        df = self.res_p
        df = df[df.pv20 <= p_thresh]
        sorted_snps = list(df.rsid)
        while len(sorted_snps) > 0:
            rand_n = str(random.random()*1e7)[0:6]
            leadSNP = sorted_snps[0]
            if len(sig_list) > topN:
                break
            try:
                subprocess.call("plink --bfile %s --r2 --ld-snp %s --allow-extra-chr --ld-window-kb %d --ld-window 100000 --ld-window-r2 0  --out %s/%s_ld_res" %(plink_bed_f, leadSNP, interval,temp_path, rand_n), shell=True)
            except IOError:
                rand_n = str(random.random()*1e7)[0:6]
                subprocess.call("plink --bfile %s --r2 --ld-snp %s --allow-extra-chr --ld-window-kb %d --ld-window 100000 --ld-window-r2 0  --out %s/%s_ld_res" %(plink_bed_f, leadSNP, interval, temp_path, rand_n), shell=True)
            ld_res = pd.read_table('%s/%s_ld_res.ld' %(temp_path, rand_n),sep='\s+')
            ld_res.index = ld_res.loc[:,'SNP_B']
            ld_res = ld_res[~ld_res.index.duplicated(keep='first')]
            nlargest_r2 = ld_res.R2.nlargest(r2_N_filter+1)
            if (nlargest_r2.iloc[1] < r2_max_filter) and (nlargest_r2.iloc[-1] < r2_N_threshold):
                droped_SNPs = ld_res[ld_res.R2 >= r2_thresh].SNP_B
                sorted_snps = [x for x in sorted_snps if x not in droped_SNPs]
                continue
            droped_SNPs = ld_res[ld_res.R2 >= r2_thresh].SNP_B
            sorted_snps = [x for x in sorted_snps if x not in droped_SNPs]
            sig_list.append(leadSNP)
        leadSNP_df = df.loc[sig_list,:]
        self.leadSNP_df = leadSNP_df
        return leadSNP_df

    def save_gwas(self, save_f,type='hdf5'):
        print("Save GWAS results ....")
        if type == 'hdf5':
            try:
                self.res_p.to_hdf(save_f, key= 'gwas_df', format='table', data_columns=['rsid','chrom','position','pv20'])
            except NameError:
                pass
            except AttributeError:
                pass

            try:
                self.leadSNP_df.to_hdf(save_f, key= 'leadSNP', format='table', data_columns=['rsid','chrom','position'])
            except NameError:
                pass
            except AttributeError:
                pass

