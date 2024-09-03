from aux_functions import *
from plotting_functions import *

from sklearn.linear_model import LinearRegression as LR
def generate_activation_score_df(deseq_effect_mat, all_intra_megaloops, tcon_h3k27ac, 
                                tcon_h3k27me3, cd4_h3k9me3, treg_h3k4me3, treg_h3k4me1, distance_cutoff = np.inf):
    rows = []
    for ind in np.where(all_intra_megaloops.sum(axis=0) > 4)[0]:
        indsoi = np.where(all_intra_megaloops[ind]>0)[0]

        indsoi = indsoi[np.abs(indsoi - ind) < distance_cutoff]

        x = deseq_effect_mat[ind][indsoi]
        bad = np.isnan(x)
        x = x[~bad]
        if len(x) < 8:
            continue
        y = tcon_h3k27ac[indsoi][~bad]
        lr_ac = LR(fit_intercept=True).fit(x.reshape(-1, 1), y).coef_[0]
        r_ac = scipy.stats.pearsonr(x, y)[0]
        
        y = tcon_h3k27me3[indsoi][~bad]
        lr_me3 = LR(fit_intercept=True).fit(x.reshape(-1, 1), y).coef_[0]
        r_me3 = scipy.stats.pearsonr(x, y)[0]
    
        y = cd4_h3k9me3[indsoi][~bad]
        lr_9me3 = LR(fit_intercept=True).fit(x.reshape(-1, 1), y).coef_[0]
        r_9me3 = scipy.stats.pearsonr(x, y)[0]
    
        y = treg_h3k4me3[indsoi][~bad]
        r_h3k4me3 = scipy.stats.pearsonr(x, y)[0]
        
        y = treg_h3k4me1[indsoi][~bad]
        r_h3k4me1 = scipy.stats.pearsonr(x, y)[0]
    
        rows.append([ind, r_ac, r_me3, lr_ac, lr_me3, lr_9me3, r_9me3, r_h3k4me3, r_h3k4me1, len(x)])
    
    megaloop_activation_score_df = pd.DataFrame(rows, columns = ['Ind', 'H3K27ac', 'H3K27me3', 
                                                                 'beta: H3K27ac', 'beta: H3K27me3', 
                                                                 'beta: lr_9me3', 'r_9me3', 'H3K4me3', 'H3K4me1', 'n']).set_index('Ind')
    return megaloop_activation_score_df