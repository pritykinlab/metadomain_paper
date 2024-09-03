import numpy as np
import scipy
import time
import cooler
from get_focal_contacts import *

import logging
import sys
# Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Set the desired log level
# Create a stream handler for logging to stdout
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)  # Set the desired log level for stream output
# Create a formatter and add it to the stream handler
formatter = logging.Formatter('%(levelname)s - %(message)s')
stream_handler.setFormatter(formatter)
# Add the stream handler to the logger
logger.addHandler(stream_handler)

def prepare_inputs(cool, chrom1, chrom2, type='intra', minweight = .0025, make_assertion=True):
    bal_intra = cool.matrix().fetch(chrom1, chrom2)
    raw_intra = cool.matrix(balance=False).fetch(chrom1, chrom2)

    w1 = cool.bins().fetch(chrom1)['weight']
    w2 = cool.bins().fetch(chrom2)['weight']
    nanind1 = np.zeros_like(w1)
    nanind2 = np.zeros_like(w2)
    nanind1[w1 <= minweight] = np.nan
    nanind2[w2 <= minweight] = np.nan

    nanind1[w1 > .0025] = 0
    nanind2[w2 > .0025] = 0
    extra_nan_filter = np.outer(nanind1, nanind2)
    if make_assertion:
        assert np.isnan(nanind1).mean() < 1e-3, np.isnan(nanind1).mean()
        assert np.isnan(nanind2).mean() < 1e-3, np.isnan(nanind2).mean()

    bal_intra[np.isnan(extra_nan_filter)] = np.nan

    if type == 'intra':
        assert (raw_intra == raw_intra.T).all()
        oe, exp = make_obs_exp_nolog(bal_intra, mat_type=type)	
    elif type == 'inter':
        assert (raw_intra.shape != raw_intra.T.shape)
        oe, exp = make_obs_exp_nolog(bal_intra, mat_type=type)	
    return bal_intra, raw_intra, oe




def compute_peaks(z, useSigma=True, sigma=.75, prominence=5):
    z = z.copy()
    z[np.isnan(z)] = np.nanmedian(z)
    if useSigma:
        z = scipy.ndimage.gaussian_filter(z, sigma=sigma)
    peak_X = np.zeros_like(z)
    peak_Y = np.zeros_like(z)
    for c, row in enumerate(z):
        peaks = scipy.signal.find_peaks(row, prominence=prominence)[0]
        peak_X[c, peaks] = 1

    for c, row in enumerate(z.T):
        peaks = scipy.signal.find_peaks(row, prominence=prominence)[0]
        peak_Y[peaks, c] = 1    
    return z, peak_X, peak_Y

def compute_prominent_peaks(oe,  useSigma=False, sigma=.75, prominence=4,):
    z, peak_X1, peak_Y1 = compute_peaks(oe, useSigma=useSigma, sigma=sigma, prominence=prominence)
    ker = np.ones((3, 3))
    peak_smooth_X1 = scipy.ndimage.correlate(peak_X1, ker)>0
    peak_smooth_Y1 = scipy.ndimage.correlate(peak_Y1, ker)>0
    return peak_smooth_X1, peak_smooth_Y1, z


def compute_interchromosomal_hub_pileup(pileup_treg, pileup_tcon, ind, indsoi, prom = .3, sigma=2):
    n = pileup_treg[ind][indsoi].shape[1]
    midpt = n//2

    delta = pileup_treg[ind][indsoi] - pileup_tcon[ind][indsoi]
    meandelta = np.nanmean(delta, axis=0)
    _1, X_treg, Y_treg = compute_peaks(meandelta, useSigma=True, sigma=sigma, prominence=prom)
    _2, X_tcon, Y_tcon = compute_peaks(-meandelta, useSigma=True, sigma=sigma, prominence=prom)
    row_treg, col_treg  = np.where(X_treg+Y_treg > 1)
    row_tcon, col_tcon  = np.where(X_tcon+Y_tcon > 1)
    
    # Check if centered
    treg_peak = False
    tcon_peak = False
    treg_peak_above_zero = (_1[row_treg, col_treg] > 0)
    tcon_peak_above_zero = (_2[row_tcon, col_tcon] > 0)

    treg_is_mid = (abs(col_treg - midpt) < 3) & (abs(row_treg - midpt) < 3)
    tcon_is_mid = (abs(col_tcon - midpt) < 3) & (abs(row_tcon - midpt) < 3)

    if (treg_is_mid & treg_peak_above_zero).any():
        treg_peak = True
    elif (tcon_is_mid & tcon_peak_above_zero).any():
        tcon_peak = True
    print("__")
    print(tcon_is_mid)
    print(treg_is_mid)

    col_treg, row_treg = col_treg[treg_peak_above_zero], row_treg[treg_peak_above_zero]
    col_tcon, row_tcon = col_tcon[tcon_peak_above_zero], row_tcon[tcon_peak_above_zero]
    return treg_peak, tcon_peak, (col_treg, row_treg), (col_tcon, row_tcon), {'delta_mat' : meandelta, 'filtered_mat' : _2,
                                                                            'treg_mat' : np.nanmean(pileup_treg[ind][indsoi], axis=0), 
                                                                            'tcon_mat' : np.nanmean(pileup_tcon[ind][indsoi], axis=0),
                                                                            }



def compute_interchromosomal_hub_pileup_baseline(pileup_treg, pileup_tcon, ind, indsoi, prom = 3, sigma=2, dist_co=3):
    n = pileup_treg[ind][indsoi].shape[1]
    midpt = n//2

    delta = pileup_treg[ind][indsoi] - pileup_tcon[ind][indsoi]
    meandelta = np.nanmean(delta, axis=0)

    treg_mat = np.nanmean(pileup_treg[ind][indsoi], axis=0)
    tcon_mat = np.nanmean(pileup_tcon[ind][indsoi], axis=0)
    
    _1, X_treg, Y_treg = compute_peaks(treg_mat, useSigma=True, sigma=sigma, prominence=prom)
    _2, X_tcon, Y_tcon = compute_peaks(tcon_mat, useSigma=True, sigma=sigma, prominence=prom)
    row_treg, col_treg  = np.where(X_treg+Y_treg > 1)
    row_tcon, col_tcon  = np.where(X_tcon+Y_tcon > 1)
    

    # Check if centered
    treg_peak = False
    tcon_peak = False
    treg_peak_above_zero = (_1[row_treg, col_treg] > 0)
    tcon_peak_above_zero = (_2[row_tcon, col_tcon] > 0)
    treg_is_mid = (abs(col_treg - midpt) < dist_co) & (abs(row_treg - midpt) < dist_co)
    tcon_is_mid = (abs(col_tcon - midpt) < dist_co) & (abs(row_tcon - midpt) < dist_co)
    
    # print("___")
    # print(row_tcon, col_tcon, midpt, tcon_peak_above_zero, tcon_is_mid)
    # print((abs(col_tcon - midpt)), (abs(row_tcon - midpt)))

    if (treg_is_mid & treg_peak_above_zero).any():
        treg_peak = True
    if (tcon_is_mid & tcon_peak_above_zero).any():
        tcon_peak = True

    col_treg, row_treg = col_treg[treg_peak_above_zero], row_treg[treg_peak_above_zero]
    col_tcon, row_tcon = col_tcon[tcon_peak_above_zero], row_tcon[tcon_peak_above_zero]
    return treg_peak, tcon_peak, (col_treg, row_treg), (col_tcon, row_tcon), {'delta_mat' : meandelta, 'filtered_mat' : _2,
                                                                            'treg_mat' : treg_mat,
                                                                            'tcon_mat' : tcon_mat,
                                                                            }



def compute_interchromosomal_hub_single_pileup(pileup, ind, indsoi):
    n = pileup[ind][indsoi].shape[1]
    midpt = n//2    
    delta = pileup[ind][indsoi]
    meandelta = np.nanmean(delta, axis=0)
    prom = .3
    _, X_treg, Y_treg = compute_peaks(meandelta, useSigma=True, sigma=2, prominence=prom)
    row_treg, col_treg  = np.where(X_treg+Y_treg > 1)
    
    # Check if centered
    treg_peak = False
    tcon_peak = False
    if ((abs(col_treg - midpt) < 5) & (abs(row_treg - midpt) < 5)).any():
        treg_peak = True
    return treg_peak, (col_treg, row_treg), {'delta_mat' : meandelta}


def get_hiccups_pvalue(raw_intra, bal_intra, type='intra', filter_n=35, filter_width=3):
    raw_intra = raw_intra.copy().astype(float)
    t1 = time.time()
    if type=='inter':
        _, full_logp_mat, _, _, _, _, _ = call_peak_with_poisson_and_peak_cutoff(raw_intra+.5, 
                                                                                 filter_n=filter_n, 
                                                                                 filter_width=filter_width)
    elif type=='intra':
        _, full_logp_mat, _, _, _, _, _ = call_peak_with_poisson_and_peak_cutoff(raw_intra, )
    elif type=='mega':
        _, full_logp_mat, _, _, _, _, _ = call_peak_with_poisson_and_peak_cutoff(raw_intra, filter_n=filter_n, filter_width=filter_width)

    full_logp_mat[np.isnan(bal_intra)] = 0
    return full_logp_mat

def collapse_hiccups_and_peak_pvalues(full_logp_mat, peak_smooth_X1, peak_smooth_Y1, oe, pco = 20):
    megaloops_of_interest = (full_logp_mat>pco) & peak_smooth_X1 & peak_smooth_Y1
    label_mat, _ = scipy.ndimage.label(megaloops_of_interest)
    obj_locations = scipy.ndimage.find_objects(label_mat)
    collapsed_logp_mat = np.zeros(megaloops_of_interest.shape)*2

    for c, i in enumerate(obj_locations):
        l, r = i
        # Get OE matrix
        submat = oe[l, r].copy() * (label_mat[l, r] == (c+1))
        S, E = np.unravel_index(np.nanargmax(np.ravel(submat)), submat.shape)
        s1, s2 = l.start, r.start
        collapsed_logp_mat[s1+S, s2+E] = full_logp_mat[s1+S, s2+E].copy()
    return collapsed_logp_mat

def set_filters_custom(n, w):
    
    inner_filter = np.zeros((n, n))
    outer_filter = np.ones((n, n))

    mid = inner_filter.shape[1]//2
    inner_filter[mid-w:mid+w+1, mid-w:mid+w+1] = 1
    outer_filter -= inner_filter

    left_filter = np.zeros((n, n))
    left_filter[mid-w:mid+w+1, :n//2] = 1
    left_filter[inner_filter==1] = 0

    right_filter = np.zeros((n, n))
    right_filter[mid-w:mid+w+1, n//2:] = 1
    right_filter[inner_filter==1] = 0

    up_filter = left_filter.T
    down_filter = right_filter.T

    side_filter = left_filter + right_filter
    vert_filter = up_filter + down_filter

    return inner_filter, outer_filter, left_filter, right_filter, up_filter, down_filter


import skimage
from skimage.feature import peak_local_max

def call_peak_with_poisson_and_peak_cutoff(image, sigma=1, pmin = 1e-300, logp_co = 10, 
                                            more_than_min_co = 4, filter_n=15, filter_width=1, 
                                            frac_min_valid=0, verbose=False):
    if verbose:
        logger.setLevel(logging.INFO)  # Set the log level to INFO
    else:
        logger.setLevel(logging.WARNING)  # Set the log level to WARNING or any desired level

    C, O, L, R, U, D = set_filters_custom(filter_n, filter_width)
    nanfilter_image = image.copy()
    nanfilter_image[np.isnan(nanfilter_image)] = 0
    filt_image = scipy.ndimage.gaussian_filter(nanfilter_image, sigma=1)
    n_C = C.sum()
    n_O = O.sum()
    
    logger.info("Correlating")
    frac_valid_C = 1-scipy.ndimage.correlate(np.isnan(image).astype(float), C)/n_C
    frac_valid_O = 1-scipy.ndimage.correlate(np.isnan(image).astype(float), O)/n_O

    counts_C = (scipy.ndimage.correlate(nanfilter_image, C)/frac_valid_C).astype(int)
    counts_O = (scipy.ndimage.correlate(nanfilter_image, O)/frac_valid_O).astype(int)
    logger.info("Done correlating")
    
    logger.info("Poisson")
    pval_mat = scipy.stats.poisson(counts_O * n_C/n_O).sf(counts_C)+pmin
    logp_mat = -np.log10(pval_mat)

    frac_valid_bool = (frac_valid_C < frac_min_valid) | (frac_valid_O < frac_min_valid)
    logp_mat[frac_valid_bool] = 0

    is_sig = (logp_mat > logp_co)# & is_peak
#     med = np.mean(filt_image[filt_image > 0])
#     is_peak = is_peak*(filt_image > med*more_than_min_co)
    logger.info("Done poisson")
    logger.info("Peak local max")
#     is_peak = peak_local_max(filt_image, indices=False) 
    logger.info("Done peak local max")
    row, col = np.where(is_sig)
    
    ind_mat = is_sig
    
    no_edge_ind_mat = np.zeros_like(ind_mat)
    sl = slice(10, -10)
    no_edge_ind_mat[sl, sl] = ind_mat[sl, sl]
    peak_row, peak_col = np.where(no_edge_ind_mat)
    
    resultsdict = {}
    resultsdict['counts_inner'] = counts_C
    resultsdict['counts_outer'] = counts_O
    return filt_image, logp_mat, row, col, peak_row, peak_col, resultsdict

