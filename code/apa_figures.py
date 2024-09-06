import numpy as np


def make_expected(df, chrom, s, e, wsz, res=400, balanced=True):
    chrom = str(chrom)
    diag = (e-s)//res
    expected = np.zeros((wsz*2+1, wsz*2+1))
    n = expected.shape[0]

    test = np.zeros((wsz*2+1, wsz*2+1))
    n = test.shape[0]

    assert type(chrom) == str
    averages = df[(df.region1 == chrom)]
    assert len(averages) > 100, print(chrom)
    diags = np.flip(np.indices((n, n))[0]) + (np.indices((n, n))[1]) - (n-1)
    assert (np.diag(diags) == 0).all()
    diags = diags + diag
    for val in np.unique(diags):
        if balanced==True:
            row = averages[averages['dist'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['balanced.avg'])
        elif balanced==False:
            row = averages[averages['dist'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['count.avg'])
    return expected

def apa_noroll(loops, cool_dict, expected_dict, wsz=20, skip_small=True, balance=True):
    mat_dict = {}
    expecteds = {}
    n=2*wsz+1

    for condition in cool_dict:
        cool = cool_dict[condition]
        res = cool.info['bin-size']
        expected_df = expected_dict[condition]
        mat_dict[condition] = []
        expecteds[condition] = []

        places = []

        for c, i in enumerate(loops):
            l1, l2 = i[:3], i[3:6]
            chrom = l1[0]
            s = ((int(l1[1]) + int(l1[2]))//2)//res*res
            e = ((int(l2[1]) + int(l2[2]))//2)//res*res
            if ((e-s)/res < wsz) and skip_small==True:
                print("Skipping small")
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                continue
            places.append([i[:6]])
            expected = make_expected(expected_df, chrom, s, e, wsz, res=res, balanced=balance).astype(float)
            new_l1 = (l1[0], int(s)-res*wsz, int(s)+res*(wsz+1))

            new_l2 = (l2[0], int(e)-res*wsz, int(e)+res*(wsz+1))
            try:
                val = np.asarray(cool.matrix(balance=balance).fetch(new_l1, new_l2)).astype(float)
            except Exception as e:
                print("Could not fetch in ", condition, ' error: ', e)
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                continue
            if val.shape != (2*wsz+1, 2*wsz+1):
                print(l1, l2)
                print(val.shape)
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                print("Wrong shape")
                continue
            diag = (e-s)//res
            inds = -np.arange(0, 2*wsz+1)
            good_inds = np.subtract.outer(inds, inds)+diag >= 2
            val[~good_inds] = np.nan
            expected[~good_inds] = np.nan

            mat_dict[condition].append(val)
            expecteds[condition].append(expected)
            c += 1

    for key in mat_dict:
        mat_dict[key] = np.asarray(mat_dict[key])
    for key in expecteds:
        expecteds[key] = np.asarray(expecteds[key])
    return mat_dict, expecteds, places



    

def apa_diag(anchors, cool_dict, expected_dict, wsz=20):
    roll_sz = 3
    print(roll_sz)
    mat_dict = {}
    expecteds = {}
    # pc_rat = .5
    pc_rat = .1
    n=2*wsz+1

    for condition in cool_dict:
        cool = cool_dict[condition]
        res = cool.info['bin-size']
        expected_df = expected_dict[condition]
        mat_dict[condition] = []
        expecteds[condition] = []

        places = []
        for c, i in enumerate(anchors):
            l1 = i[:3]
            chrom = l1[0]
            s = ((int(l1[1]) + int(l1[2]))//2)//res*res
            places.append([i[:3]])
            expected = make_expected(expected_df, chrom, s, s, wsz, res=res, balanced=True)
            new_l1 = (l1[0], int(s)-res*wsz, int(s)+res*(wsz+1))
            try:
                val = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l1))
            except Exception as e:
                print("Could not fetch in ", condition, ' error: ', e)
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                continue
            if val.shape != (2*wsz+1, 2*wsz+1):
                print(val.shape)
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                print("Wrong shape")
                continue
            diag = 0

            inds = -np.arange(0, 2*wsz+1)
            good_inds = np.subtract.outer(inds, inds)+diag >= 2

            val[~good_inds] = np.nan
            expected[~good_inds] = np.nan

            mat_dict[condition].append(val)
            expecteds[condition].append(expected)
            c += 1
        print("Done with", condition)

    for key in mat_dict:
        mat_dict[key] = np.asarray(mat_dict[key])
    for key in expecteds:
        expecteds[key] = np.asarray(expecteds[key])
    return mat_dict, expecteds, places




def apa_diag_nonorm(anchors, cool_dict, wsz=20):
    mat_dict = {}

    for condition in cool_dict:
        cool = cool_dict[condition]
        res = cool.info['bin-size']
        mat_dict[condition] = []

        places = []
        for c, i in enumerate(anchors):
            l1 = i[:3]
            s = ((int(l1[1]) + int(l1[2]))//2)//res*res
            places.append([i[:3]])
            new_l1 = (l1[0], int(s)-res*wsz, int(s)+res*(wsz+1))
            try:
                val = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l1))
            except Exception as e:
                print("Could not fetch in ", condition, ' error: ', e)
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                continue
            if val.shape != (2*wsz+1, 2*wsz+1):
                print(val.shape)
                bad = np.zeros((2*wsz+1, 2*wsz+1))*np.nan
                mat_dict[condition].append(bad)
                print("Wrong shape")
                continue
            diag = 0

            inds = -np.arange(0, 2*wsz+1)
            good_inds = np.subtract.outer(inds, inds)+diag >= 2

            val[~good_inds] = np.nan
            mat_dict[condition].append(val)
            c += 1
        print("Done with", condition)

    for key in mat_dict:
        mat_dict[key] = np.asarray(mat_dict[key])
    return mat_dict, places


from aux_functions import make_int
import cooler
def process_anchors(anchors_chunk, cool_dict, wsz, res):
    mat_dict_chunk = {}
    places_chunk = []
    for condition in cool_dict:
        logging.info(f"Processing {condition}")
        cool = cooler.Cooler(cool_dict[condition])
        logging.info(f"Created Cooler object for {condition}")
        mat_dict_chunk[condition] = []
        for i in anchors_chunk:
            l1 = make_int(i[:3])
            s = ((int(l1[1]) + int(l1[2])) // 2) // res * res
            places_chunk.append([i[:3]])
            new_l1 = (l1[0], int(s) - res * wsz, int(s) + res * (wsz + 1))
            try:
                val = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l1))
            except Exception as e:
                logging.error(f"Could not fetch in {condition}, error: {e}")
                bad = np.zeros((2 * wsz + 1, 2 * wsz + 1)) * np.nan
                mat_dict_chunk[condition].append(bad)
                continue
            if val.shape != (2 * wsz + 1, 2 * wsz + 1):
                logging.warning(f"Wrong shape: {val.shape} in {condition}")
                bad = np.zeros((2 * wsz + 1, 2 * wsz + 1)) * np.nan
                mat_dict_chunk[condition].append(bad)
                continue
            diag = 0
            inds = -np.arange(0, 2 * wsz + 1)
            good_inds = np.subtract.outer(inds, inds) + diag >= 2
            val[~good_inds] = np.nan
            mat_dict_chunk[condition].append(val)
    return mat_dict_chunk, places_chunk


import pybedtools as pbt
def split_bedtool(bedtool, num_splits):
    n = len(bedtool)
    chunk_size = n // num_splits
    chunks = []
    for i in range(num_splits):
        if i == num_splits - 1:
            chunks.append(list(bedtool[i * chunk_size :]))
        else:
            chunks.append(list(bedtool[i * chunk_size : (i + 1) * chunk_size]))
    return chunks

import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set the logging level
    format='%(asctime)s - %(levelname)s - %(message)s',  # Log message format
    datefmt='%Y-%m-%d %H:%M:%S',  # Date format
    handlers=[
        logging.FileHandler("apa_diag_nonorm_parallel.log"),  # Log to a file
        logging.StreamHandler()  # Also log to console
    ]
)

import multiprocessing as mp
def apa_diag_nonorm_parallel(anchors, cool_dict, wsz=20, num_processes=30):
    logging.info("Starting apa_diag_nonorm_parallel")
    res = int(next(iter(cool_dict.values())).split("/")[-1])
    logging.info(f"Using bin size: {res}")

    # Split the anchors into chunks for each process
    anchors_chunks = split_bedtool(anchors, num_processes)
    # Create a multiprocessing pool
    with mp.Pool(processes=num_processes) as pool:
        logging.info("Starting multiprocessing pool")
        results = [pool.apply_async(process_anchors, args=(chunk, cool_dict, wsz, res)) for chunk in anchors_chunks]
        output = [p.get() for p in results]
        logging.info("Finished processing all chunks")

    # Combine the results from all processes
    mat_dict_combined = {}
    places_combined = []

    for mat_dict_chunk, places_chunk in output:
        places_combined.extend(places_chunk)
        for condition in mat_dict_chunk:
            if condition not in mat_dict_combined:
                mat_dict_combined[condition] = []
            mat_dict_combined[condition].extend(mat_dict_chunk[condition])

    # Convert lists to arrays
    for key in mat_dict_combined:
        mat_dict_combined[key] = np.asarray(mat_dict_combined[key])

    logging.info("apa_diag_nonorm_parallel completed")
    return mat_dict_combined, places_combined
