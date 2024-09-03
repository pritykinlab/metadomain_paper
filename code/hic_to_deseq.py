import cooler 
import numpy as np
import pandas as pd
import pybedtools as pbt

wsz = 2
shift = 10
res = 5000
hic_reads = pd.DataFrame()

values = {}
shifted_values = {}
places = []
shifted_places = []


not_shifted = []
is_shifted = []

cools = {}
files = ['Tn_rep1_MAPQ30_raw_5000.cool', 'Tn_rep3_MAPQ30_raw_5000.cool', 'Treg_rep2_MAPQ30_raw_5000.cool', 'Tn_rep2_MAPQ30_raw_5000.cool', 'Treg_rep1_MAPQ30_raw_5000.cool', 'Treg_rep3_MAPQ30_raw_5000.cool']
for file in files:
	file_name = "".join(file.split("_")[:2])
	cools[file_name] = cooler.Cooler(file)

for cool_name in cools:
	values[cool_name] = []
	shifted_values[cool_name] = []

import scipy 
peaks_loops = pbt.BedTool("diff_loops.csv")
for i in peaks_loops:
	l1, l2 = i[:3], i[3:6]
	
	l1 = (l1[0], int(l1[1])-res*wsz, int(l1[2])+res*wsz)
	l2 = (l2[0], int(l2[1])-res*wsz, int(l2[2])+res*wsz)
	

	for cool_name in cools:
		cool = cools[cool_name]

		vals = np.asarray(cool.matrix(balance=False).fetch(l1, l2))
		vals[np.isnan(vals)] = 0
		vals = vals+3
		vals = scipy.ndimage.gaussian_filter(vals, sigma=1)
		val = vals.sum()

		try:
			shifted_l1 = (l1[0], int(l1[1])-shift*res-res*wsz, int(l1[2])-shift*res+res*wsz)
			shifted_l2 = (l2[0], int(l2[1])-shift*res-res*wsz, int(l2[2])-shift*res+res*wsz)
			
			shifted_vals = np.asarray(cool.matrix(balance=False).fetch(shifted_l1, shifted_l2))
			shifted_vals[np.isnan(shifted_vals)] = 0
			shifted_val = shifted_val.sum()

		except:
			shifted_l1 = (l1[0], int(l1[1])+shift*res-res*wsz, int(l1[2])+shift*res+res*wsz)
			shifted_l2 = (l2[0], int(l2[1])+shift*res-res*wsz, int(l2[2])+shift*res+res*wsz)


			shifted_vals = np.asarray(cool.matrix(balance=False).fetch(shifted_l1, shifted_l2))
			shifted_vals[np.isnan(shifted_vals)] = 0
			shifted_val = shifted_vals.sum()

		values[cool_name].append(val)
		shifted_values[cool_name].append(shifted_val)

	place = "_".join(i[:6])
	shifted_place = "_".join([str(x) for x in shifted_l1]) + "_" + "_".join([str(x) for x in shifted_l2])

	places.append(place)
	shifted_places.append(shifted_place)
	not_shifted.append(0)
	is_shifted.append(1)


for cool_name in cools:
	hic_reads[cool_name] = (values[cool_name] + shifted_values[cool_name])

hic_reads["shifted"] = not_shifted + is_shifted
hic_reads['places'] = places + shifted_places

hic_reads.to_csv('hic_values.csv')

