from aux_functions import *

def process_liftover_output(path='./hg38_megaloop_liftover/human_regions.bed'):
	liftover_df = pd.read_csv(path, sep='\t', header=None, )
	liftover_df.columns=['chrom', 'start', 'end', 'mouse_start', 'number']
	liftover_df['mouse_start'] = liftover_df['mouse_start'].str.replace("0001-", "0000-")
	liftover_df['size'] = liftover_df['end'] - liftover_df['start']
	liftover_df = liftover_df[liftover_df['size'] < 200_000 ]
	liftover_df = liftover_df[liftover_df['size'] > 10_000 ]
	liftover_df = liftover_df.sort_values('size', ascending=False)
	liftover_df = liftover_df.drop_duplicates('mouse_start')
	return liftover_df

def run_full_mouse_conversion(mouse_megaloops, path='./hg38_megaloop_liftover/human_regions.bed'):
	liftover_df = process_liftover_output(path)
	human_megaloops_50kb = convert_mouse_megaloops_to_human(liftover_df, mouse_megaloops)
	return human_megaloops_50kb, liftover_df

def convert_mouse_megaloops_to_human(liftover_df, mouse_megaloops):
	human_megaloops_50kb = []
	starts = set(liftover_df['mouse_start'].values)
	for row in mouse_megaloops:
		row1, row2 = row[:3], row[3:6]
		row1, row2 = map(add_chr_to_anc, [row1, row2])
		row1, row2 = tuple_to_grange(*row1), tuple_to_grange(*row2)
		if (row1 not in starts) or (row2 not in starts):
			continue
		else:
			human_row1 = liftover_df[liftover_df['mouse_start'] == row1]
			human_row2 = liftover_df[liftover_df['mouse_start'] == row2]
			human_row1 = human_row1[['chrom', 'start', 'end']].values[0]
			human_row2 = human_row2[['chrom', 'start', 'end']].values[0]
			human_megaloops_50kb.append(list(human_row1) + list(human_row2) + list(grange_to_tuple(row1)) + list(grange_to_tuple(row2)))
	return human_megaloops_50kb

def filter_to_similar_distances(random_pair_df, human_megaloop_df, n=1):
	random_pair_df['distance'] = (random_pair_df['start1'] - random_pair_df['start2']).abs()
	random_pair_df = random_pair_df[random_pair_df['distance'] < 1.5*1e8]
	distance_matched = []
	for _, i in human_megaloop_df.iterrows():
		js = (random_pair_df['distance'] - i['distance']).abs().sort_values().index[:n]
		for j in js:
			distance_matched.append(random_pair_df.loc[j])
	return pd.DataFrame(distance_matched).drop_duplicates()


def sample_megaloops_from_nonconnected_megaloop_anchors(human_megaloop_df):
	distance_matched_loop_df = []
	for chrom in human_megaloop_df['chrom1'].unique():
		### Create all tuples that are in the dataframe
		subdf = human_megaloop_df[(human_megaloop_df['chrom1'] == chrom)]
		existing_tuples = []
		for i, j in zip(subdf.iloc[:, :3].values, subdf.iloc[:, 3:6].values):
			existing_tuples.append((i, j))
			existing_tuples.append((j, i))
		
		### Create all possible tuples that are not in the dataframe
		random_pairs = set()
		all_anchors = pd.DataFrame(np.concatenate([subdf.iloc[:, :3].values, subdf.iloc[:, 3:6].values], axis=0))
		for u1 in all_anchors.values:
			for u2 in all_anchors.values:
				u1 = list(u1)
				u2 = list(u2)
				L = tuple(list(u1) + list(u2))
				R = tuple(list(u2) + list(u1))
				if (R not in random_pairs) and (L not in random_pairs):
					random_pairs.add(L)

		### Filter all tuples to be ones that share a similar distance distribution
		random_pair_df = pd.DataFrame(random_pairs, columns = human_megaloop_df.columns[:6])
		distance_matched = filter_to_similar_distances(random_pair_df, human_megaloop_df)
		distance_matched_loop_df.append(distance_matched)
	distance_matched_loop_df = pd.concat(distance_matched_loop_df, axis=0).drop_duplicates()
	return distance_matched_loop_df

def sample_megaloops_from_nonconnected_megaloop_anchors_using_anchors(megaloop_inds, megaloop_mat, all_ind_to_region):
	real_megaloops = []
	fake_megaloops = []
	for u1 in megaloop_inds:
		for u2 in megaloop_inds:
			l1, l2 = all_ind_to_region[u1], all_ind_to_region[u2]
			if megaloop_mat[u1, u2]:
				real_megaloops.append(list(l1) + list(l2))
			elif l1[0] == l2[0]:
				fake_megaloops.append(list(l1) + list(l2))
	real_megaloops = pd.DataFrame(real_megaloops, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
	real_megaloops['distance'] = (real_megaloops['start1'] - real_megaloops['start2']).abs()
	fake_megaloops = pd.DataFrame(fake_megaloops, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
	fake_megaloops['distance'] = (fake_megaloops['start1'] - fake_megaloops['start2']).abs()
	distance_matched_fake_megaloops = filter_to_similar_distances(fake_megaloops, real_megaloops, n = 5)
	return real_megaloops, distance_matched_fake_megaloops