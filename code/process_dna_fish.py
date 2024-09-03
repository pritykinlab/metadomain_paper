from plotting_functions import init_subplots
import matplotlib.pyplot as plt
from scipy.ndimage import label
from skimage.feature import peak_local_max
from skimage.measure import find_contours
import numpy as np
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
print("Imported all modules")
import pandas as pd
from skimage import measure
from skimage.filters import difference_of_gaussians as dog

import tifffile as tiff

class FishExperiment():
	def __init__(self, path_to_dapi, path_to_fluorophore_dict, lowco = 800, highco = 3000, threshold_rel = 0.5,
			  segment_co = 40, filter_ctla4=False):
		self.threshold_rel = threshold_rel
		self.lowco = lowco
		self.highco = highco
		self.dapi_img = tiff.imread(path_to_dapi)
		self.dapi_array = np.array(self.dapi_img).sum(axis=2)
		self.flurophore_dict = {}
		for key, value in path_to_fluorophore_dict.items():
			img = tiff.imread(value)
			array = np.array(img).sum(axis=2)
			setattr(self, key + '_array', array)
			self.flurophore_dict[key] = array
			if filter_ctla4:
				if key == 'Ctla4':
					self.flurophore_dict[key] = self.flurophore_dict[key] - self.flurophore_dict[key].min()
					bad_inds = self.flurophore_dict[key] > 40_000
					struct2 = ndi.generate_binary_structure(2, 2)
					dilated_structure = ndi.morphology.binary_dilation(bad_inds, structure=struct2, iterations=30)
					self.flurophore_dict[key][dilated_structure] = 0

		self.cell_labels = self.segment_cells(co=segment_co)
		self.cell_labels = self.permute_labels()
		self.all_original_cells = self.cell_labels.copy()
		self.cell_ids, self.cell_counts = self.set_label_counts()
		self.cell_labels = self.remove_outlier_cells()
		self.peak_dict, self.dogmat_dict = self.call_peaks(threshold_rel=threshold_rel)
		self.peak_cell_df = self.make_peak_cell_df()
		self.valid_cell_ids = np.unique(self.cell_labels[~np.isnan(self.cell_labels)])

	def remove_outlier_cells(self, ):
		small_cells_idx = (self.cell_counts < self.lowco) 
		large_cells_idx = (self.cell_counts > self.highco)
		good_size_idx = (~small_cells_idx) & (~large_cells_idx)
		good_idx = np.where(good_size_idx)[0]
		filtered_array = self.cell_labels.copy()
		filtered_array = filtered_array.copy()
		self.small_cells = self.cell_labels.copy()
		self.large_cells = self.cell_labels.copy()

		filtered_array[~np.isin(filtered_array, self.cell_ids[good_idx])] = np.nan
		self.small_cells[~np.isin(self.cell_labels, self.cell_ids[small_cells_idx])] = np.nan
		self.large_cells[~np.isin(self.cell_labels, self.cell_ids[large_cells_idx])] = np.nan
		return filtered_array

	def segment_cells(self, co=40):
		image = self.dapi_array > co
		distance = ndi.distance_transform_edt(image)
		coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=image, min_distance=10)
		mask = np.zeros(distance.shape, dtype=bool)
		mask[tuple(coords.T)] = True
		markers, _ = ndi.label(mask)
		cell_labels = watershed(-distance, markers, mask=image).astype(float)
		return cell_labels

	def call_peaks(self, **kwargs):
		peaks = {}
		dogmats = {}
		for key, value in self.flurophore_dict.items():
			dogmat, (row, col) = get_peaks(value, **kwargs)	
			peaks[key] = (row, col)
			dogmats[key] = dogmat
		return peaks, dogmats
	
	def make_peak_cell_df(self):
		rowdata = []
		for key, (rows, cols) in self.peak_dict.items():
			cell_labels = self.cell_labels[rows, cols]
			rowdata.extend(zip([key]*len(rows), cell_labels, rows, cols))
		peak_cell_df = pd.DataFrame(rowdata, columns=['fluorophore', 'cell_label', 'row', 'col'])
		return peak_cell_df

	def set_label_counts(self):
		values, counts = np.unique(self.cell_labels, return_counts=True)
		values, counts = values[:-1], counts[:-1]
		return values, counts

	def permute_labels(self):
		unique_values = np.unique(self.cell_labels)
		permuted_values = np.random.permutation(unique_values)
		value_mapping = {original: permuted for original, permuted in zip(unique_values, permuted_values)}
		value_mapping[np.nan] = np.nan
		permuted_array = np.vectorize(value_mapping.get)(self.cell_labels)
		permuted_array[self.cell_labels==0] = np.nan
		return permuted_array




	def make_qc_plots(n = 5):
		pass
		# for i
	def plot_cell(self, cell_id):
		X, Y = np.where(self.cell_labels == cell_id)
		sl1, sl2 = slice(min(X), max(X)), slice(min(Y), max(Y))
		self.make_qc_plot(sl1, sl2)

	def make_qc_plot(self, sl1, sl2, fgsz=(2, 2)):
		s1, e1, s2, e2 = sl1.start, sl1.stop, sl2.start, sl2.stop 

		fig, axs = init_subplots(5, 1, space=.2, fgsz=(2, 2))
		axs[0].matshow(self.dapi_array[sl1, sl2], cmap='coolwarm', extent=[s2, e2, e1, s1])
		axs[1].matshow(self.all_original_cells[sl1, sl2], cmap='tab20', extent=[s2, e2, e1, s1])
		axs[2].matshow(self.small_cells[sl1, sl2], cmap='tab20', extent=[s2, e2, e1, s1])
		axs[3].matshow(self.large_cells[sl1, sl2], cmap='tab20', extent=[s2, e2, e1, s1])
		axs[4].matshow(self.cell_labels[sl1, sl2], cmap='tab20', extent=[s2, e2, e1, s1])

		axs[0].set_title("DAPI")
		axs[1].set_title("All Nuclei")
		axs[2].set_title("Too Small Nuclei")
		axs[3].set_title("Too Large Nuclei")
		axs[4].set_title("Valid Nuclei")
		for ax in axs:
			ax.axhline((s1+e1)//2, linestyle='--', color='lightgray')
			ax.axvline((s2+e2)//2, linestyle='--', color='lightgray')
			ax.set_xticks([])
			ax.set_yticks([])
			
		tmp2 = self.cell_labels[sl1, sl2].copy()
		tmp2[np.isnan(tmp2)] = 0
		num_features = len(np.unique(tmp2))
		us = np.unique(tmp2)
		for region_label in range(1, num_features):
			region = (tmp2 == us[region_label]).astype(int)
			contours = find_contours(region, 0.5)
			for contour in contours:
				for ax in axs:
					ax.plot(contour[:, 1]+s2, contour[:, 0]+s1, linewidth=1, color='black')
		
		n = len(self.flurophore_dict)
		fig, axs = init_subplots(2*n, 2, space=.4, fgsz=fgsz)
		for c, (name, array) in enumerate(self.flurophore_dict.items()):
			axs[c].matshow(array[sl1, sl2], cmap='coolwarm', extent=[s2, e2, e1, s1])
			axs[c].set_title(f'{name} (raw)')
		
			peaks_y, peaks_x = self.peak_dict[name]	
			good_peaks = (peaks_y < e1) & (peaks_y > s1) & (peaks_x < e2) & (peaks_x > s2)
			peaks_y, peaks_x = peaks_y[good_peaks], peaks_x[good_peaks]
			dogmat = self.dogmat_dict[name]
			axs[c+n].matshow(dogmat[sl1, sl2], cmap='coolwarm', extent=[s2, e2, e1, s1])
			axs[c+n].scatter(peaks_x, peaks_y-.5, edgecolor='black', facecolor='none', marker='s', s=200, zorder=3)		
			axs[c+n].set_title(f'{name} (DOG)')
		for ax in axs:
			ax.axhline((s1+e1)//2, linestyle='--', color='lightgray')
			ax.axvline((s2+e2)//2, linestyle='--', color='lightgray')
			ax.set_xticks([s2, e2])
			ax.set_yticks([s1, e1])
		for region_label in range(1, num_features):
			region = (tmp2 == us[region_label]).astype(int)
			contours = find_contours(region, 0.5)
			for contour in contours:
				for ax in axs:
					ax.plot(contour[:, 1]+s2, contour[:, 0]+s1, linewidth=1, color='black')


def get_peaks(mat,threshold_rel=.5):
	dogmat = dog(mat, 1, 2)
	maxes = peak_local_max(dogmat, min_distance=4, threshold_rel=threshold_rel)
	if len(maxes) > 0:
		row, col = list(zip(*maxes))
		row, col = list(map(np.array, [row, col]))
	else:
		row, col = np.array([]), np.array([])
	return dogmat, (row, col)

def get_distances(self):
	pass




class ControlFishExperiment(FishExperiment):
	def __init__(self, path_to_image, lowco = 3200, highco = 30_000, threshold_rel = 0.5):
		self.threshold_rel = threshold_rel
		self.lowco = lowco
		self.highco = highco
		self.img = tiff.imread(path_to_image)
		self.dapi_array = np.array(self.img)[2]
		self.dapi_array = self.dapi_array - np.median(self.dapi_array)
		self.dapi_array[self.dapi_array < 0] = 0

		self.flurophore_dict = {}

		self.flurophore_dict["Arl4c"] = np.array(self.img)[1]
		self.flurophore_dict["Ctla4"] = np.array(self.img)[0]

		self.cell_labels = self.segment_cells(co=1000)
		self.cell_labels = self.permute_labels()
		self.all_original_cells = self.cell_labels.copy()
		self.cell_ids, self.cell_counts = self.set_label_counts()
		self.cell_labels = self.remove_outlier_cells()
		self.peak_dict, self.dogmat_dict = self.call_peaks(threshold_rel=threshold_rel)
		self.peak_cell_df = self.make_peak_cell_df()






def get_distance_df(fish_expt):
    fluor_neighbors = []
    n_cells = len(fish_expt.peak_cell_df['cell_label'].unique())
    for cell_label, i in fish_expt.peak_cell_df.groupby('cell_label'):
        us = i['fluorophore'].unique()
        if (fish_expt.valid_cell_ids == cell_label).sum() == 0:
            continue
        if cell_label == 0:
            continue
        if len(us) <= 1:
            continue
        for u1 in us:
            for u2 in us:
                if u1 == u2:
                    continue
                idx1 = i['fluorophore'] == u1
                idx2 = i['fluorophore'] == u2
                if (idx1.sum() != 2) or (idx2.sum() != 2):
                    continue
                sub1 = i.loc[idx1]
                sub2 = i.loc[idx2]
                for _ in range(2):
                    if (len(sub1) == 0) or (len(sub2) == 0):
                        continue
                
                    d_row = np.add.outer(sub1['row'].values, -sub2['row'].values)
                    d_col = np.add.outer(sub1['col'].values, -sub2['col'].values)
                
                    dists = np.sqrt(d_row * d_row + d_col * d_col)
                    ind1, ind2 = np.unravel_index(np.argmin(dists), dists.shape)
                    nucleus_size = fish_expt.cell_counts[fish_expt.cell_ids == cell_label][0]
                    
                    row_data = [u1, u2, cell_label, nucleus_size, dists[ind1, ind2]]
                    fluor_neighbors.append(row_data)
                    sub1 = sub1.drop(sub1.index[ind1])
                    sub2 = sub2.drop(sub2.index[ind2])
    df = pd.DataFrame(fluor_neighbors, columns = ['from', 'to', 'cell', 'nucleus_size', 'distance']
            )
    n_cells_with_fluors = len(np.unique(df['cell']))
    print('# Cells:', n_cells, '; % with fluor: ', n_cells_with_fluors/n_cells)
    return df

