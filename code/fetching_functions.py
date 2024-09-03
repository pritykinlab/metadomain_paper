import cooler
from aux_functions import remove_chr, add_chr_to_anc

def fetch_from_cooler(cool, l1, l2):
	if 'chr' in cool.chromnames[0]:
		l1, l2 = map(add_chr_to_anc, [l1, l2])
	else:
		l1, l2 = map(remove_chr, [l1, l2])
	m = cool.matrix().fetch(l1, l2)
	return m