from pyjaspar import jaspardb
import pandas as pd

def get_jaspar_motif_metadata(motif_ids):
    # Create a JASPAR database object
    jdb = jaspardb()
    
    motif_metadata = {}
    
    for motif_id in motif_ids:
        try:
            # Fetch motif by ID
            motif = jdb.fetch_motif_by_id(motif_id)
            
            # Collect various metadata
            metadata = {
                'name': motif.name if hasattr(motif, 'name') else 'Unknown',
                'family': motif.tf_family if hasattr(motif, 'tf_family') else 'Unknown',
                'class': motif.tf_class if hasattr(motif, 'tf_class') else 'Unknown',
                'species': motif.species if hasattr(motif, 'species') else 'Unknown',
                'matrix': motif.matrix if hasattr(motif, 'matrix') else 'Unknown',
                'sequence_logo': motif.sequence_logo() if hasattr(motif, 'sequence_logo') else 'Unknown'
            }
            motif_metadata[motif_id] = metadata
        except Exception as e:
            motif_metadata[motif_id] = f'Error: {e}'
    
    return motif_metadata

def parse_name(parsed_name, row):
    if parsed_name == 'Regulators of differentiation':
        if "MEF" in row.name:
            parsed_name = 'MEF'
        else:
            raise Exception
    elif parsed_name == 'More than 3 adjacent zinc fingers':
        parsed_name = row.name
    elif parsed_name == 'Other  with up to three adjacent zinc fingers':
        parsed_name = row.name
    elif parsed_name == 'Interferon-regulatory':
        parsed_name = 'IRF'
    elif parsed_name == 'Heteromeric CCAAT-binding':
        parsed_name = 'NFYB/C'
    elif parsed_name == 'Nuclear factor 1':
        parsed_name = 'NFI'
    elif parsed_name == 'Factors with multiple dispersed zinc fingers':
        parsed_name = row.name
    elif parsed_name == 'Three-zinc finger Kruppel':
        parsed_name = row.name
    elif parsed_name == 'Early B-Cell Factor':
        if "ebf" in row.name.lower():
            parsed_name = 'EBF'
        else:
            raise Exception
    elif parsed_name == 'FTZF1related(NR5A)':
        if "nr5" in row.name.lower():
            parsed_name = 'NR5A1/2'
        else:
            raise Exception
    elif parsed_name == 'GCNF receptors (NR6)':
        if "nr6" in row.name.lower():
            parsed_name = 'NR6A1'
        else:
            raise Exception
    elif parsed_name == 'FTZF1related(NR5A)':
        if "nr5" in row.name.lower():
            parsed_name = 'NR5A1/2'
        else:
            raise Exception
    elif parsed_name == 'RXR receptors (NR2)':
        if "nr2" in row.name.lower():
            parsed_name = 'NR2'
        elif 'hnf' in row.name.lower():
            parsed_name = 'HNF4G'
        elif "rxr" in row.name.lower():
            parsed_name = 'RXR'
        else:
            raise Exception   
    elif parsed_name == 'Steroid hormone receptors (NR3)':
        if "nr3" in row.name.lower():
            parsed_name = 'NR3'
        elif "esr" in row.name.lower():
            parsed_name = 'ESSR'
        elif "ar" in row.name.lower():
            parsed_name = 'Ar'
        else:
            raise Exception
    elif parsed_name == 'Thyroid hormone receptor  (NR1)':
        if "nr1" in row.name.lower():
            parsed_name = 'NR1'
        elif "rar" in row.name.lower():
            parsed_name = 'RARA'
        elif "ppar" in row.name.lower():
            parsed_name = 'PPAR'
        elif "thrb" in row.name.lower():
            parsed_name = 'THRB'
        elif "vdr" in row.name.lower():
            parsed_name = 'VDR'
        elif "ror" in row.name.lower():
            parsed_name = 'ROR'
        elif "thra" in row.name.lower():
            parsed_name = 'THRA'
        else:
            raise Exception

    return parsed_name


def parse_coarse_name(motif_metadata):
    motif_metadata['coarse_parsed_family'] = motif_metadata['parsed_family'].copy()
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("KLF"), 'coarse_parsed_family'] = 'KLF'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("ZNF"), 'coarse_parsed_family'] = 'ZNF'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("ZBTB"), 'coarse_parsed_family'] = 'ZBTB'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("EGR"), 'coarse_parsed_family'] = 'EGR'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("PLAG"), 'coarse_parsed_family'] = 'PLAG'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("SP"), 'coarse_parsed_family'] = 'SP'
    motif_metadata.loc[motif_metadata['parsed_family'].str.lower().str.contains("gli"), 'coarse_parsed_family'] = 'GLI'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("CTCF"), 'coarse_parsed_family'] = 'CTCF'
    motif_metadata.loc[motif_metadata['parsed_family'].str.contains("OSR"), 'coarse_parsed_family'] = 'OSR'
    motif_metadata.loc[motif_metadata['parsed_family'].str.lower().str.contains("zfp"), 'coarse_parsed_family'] = 'ZFP'
    return motif_metadata

def get_motif_metadata(motif_id_to_name_dict):
    ''' JASPAR ID to deduplicated name dictionary '''
    motif_metadata = get_jaspar_motif_metadata(motif_id_to_name_dict.keys())
    motif_metadata = pd.DataFrame(motif_metadata).T
    motif_metadata['dedup_name'] = [motif_id_to_name_dict[x] for x in motif_metadata.index]
    motif_metadata = motif_metadata.set_index('dedup_name')
    motif_metadata['n_family'] = motif_metadata['family'].apply(len).sort_values()
    
    parsed_family = []
    for _, row in motif_metadata.iterrows():
        i = row.family
        if i == 'Unknown':
            parsed_name = 'Unknown'
        elif len(i) == 1:
            parsed_name = i[0].replace("-related", "").replace("factors", "").strip()
            parsed_name = parse_name(parsed_name, row)
        else:
            parsed_names = [_.replace("-related", "").replace("factors", "").strip() for _ in i]
            parsed_names = [parse_name(_, row) for _ in parsed_names]
            parsed_name = '-'.join(sorted(parsed_names))
        parsed_family.append(parsed_name)
    
    motif_metadata['parsed_family'] = parsed_family
    motif_metadata = parse_coarse_name(motif_metadata)
    return motif_metadata

