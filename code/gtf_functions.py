import pandas as pd
def make_tss_df(transcripts_file, PARSED_CHROMS = {'chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 
                                                   'chr18', 'chr19', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 
                                                   'chrX'}):
    transcript2gene = {}
    transcript2info = {}
    transcript_rows = []
    tss_rows = []

    for i in transcripts_file :
        if 'transcript_id' not in i.attrs:
            continue
        if 'misc_RNA' in i.attrs['gene_type']:
            continue
        if 'rRNA' in i.attrs['gene_type']:
            continue
        if 'transcript_support_level' not in i.attrs:
            continue

        transcript = i.attrs['transcript_id']
        gene = i.attrs['gene_id']
        gene_name = i.attrs['gene_name']
        transcript2gene[transcript] = gene
        chrom = i.chrom
        start = i.start
        end = i.end
        strand = i.strand
        transcript_type = i.attrs['transcript_type']
        transcript_name = i.attrs['transcript_name']
        transcript_id = i.attrs['transcript_id']
        support_level = i.attrs['transcript_support_level']
        if strand == '-':
            row = [gene_name, chrom, end-1, end, strand, support_level, transcript_type, transcript_name, transcript_id]
            tss_rows.append(row)
        elif strand == '+':
            row = [gene_name, chrom, start, start+1, strand, support_level, transcript_type, transcript_name, transcript_id]
            tss_rows.append(row)
        transcript_rows.append([gene_name, chrom, start, end, strand, support_level, transcript_type, transcript_name, transcript_id])
        
    my_tss_df = pd.DataFrame(tss_rows)
    my_tss_df = my_tss_df[my_tss_df[1].isin(PARSED_CHROMS.union({"chrY"}))]
    my_tss_df.columns = ['gene_name', 'chrom', 'start', 'end', 'strand', 'support_level', 'transcript_type', 'transcript_name', 'transcript_id']
    my_tss_df = my_tss_df[['chrom', 'start', 'end', 'gene_name', 'strand', 'support_level', 'transcript_type', 'transcript_name', 'transcript_id']]
    
    my_transcript_df = pd.DataFrame(transcript_rows)
    my_transcript_df = my_transcript_df[my_transcript_df[1].isin(PARSED_CHROMS.union({"chrY"}))]
    my_transcript_df.columns = ['gene_name', 'chrom', 'start', 'end', 'strand', 'support_level', 'transcript_type', 'transcript_name', 'transcript_id']
    my_transcript_df = my_transcript_df[['chrom', 'start', 'end', 'gene_name', 'strand', 'support_level', 'transcript_type', 'transcript_name', 'transcript_id']]
    
    return my_tss_df, my_transcript_df