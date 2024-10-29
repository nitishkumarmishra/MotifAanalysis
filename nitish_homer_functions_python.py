import os
import pandas as pd
import gzip
import re
import pandas as pd
from scipy.stats import hypergeom




def get_transcriptome_feature_fasta_filename(species, feat):
    items = feat
    for i in range(len(items), 0, -1):
        file = "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/m39/Mus_musculus.GRCm39.104.rdna_rn18s.transcriptome.gtf.gz"
        basef = file.replace(".gtf.gz", "")
        featf = f"{basef}.{ '_'.join(items[:i]) }.fa.gz"
        if os.path.exists(featf):
            break
    return featf




def read_gtf(file, header=False):
    with gzip.open(file, 'rt') as f:
        in_gtf = pd.read_csv(f, sep='\t', header=None, comment='#', names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    return in_gtf



def load_transcriptome_gtf(species, verbose=False):
    if species == "m39":
        fname = "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/m39/Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz"
    else:
        raise ValueError(f"ERROR! unrecognized species=[{species}]!!!")

    if verbose:
        print(f"\t\t\tloading transcriptome gtf: [{fname}]")

    indf = read_gtf(file=fname)
    return indf



def get_gtf_attribute_field(gtf, field, num=False):
    pattern = re.compile(f'{field} "([^"]+)"')
    extracted = gtf['attribute'].apply(lambda x: pattern.search(x).group(1) if pattern.search(x) else None)
    if num:
        extracted = pd.to_numeric(extracted)
    return extracted



def prepare_transcriptome_gtf(gtf):
    if 'transcript_id' not in gtf.columns:
        gtf['transcript_id'] = get_gtf_attribute_field(gtf, 'transcript_id')
    if 'gene_name' not in gtf.columns:
        gtf['gene_name'] = get_gtf_attribute_field(gtf, 'gene_name')
    if 'gene_id' not in gtf.columns:
        gtf['gene_id'] = get_gtf_attribute_field(gtf, 'gene_id')
    return gtf







def load_homer_motif(motif):
    # Read the header
    with open(motif, 'r') as file:
        header = file.readline().strip()

    # Read the matrix
    inmat = pd.read_csv(motif, sep='\t', header=None, skiprows=1, quoting=3).values

    inmotif = {'header': header, 'matrix': inmat}

    return inmotif





def load_homer_motif_annotation(hits, bed=False):
    try:
        indf = pd.read_csv(hits, sep="\t", header=None, quoting=3, names=["chr", "offset", "sequence", "motif", "strand", "score"])
    except Exception as e:
        return None

    if bed:
        indf = homer_motif_annotation_to_bed(indf)

    return indf

def homer_motif_annotation_to_bed(annot):
    # Implement the conversion logic here
    pass






def test_hypergeometric_enrichment(universe, classes, groups, verbose=False):
    def verb(message):
        if verbose:
            print(message)

    ##### check
    verb("\tchecking.\n")

    assert isinstance(groups, pd.DataFrame)
    assert isinstance(classes, pd.DataFrame)
    assert groups.notnull().all().all()
    assert classes.notnull().all().all()
    assert all(col in classes.columns for col in ["element", "group"])
    assert all(col in groups.columns for col in ["element", "group"])
    assert set(groups['element']).issubset(set(universe))
    assert set(classes['element']).issubset(set(universe))

    ##### prepare
    verb("\tpreparing.\n")

    # unique
    classes = classes[['group', 'element']].drop_duplicates()
    groups = groups[['group', 'element']].drop_duplicates()

    # prepare
    pop_size = len(set(universe))

    classes['group'] = pd.Categorical(classes['group'], categories=classes['group'].unique())

    ctab = classes['group'].value_counts()
    classdf = pd.DataFrame({'group': ctab.index, 'pop_success': ctab.values})

    ##### check again
    verb("\tchecking again.\n")

    assert not groups.isnull().any().any()
    assert not classes.isnull().any().any()
    assert not classdf.isnull().any().any()

    resdf = pd.DataFrame()
    for groupx in groups['group'].unique():
        verb(f"\t\t{groupx}\n")

        genes = groups.loc[groups['group'] == groupx, 'element']

        subclass = classes[classes['element'].isin(genes)]
        stab = subclass['group'].value_counts()

        scdf = pd.DataFrame({'group': stab.index, 'sample_success': stab.values})
        overdf = pd.merge(classdf, scdf, on='group', how='outer').fillna(0)
        assert not overdf.isnull().any().any()

        q = overdf['sample_success']
        m = overdf['pop_success']
        n = pop_size - overdf['pop_success']
        k = len(genes)

        # Recall that:
        #   lower.tail = TRUE  -->  P[X <= x]
        #   lower.tail = FALSE -->  P[X  > x]

        utail = hypergeom.sf(q - 1, pop_size, m, k)
        dens_x = hypergeom.pmf(q, pop_size, m, k)
        p_over = utail + dens_x
        p_under = 1.0 - utail

        rdf = pd.DataFrame({
            'class': overdf['group'],
            'group': groupx,
            'success_sample': q,
            'sample_size': k,
            'success_population': m,
            'population_size': pop_size,
            'frac_success_sample': q / k,
            'frac_success_population': m / pop_size,
            'p_overrepresented': p_over,
            'p_underrepresented': p_under
        })
        resdf = pd.concat([resdf, rdf], ignore_index=True)

    resdf = resdf.sort_values(by=['group', 'p_overrepresented', 'success_population', 'sample_size'])

    return resdf





def de_analysis_type_to_directory_and_tag(type, project):
    if type == "mrna":
        bdir = "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\\ codes/out/limma-voom.mrna"
        btag = "limma-voom.mrna"
    elif type == "rrna-var":
        bdir = "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\\ codes/out/limma-voom.mrna"
        btag = "limma-voom.rrna-var"
    else:
        raise ValueError(f"ERROR! unrecognized type=[{type}]!!!")

    res = {"dir": bdir, "tag": btag}
    return res




def load_project_DE_test_tables(file_names, project, individual):
    df_de = pd.DataFrame()
    for filename in file_names:
        df1 = pd.read_csv(filename, sep="\t", header=0, index_col=0, comment="#", quotechar='"')
        cols = df1.columns
        condition1 = cols
        condition2 = cols
        condition1 = "tgfb48" if condition1 == "tgfb" else "unt48" if condition1 == "unt" else condition1
        condition2 = "tgfb48" if condition2 == "tgfb" else "tgfbCX5461" if condition2 in ["tgfbCX", "tgfbCX5461100nm"] else condition2
        df1.columns = ['expression1', 'expression2', 'log2FC', 'FDR', 'p.value']
        df1['project'] = project
        df1['individual'] = individual
        df0 = pd.DataFrame({'gene_name': df1.index, 'condition1': condition1, 'condition2': condition2})
        df1.reset_index(drop=True, inplace=True)
        df1 = pd.concat([df0, df1.reset_index(drop=True)], axis=1)
        df_de = pd.concat([df_de, df1], ignore_index=True)
    return df_de





def build_gene_set(de, project, conditions, reversible=None, include_notSig=None, FDR_thresh=0.05, element="gene_name", gene_set=None, gname='default', verbose=False):
    if reversible is None:
        reversible = len(conditions) > 2
    if include_notSig is None:
        include_notSig = not reversible

    if gene_set is None:
        gene_set = {"sets": {}, "factors": {}, "tables": {}}

    if isinstance(project, pd.DataFrame):
        project = project['project'].unique()

    gset = {}

    if not reversible:
        if verbose:
            print("\tpairwise.")

        assert len(conditions) == 2

        ### DE
        gset['up'] = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] <= FDR_thresh) & (de['log2FC'] > 0) & (de['project'].isin(project))]
        gset['down'] = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] <= FDR_thresh) & (de['log2FC'] < 0) & (de['project'].isin(project))]

        ### factor and notsig
        if include_notSig:
            gset['notSig'] = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] > FDR_thresh) & (de['project'].isin(project))]
            gfac = pd.Categorical(list(gset.keys()), categories=["notSig", "down", "up"])
        else:
            gfac = pd.Categorical(list(gset.keys()), categories=["down", "up"])

        ### table
        gtab = de[(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['project'].isin(project))]

    else:
        if verbose:
            print("\treversible.")

        assert len(conditions) == 3
        assert not include_notSig

        ### rev
        up1 = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] <= FDR_thresh) & (de['log2FC'] > 0) & (de['project'].isin(project))]
        down1 = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] <= FDR_thresh) & (de['log2FC'] < 0) & (de['project'].isin(project))]

        up2 = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] <= FDR_thresh) & (de['log2FC'] > 0) & (de['project'].isin(project))]
        down2 = de[element][(de['condition1'] == conditions) & (de['condition2'] == conditions) & (de['FDR'] <= FDR_thresh) & (de['log2FC'] < 0) & (de['project'].isin(project))]

        gset['upDown'] = set(up1).intersection(down2)
        gset['downUp'] = set(down1).intersection(up2)

        ### factor
        gfac = pd.Categorical(list(gset.keys()), categories=["downUp", "upDown"])

        ### table
        gtab = de[((de['condition1'] == conditions) & (de['condition2'] == conditions)) | ((de['condition1'] == conditions) & (de['condition2'] == conditions)) & (de['project'].isin(project))]

    gene_set['sets'][gname] = gset
    gene_set['factors'][gname] = gfac
    gene_set['tables'][gname] = gtab
    return gene_set








