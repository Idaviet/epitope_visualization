import pandas as pd
import numpy as np
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
import os


def map_epitopes_to_abdab(abdab_file, epitope_file):
    """Map IEDB epitope entries to CoV-AbDab antibody entries.

    Matches antibodies from CoV-AbDab to epitopes from IEDB by comparing
    CDRH3, CDRL3, VH, VL, V-genes, and J-genes. For each antibody in the
    CoV-AbDab database, finds matching epitopes in the IEDB database and
    populates the CoV-AbDab DataFrame with matched epitope information.
    Also computes consensus epitopes when multiple chain-level matches agree.

    Args:
        abdab_file (str): Path to the CoV-AbDab CSV or Excel file.
        epitope_file (str): Path to the IEDB epitope CSV or Excel file.

    Returns:
        pd.DataFrame: CoV-AbDab DataFrame augmented with epitope mapping
            columns including per-chain epitope matches, gene-level matches,
            consensus epitopes, and match counts.
    """

    
    abdab_df = pd.read_excel(abdab_file) if abdab_file.endswith('.xlsx') else pd.read_csv(abdab_file)
    
    epitope_df = pd.read_csv(epitope_file) if epitope_file.endswith('.csv') else pd.read_excel(epitope_file)
    epitope_df.rename(columns={"Epitope - Name": "epitope"}, inplace = True)


    epitope_df['cdrh3'] = epitope_df['Chain 1 - CDR3 Curated'].combine_first(epitope_df['Chain 1 - CDR3 Calculated'])
    epitope_df['cdrl3'] = epitope_df['Chain 2 - CDR3 Curated'].combine_first(epitope_df['Chain 2 - CDR3 Calculated'])

    abdab_df['cdrh3_epitopes'] = np.nan
    abdab_df['n_cdrh3_epi'] = np.nan

    abdab_df['cdrl3_epitopes'] = np.nan
    abdab_df['n_cdrl3_epi'] = np.nan

    abdab_df['cdr3_hl_epitopes'] = np.nan
    abdab_df['n_cdr3_hl_epi'] = np.nan

    abdab_df['vh_epitopes'] = np.nan
    abdab_df['n_vh_epi'] = np.nan

    abdab_df['vl_epitopes'] = np.nan
    abdab_df['n_vl_epi'] = np.nan


    ### to add
    abdab_df['Hv_gene_epitopes'] = np.nan
    abdab_df['n_Hv_gene_epi'] = np.nan

    abdab_df['Lv_gene_epitopes'] = np.nan
    abdab_df['n_Lv_gene_epi'] = np.nan

    abdab_df['Hj_gene_epitopes'] = np.nan
    abdab_df['n_Hj_gene_epi'] = np.nan

    abdab_df['Lj_gene_epitopes'] = np.nan
    abdab_df['n_Lj_gene_epi'] = np.nan

    abdab_df['cdr3_cons_epitopes'] = np.nan
    abdab_df['n_cdr3_epi'] = np.nan

    abdab_df['vhl_cons_epitopes'] = np.nan
    abdab_df['n_vhl_epi'] = np.nan

    
    # ### to add
    # abdab_df['Hvj_cons_epitopes'] = np.nan
    # abdab_df['n_Hvj_epi'] = np.nan

    # abdab_df['Lvj_cons_epitopes'] = np.nan
    # abdab_df['n_Lvj_epi'] = np.nan

    # abdab_df['HLvj_cons_epitopes'] = np.nan
    # abdab_df['n_HLvj_epi'] = np.nan



    abdab_df['all_epitopes'] = np.nan
    abdab_df['n_all_epis'] = np.nan

    abdab_df['consensus_epitope'] = np.nan

    for idx, row in epitope_df.iterrows():
        epitope = row['epitope']
        epitope_misc = f"{row['Receptor - IEDB Receptor ID']} | {row['Receptor - Reference Name']} | {row['Receptor - Type']}"

        cdrh3 = row['cdrh3']
        cdrl3 = row['cdrl3']
        vh = row['Chain 1 - Protein Sequence']
        vl = row['Chain 2 - Protein Sequence']

        hv_gene = np.nan
        lv_gene = np.nan

        hj_gene = np.nan
        lj_gene = np.nan

        
        if row['Chain 1 - Curated V Gene'] and isinstance(row['Chain 1 - Curated V Gene'], str) and not pd.isna(row['Chain 1 - Curated V Gene']):
            hv_gene = row['Chain 1 - Curated V Gene']
        elif row['Chain 1 - Calculated V Gene'] and isinstance(row['Chain 1 - Calculated V Gene'], str) and not pd.isna(row['Chain 1 - Calculated V Gene']):
            hv_gene = row['Chain 1 - Calculated V Gene']
            hv_gene = hv_gene.split('*')[0] if '*' in hv_gene else hv_gene

        if row['Chain 1 - Curated J Gene'] and isinstance(row['Chain 1 - Curated J Gene'], str) and not pd.isna(row['Chain 1 - Curated J Gene']):
            hj_gene = row['Chain 1 - Curated J Gene']
        elif row['Chain 1 - Calculated J Gene'] and isinstance(row['Chain 1 - Calculated J Gene'], str) and not pd.isna(row['Chain 1 - Calculated J Gene']):
            hj_gene = row['Chain 1 - Calculated J Gene']
            hj_gene = hj_gene.split('*')[0] if '*' in hj_gene else hj_gene


        
        if row['Chain 2 - Curated V Gene'] and isinstance(row['Chain 2 - Curated V Gene'], str) and not pd.isna(row['Chain 2 - Curated V Gene']):
            lv_gene = row['Chain 2 - Curated V Gene']
        elif row['Chain 2 - Calculated V Gene'] and isinstance(row['Chain 2 - Calculated V Gene'], str) and not pd.isna(row['Chain 2 - Calculated V Gene']):
            lv_gene = row['Chain 2 - Calculated V Gene']
            lv_gene = lv_gene.split('*')[0] if '*' in lv_gene else lv_gene

        if row['Chain 2 - Curated J Gene'] and isinstance(row['Chain 2 - Curated J Gene'], str) and not pd.isna(row['Chain 2 - Curated J Gene']):
            lj_gene = row['Chain 2 - Curated J Gene']
        elif row['Chain 2 - Calculated J Gene'] and isinstance(row['Chain 2 - Calculated J Gene'], str) and not pd.isna(row['Chain 2 - Calculated J Gene']):
            lj_gene = row['Chain 2 - Calculated J Gene']
            lj_gene = lj_gene.split('*')[0] if '*' in lj_gene else lj_gene

        iter_list = [[cdrh3, 'CDRH3', 'cdrh3_epitopes'], [cdrl3, 'CDRL3', 'cdrl3_epitopes'], [vh, 'VHorVHH', 'vh_epitopes'], [vl, 'VL', 'vl_epitopes'], [hv_gene, 'Heavy V Gene', 'Hv_gene_epitopes'], [lv_gene, 'Light V Gene', 'Lv_gene_epitopes'], [hj_gene, 'Heavy J Gene', 'Hj_gene_epitopes'], [lj_gene, 'Light J Gene', 'Lj_gene_epitopes']]

        for i in iter_list:
            epi_seq, ab_col, epi_col = i
            # print(epi_seq, ab_col, epi_col)

            if isinstance(epi_seq, str) and epi_seq.strip() != '':
                mask = abdab_df[ab_col].str.contains(epi_seq, na=False)
                # n_epi = 0
                for idx in abdab_df[mask].index:
                    current = abdab_df.at[idx, epi_col]
                    if pd.isna(current) or current == '':
                        abdab_df.at[idx, epi_col] = epitope
                        # n_epi += 1 
                        # abdab_df.at[idx, 'n_cdrh3_epi'] = n_epi

                    else:
                        if epitope not in str(current).split(' | '):
                            abdab_df.at[idx, epi_col] = str(current) + ' | ' + epitope
                            # n_epi += 1
                            # abdab_df.at[idx, 'n_cdrh3_epi'] = n_ep
                            
    for idx, row in abdab_df.iterrows():
        cdrh3_epi = row['cdrh3_epitopes']
        n_cdrh3_epi = int(0) if pd.isna(cdrh3_epi) else len(str(cdrh3_epi).split(' | '))
        abdab_df.at[idx, 'n_cdrh3_epi'] = n_cdrh3_epi

        cdrl3_epi = row['cdrl3_epitopes']
        n_cdrl3_epi = int(0) if pd.isna(cdrl3_epi) else len(str(cdrl3_epi).split(' | '))
        abdab_df.at[idx, 'n_cdrl3_epi'] = n_cdrl3_epi

        vh_epi = row['vh_epitopes']
        n_vh_epi = int(0) if pd.isna(vh_epi) else len(str(vh_epi).split(' | '))
        abdab_df.at[idx, 'n_vh_epi'] = n_vh_epi
        
        vl_epi = row['vl_epitopes']
        n_vl_epi = int(0) if pd.isna(vl_epi) else len(str(vl_epi).split(' | '))
        abdab_df.at[idx, 'n_vl_epi'] = n_vl_epi

        hv_gene_epi = row['Hv_gene_epitopes']
        n_hv_gene_epi = int(0) if pd.isna(hv_gene_epi) else len(str(hv_gene_epi).split(' | '))
        abdab_df.at[idx, 'n_Hv_gene_epi'] = n_hv_gene_epi

        lv_gene_epi = row['Lv_gene_epitopes']
        n_lv_gene_epi = int(0) if pd.isna(lv_gene_epi) else len(str(lv_gene_epi).split(' | '))
        abdab_df.at[idx, 'n_Lv_gene_epi'] = n_lv_gene_epi

        hj_gene_epi = row['Hj_gene_epitopes']
        n_hj_gene_epi = int(0) if pd.isna(hj_gene_epi) else len(str(hj_gene_epi).split(' | '))
        abdab_df.at[idx, 'n_Hj_gene_epi'] = n_hj_gene_epi

        lj_gene_epi = row['Lj_gene_epitopes']
        n_lj_gene_epi = int(0) if pd.isna(lj_gene_epi) else len(str(lj_gene_epi).split(' | '))
        abdab_df.at[idx, 'n_Lj_gene_epi'] = n_lj_gene_epi


        n_cdr3_epi = n_cdrh3_epi + n_cdrl3_epi
        abdab_df.at[idx, 'n_cdr3_hl_epi'] = n_cdr3_epi
        n_vhl_epi = n_vh_epi + n_vl_epi
        abdab_df.at[idx, 'n_vhl_epi'] = n_vhl_epi
        
        n_hvj_epi = n_hv_gene_epi + n_hj_gene_epi
        abdab_df.at[idx, 'n_Hvj_epi'] = n_hvj_epi
        n_lvj_epi = n_lv_gene_epi + n_lj_gene_epi
        abdab_df.at[idx, 'n_Lvj_epi'] = n_lvj_epi
        n_hlvj_epi = n_hvj_epi + n_lvj_epi
        abdab_df.at[idx, 'n_HLvj_epi'] = n_hlvj_epi


        n_total_epi = n_cdr3_epi + n_vhl_epi
        abdab_df.at[idx, 'n_all_epis'] = n_total_epi


        if cdrh3_epi == cdrl3_epi and cdrh3_epi == vh_epi and cdrh3_epi == vl_epi and pd.notna(cdrh3_epi):
            abdab_df.at[idx, 'consensus_epitope'] = cdrh3_epi

        if cdrh3_epi == cdrl3_epi and pd.notna(cdrh3_epi):
            abdab_df.at[idx, 'cdr3_cons_epitopes'] = cdrh3_epi

        if vh_epi == vl_epi and pd.notna(vh_epi):
            abdab_df.at[idx, 'vhl_cons_epitopes'] = vh_epi

        all_epitopes = f'{cdrh3_epi if pd.notna(cdrh3_epi) else ""} | {cdrl3_epi if pd.notna(cdrl3_epi) and cdrl3_epi not in [cdrh3] else ""} | {vh_epi if pd.notna(vh_epi) else ""} | {vl_epi if pd.notna(vl_epi) else ""}'

        all_epitopes = all_epitopes.strip(' | ')

        abdab_df.at[idx, 'all_epitopes'] = all_epitopes


    print(len(abdab_df['all_epitopes'].value_counts(dropna=True)))

    return abdab_df

### EXAMPLE USE ###
# abdab_file = "/Users/isaacdaviet/Downloads/CoV-AbDab_080224.csv"
# epitope_file = "/Users/isaacdaviet/Downloads/IEDB_ALL_huEpi_wBCRseqs.csv"

# abepi_df = map_epitopes_to_abdab(abdab_file, epitope_file)

# abepi_df.to_csv("/Users/isaacdaviet/Desktop/comp/CovAbDab_IEDB_AllHuEpi_epitopes.csv", index=False)

# column = 'all_epitopes'
# print(abepi_df[column].value_counts(dropna=True))
# print(len(abepi_df[column].value_counts(dropna=True)))



def get_pos_dict(epitope_file):
    """
    Builds a frequency dictionary of amino-acid/position tokens from an IEDB epitope CSV.

    Parses the ``'Epitope - Name'`` column of an IEDB export, handling semicolon-delimited
    alternate epitopes, colon-stripped prefixes, and comma-separated residue tokens. Each
    token is decomposed into an amino acid letter and a position number via regex; tokens
    that fail to yield both are silently skipped. Duplicate AApos strings across rows are
    aggregated by incrementing their count.

    Args:
        epitope_file (str): Path to the IEDB CSV file. Must contain an ``'Epitope - Name'``
            column.

    Returns:
        dict: Keys are ``AApos`` strings (e.g. ``'F486'``). Values are dicts with keys:
            ``'aa'`` (str) — amino acid letter(s), ``'pos'`` (int) — sequence position,
            ``'count'`` (int) — occurrences across all rows.
    """

    epitope_df = pd.read_csv(epitope_file)
    pos_dict = {}

    for idx, row in epitope_df.iterrows():
        epitope_str = row['Epitope - Name']

        if ';' in epitope_str:
            hold = []
            new_epitope_str = epitope_str.split(';')

            for i in new_epitope_str:
                if ':' in i:
                    i = i.split(':')[1]
                hold.append(i.strip())
            epitope_str = ', '.join(hold)

        epi_aa = epitope_str.split(',')
        epi_aa = [e.strip() for e in epi_aa]
        # print(row['Epitope ID - IEDB IRI'])
        
        for aa_pos in epi_aa:
            aa = re.findall(r'[A-Za-z]+', aa_pos)
            pos = re.findall(r'\d+', aa_pos)
            
            # print(pos, aa)
            # print(int(pos[0]))

            if aa_pos not in pos_dict.keys() and len(aa) > 0 and len(pos) > 0:
                # print(aa_pos, aa, pos)
                pos_dict[aa_pos] = {'aa': str(aa[0]), 'pos': int(pos[0]), 'count': int(1)}

            elif len(aa) > 0 and len(pos) > 0:
                pos_dict[aa_pos]['count'] += 1

    return pos_dict

def start_end_pos_dict_filter(pos_dict, start, end):
    """
    Filters a position dictionary to entries within a given residue range (inclusive).

    Iterates over the dictionary returned by :func:`get_pos_dict` and retains only
    entries whose ``'pos'`` value falls between ``start`` and ``end``.

    Args:
        pos_dict (dict): Dictionary from :func:`get_pos_dict`, keyed by ``AApos`` strings
            with values containing ``'pos'`` (int).
        start (int): Lower bound of the range (inclusive).
        end (int): Upper bound of the range (inclusive).

    Returns:
        dict: Subset of ``pos_dict`` containing only entries where ``start <= pos <= end``.
    """
    
    filtered_dict = {}
    for pos_name, vals in pos_dict.items():
        # print(vals)
        pos = vals['pos']
        if start and end:
            if pos >= start and pos <= end:
                filtered_dict[pos_name] = vals
        else:
            filtered_dict[pos_name] = vals

    return filtered_dict

def color_epitope_plot_by_property(plot_df, aa_col = 'aa', property_col= 'hydrophobicity', properties_table = '/Users/isaacdaviet/Desktop/comp/iedb/aa_physicochemical_properties.csv'):
    """Maps each row in a plot DataFrame to a physicochemical property value, averaging across multi-AA entries.

    Looks up each amino acid in a properties table, and for rows containing multiple AAs
    separated by ``'/'``, assigns the mean of their property values.

    Args:
        plot_df (pd.DataFrame): DataFrame with an amino acid column (see ``aa_col``).
        aa_col (str): Column name containing amino acid letters, optionally ``'/'``-delimited
            for multiple residues. Defaults to ``'aa'``.
        property_col (str): Column name in the properties table to use for the lookup.
            Defaults to ``'hydrophobicity'``.
        properties_table (str or pd.DataFrame): Path to a CSV file or an existing DataFrame
            with amino acid properties. Must contain an ``'AA'`` column and a column matching
            ``property_col``.

    Returns:
        pd.DataFrame: ``plot_df`` with a new column (named after ``property_col``) containing
        the averaged property value for each row.
    """
     
    prop_df = pd.read_csv(properties_table) if properties_table.endswith('.csv') else properties_table
    
    plot_df[property_col]= np.nan
    for idx, row in plot_df.iterrows():
        aa_list = row[aa_col].split('/')
        # print(aa_list)
        prop_vals = []

        for aa in aa_list:
            # print(aa)
            if aa in prop_df['AA'].values:
                # print(prop_df['AA'])
                val = prop_df.loc[prop_df['AA'] == aa, property_col].values[0]
                # print(val)

                prop_vals.append(val)
        
        if len(prop_vals) > 0:
            avg_val = sum(prop_vals) / len(prop_vals)
            plot_df.at[idx, property_col] = avg_val

    return plot_df

def plot_by_property_value(plot_df, property_col):
    """
    Renders a bar chart of epitope counts per position, colored by a physicochemical property.

    Bars are colored on a coolwarm gradient scaled to the min/max range of the property
    column.

    Args:
        plot_df (pd.DataFrame): Must contain ``'pos_label'`` (str, x-axis labels),
            ``'count'`` (int, bar heights), and a column matching ``property_col`` (float,
            used for bar coloring).
        property_col (str): Name of the column in ``plot_df`` whose values drive the
            color mapping.

    Returns:
        None: Displays the plot via :func:`matplotlib.pyplot.show`.
    """

    plt.figure(figsize=(50, 6))

    norm = plt.Normalize(plot_df[property_col].min(), plot_df[property_col].max())
    colors = cm.coolwarm(norm(plot_df[property_col]))
    
    plt.bar(plot_df['pos_label'], plot_df['count'], color=colors)
    plt.xticks(rotation=90)  # rotate x labels if there are many
    plt.xlabel('Position + Amino Acids')
    plt.ylabel('Count')
    plt.title(f'Counts per Position with Amino Acids _ colored by {property_col}')
    plt.tight_layout()
    plt.show()

def plot_all_physical_props(plot_df, aa_col = 'aa', properties_table = '/Users/isaacdaviet/Desktop/comp/iedb/aa_physicochemical_properties.csv'):
    """
    Generates a color-coded bar chart for every physicochemical property in the AA properties table.

    Iterates over all property columns (excluding ``'AA'``), calling
    :func:`color_epitope_plot_by_property` and :func:`plot_by_property_value` for each.

    Args:
        plot_df (pd.DataFrame): DataFrame containing amino acid and position data.
        aa_col (str): Column name with amino acid letters. Defaults to ``'aa'``.
        properties_table (str or pd.DataFrame): Path to CSV or DataFrame with amino acid
            properties. Must contain an ``'AA'`` column plus one or more property columns.

    Returns:
        pd.DataFrame: ``plot_df`` updated with property value columns from the final
        iteration.
    """

    prop_df = pd.read_csv(properties_table) if properties_table.endswith('.csv') else properties_table
    prop_cols = [col for col in prop_df.columns if col != 'AA']


    for prop in prop_cols:
        plot_df = color_epitope_plot_by_property(plot_df, aa_col = aa_col, property_col=prop, properties_table=properties_table)
        plot_by_property_value(plot_df, prop)
        plot_df = plot_df.rename(columns={'prop_val': prop})

    return plot_df 

def generate_epi_df(epitope_file, start_pos, end_pos):
    """
    Loads an epitope CSV, parses positions, filters by range, and returns a sorted DataFrame.

    Chains :func:`get_pos_dict`, :func:`start_end_pos_dict_filter`, and transposition into
    a DataFrame sorted by descending count.

    Args:
        epitope_file (str): Path to the IEDB CSV file.
        start_pos (int): Lower bound of the residue range (inclusive).
        end_pos (int): Upper bound of the residue range (inclusive).

    Returns:
        pd.DataFrame: Position-level DataFrame with columns ``'aa'``, ``'pos'``, ``'count'``,
        sorted by ``'count'`` descending.
    """

    pos_dict = get_pos_dict(epitope_file)
    filtered_dict = start_end_pos_dict_filter(pos_dict, start_pos, end_pos)
    epi_df = pd.DataFrame(filtered_dict).T
    return epi_df.sort_values(by='count', ascending=False)

def aggregate_plot_df(epi_df):
    """
    Aggregates an epitope DataFrame by position, summing counts and merging amino acids.

    Groups rows sharing the same ``'pos'`` value, sums their ``'count'`` columns, and
    joins their ``'aa'`` values into a ``'/'``-delimited string of sorted unique residues.
    Appends a ``'pos_label'`` column combining position and amino acids for use as x-axis
    labels in downstream plots.

    Args:
        epi_df (pd.DataFrame): DataFrame from :func:`generate_epi_df` with columns
            ``'aa'``, ``'pos'``, ``'count'``.

    Returns:
        pd.DataFrame: Aggregated DataFrame with columns ``'pos'``, ``'count'``, ``'aa'``,
        ``'pos_label'``.
    """
    
    plot_df = epi_df.reset_index()
    plot_df = plot_df.groupby('pos').agg({
        'count': 'sum',
        'aa': lambda x: '/'.join(sorted(set(x)))
    }).reset_index()
    plot_df['pos_label'] = plot_df['pos'].astype(str) + plot_df['aa']
    return plot_df

def plot_counts_bar(plot_df, title="Counts per Position with Amino Acids", figsize=(60, 6)):
    """
    Renders a simple bar chart of epitope counts per position.

    Args:
        plot_df (pd.DataFrame): DataFrame from :func:`aggregate_plot_df` with
            ``'pos_label'`` and ``'count'`` columns.
        title (str): Chart title. Defaults to ``'Counts per Position with Amino Acids'``.
        figsize (tuple): Figure dimensions as ``(width, height)`` in inches.
            Defaults to ``(60, 6)``.

    Returns:
        None: Displays the plot via :func:`matplotlib.pyplot.show`.
    """
    
    plt.figure(figsize=figsize)
    plt.bar(plot_df['pos_label'], plot_df['count'], color='skyblue')
    plt.xticks(rotation=90)
    plt.xlabel('Position + Amino Acids')
    plt.ylabel('Count')
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_max_normalized_property(res_props_df, prop_cols, figsize=(60, 6)):
    """
    Renders a bar chart where each position is colored by its dominant normalized property.

    For each property in ``prop_cols``, computes a min-max normalized value. Each bar is
    assigned the property with the highest normalized score at that position, and colored
    using a tab20 categorical palette with a corresponding legend.

    Args:
        res_props_df (pd.DataFrame): DataFrame containing ``'pos_label'``, ``'count'``,
            and a numeric column for each property in ``prop_cols``.
        prop_cols (list of str): Property column names to normalize and compare.
        figsize (tuple): Figure dimensions as ``(width, height)`` in inches.
            Defaults to ``(60, 6)``.

    Returns:
        pd.DataFrame: ``res_props_df`` augmented with ``'{prop}_norm'`` columns,
        a ``'max_norm'`` column, and a ``'max_property'`` column indicating the
        dominant property per row.
    """
    
    hold_df = res_props_df.copy()

    for prop in prop_cols:
        hold_df[f'{prop}_norm'] = (hold_df[prop] - hold_df[prop].min()) / (hold_df[prop].max() - hold_df[prop].min())

    hold_df['max_norm'] = hold_df[[f"{p}_norm" for p in prop_cols]].max(axis=1)
    hold_df['max_property'] = hold_df[[f"{p}_norm" for p in prop_cols]].idxmax(axis=1)

    unique_props = hold_df['max_property'].unique()
    colors_map = {prop: color for prop, color in zip(unique_props, cm.tab20.colors)}
    bar_colors = hold_df['max_property'].map(colors_map)

    plt.figure(figsize=figsize)
    plt.bar(hold_df['pos_label'], hold_df['count'], color=bar_colors)
    plt.xticks(rotation=90)
    plt.xlabel('Position + Amino Acids')
    plt.ylabel('Count')
    plt.title('Counts per Position with Amino Acids')
    legend_elements = [Patch(facecolor=colors_map[prop], label=prop) for prop in unique_props]
    plt.legend(handles=legend_elements, title='Max Property', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

    return hold_df

def plot_antigen_property_distribution(plot_df, ag_aa_seq, prop_cols, properties_table):
    """
    Generates a property-colored bar chart using antigen residue identities instead of epitope residues.

    Maps each position in ``plot_df`` to the amino acid at the corresponding index in
    ``ag_aa_seq`` (1-based), then runs the full property pipeline
    (:func:`plot_all_physical_props` followed by :func:`plot_max_normalized_property`).

    Args:
        plot_df (pd.DataFrame): DataFrame from :func:`aggregate_plot_df` with a ``'pos'`` column.
        ag_aa_seq (str): Full antigen amino acid sequence. Position ``p`` is looked up as
            ``ag_aa_seq[p-1]``; positions outside the sequence length are assigned ``None``.
        prop_cols (list of str): Property column names passed to
            :func:`plot_max_normalized_property`.
        properties_table (str or pd.DataFrame): Path to CSV or DataFrame with amino acid
            properties.

    Returns:
        pd.DataFrame: Result of :func:`plot_max_normalized_property` with normalized
        property columns and max-property assignments.
    """
    
    plot_df['aa_original'] = plot_df['pos'].apply(
        lambda p: ag_aa_seq[p-1] if 1 <= p <= len(ag_aa_seq) else None
    )
    res_props_df = plot_all_physical_props(plot_df, aa_col='aa_original', properties_table=properties_table)
    return plot_max_normalized_property(res_props_df, prop_cols, figsize=(50, 6))


### --- PIPELINE RUNNER --- ###

def run_epitope_analysis(epitope_file, ag_aa_seq, properties_table, start_pos=None, end_pos=None):
    """
    Runs the full epitope visualization pipeline from raw CSV to property-colored bar charts.

    Chains all pipeline functions in order:

        1. :func:`generate_epi_df` — parse and filter epitopes.
        2. :func:`aggregate_plot_df` — aggregate counts by position.
        3. :func:`plot_all_physical_props` — color by every available property.
        4. :func:`plot_max_normalized_property` — render final chart with max-normalized
           coloring across a curated subset of properties.

    The curated property subset includes hydrophobicity, hydrophilicity, hydrogen bonds,
    side chain volume, polarity, polarizability, solvent-accessible surface area, net
    charge index, and amino acid mass.

    Args:
        epitope_file (str): Path to the IEDB CSV file.
        ag_aa_seq (str): Full antigen amino acid sequence for residue lookups.
        properties_table (str or pd.DataFrame): Path to CSV or DataFrame with amino acid
            properties.
        start_pos (int): Lower bound of the residue range (inclusive).
        end_pos (int): Upper bound of the residue range (inclusive).

    Returns:
        pd.DataFrame: Result of :func:`plot_max_normalized_property` with normalized
        property columns and max-property assignments.
    """

    prop_cols = ['hydrophobicity','hydrophilicity','h-bonds','side_chain_vol','polarity','polarizability','solv_accessible_surface_area','net_charge_index','aa_mass']

    # 1 Epitope DF
    epi_df = generate_epi_df(epitope_file, start_pos, end_pos)
    # print(f"Filtered epitope count: {len(epi_df)}")
    # print(epi_df)

    # 2 Aggregate for plotting
    plot_df = aggregate_plot_df(epi_df)
    # print(plot_df.sort_values(by='count', ascending=False))

    # 3 Simple count plot
    # plot_counts_bar(plot_df)

    # 4 Plot by all properties
    res_props_df = plot_all_physical_props(plot_df, properties_table=properties_table)

    # 5 Plot by max normalized property (epitope AAs)
    hold_df = plot_max_normalized_property(res_props_df, prop_cols)

    # 6 Plot by max normalized property (antigen AAs)
    plot_antigen_property_distribution(plot_df, ag_aa_seq, prop_cols, properties_table)

    return epi_df, plot_df, res_props_df, hold_df



### FULL antigen/protein sequence from Uniprot ID described in IEDB search. Specify start and end positions of the region of interest in the sequence using start_pos and end_pos variables below
### Example used: SARS-CoV-2 Spike protein (UNIPROT P0DTC2), RBD region (residues 319-541)

# physicochemical_properties_table ='/Users/isaacdaviet/Desktop/comp/iedb/aa_physicochemical_properties.csv'

# epitope_file = '/Users/isaacdaviet/Desktop/comp/iedb/IEDB_BCRepitopes_SarsCov_disc_neut.csv'

# ag_aa_seq = 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'

# start_pos = 319
# end_pos = 541

# epi_df, plot_df, res_props_df, hold_df = run_epitope_analysis(
#     epitope_file=epitope_file,
#     ag_aa_seq=ag_aa_seq,
#     properties_table=physicochemical_properties_table
# ) ### dfs are a bit all over the place. I might clean them up at some point, but if this message is still here, it means the dfs of interestshould probably be generated from scratch using the functions above
