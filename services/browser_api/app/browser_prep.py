import pandas as pd
import numpy as np
import plotly.graph_objects as go
import subprocess
import json
import os


### will need to have create_master_keys function when data structure is figured out for master keys

def create_master(rel, master, out_path):
    # Simplify columns
    master.rename(columns = {'GP2ID': 'IID', 'age_at_sample_collection': 'age', 'baseline_GP2_phenotype':'pheno', 'biological_sex_for_qc': 'sex'}, inplace = True)
    clean_key = master[['IID', 'sex', 'pheno', 'age', 'study', 'nba_prune_reason', 'nba_related', 'nba_label', 'nba', 'wgs_prune_reason', 'wgs_label', 'wgs']]
    clean_key['release'] = rel

    # Simplify pruned_reason
    clean_key.nba_prune_reason = clean_key.nba_prune_reason.str.split('-').str[0]
    # clean_key.wgs_prune_reason = clean_key.wgs_prune_reason.str.split('-').str[0]

    # Remove trailing spaces
    clean_key.nba_prune_reason = clean_key.nba_prune_reason.str.split(' ').str[0]
    # clean_key.wgs_prune_reason = clean_key.wgs_prune_reason.str.split(' ').str[0]

    # Split into NBA and WGS keys
    # wgs = clean_key[(clean_key.wgs == 1)]

    clean_nba = clean_key[['IID','sex','pheno','age', 'study', 'nba_prune_reason','nba_related','nba_label', 'release']]
    # clean_wgs = wgs[['IID','sex','pheno','age', 'study', 'wgs_prune_reason','wgs_label', 'release']]

    # Make consistent cols
    clean_nba.rename(columns = {'nba_prune_reason': 'prune_reason', 'nba_related': 'related', 'nba_label': 'label'}, inplace = True)
    # clean_wgs.rename(columns = {'wgs_prune_reason': 'prune_reason', 'wgs_label': 'label'}, inplace = True)

    # Save files
    key_path = 'nba_app_key.csv'
    clean_nba.to_csv(f'{out_path}/{key_path}', index = False)
    # clean_wgs.to_csv(f'{out_path}/wgs_app_key.csv', index = False)

def plot_funnel(funnel_df):
    funnel_plot = go.Figure(go.Funnelarea(
        text=[f'<b>{i}</b>' for i in funnel_df['step_name']],
        values=funnel_df['remaining_samples'],
        marker={
            "colors": [
                   "#999999",  # Gray  
                    "#E69F00",  # Orange  
                    "#56B4E9",  # Sky Blue  
                    "#009E73",  # Green  
                    "#AA4499",  # Purple  
                    "#F0E442",  # Yellow  
                    "#0072B2",  # Blue  
                    "#D55E00",  # Dark Orange  
                    "#CC79A7",  # Pink  
                    "#882255"   # Deep Red  
            ]
        },
        opacity=1.0, textinfo='text',
        customdata=funnel_df['remaining_samples'],
        hovertemplate='Remaining Samples:<br>%{customdata[0]:.f}'+'<extra></extra>'
    ))

    funnel_plot.update_layout(showlegend=False, margin=dict(l=0, r=0, t=10, b=10))

    return funnel_plot

def plot_variants(variant_df):
    variant_plot = go.Figure()
    for col, color in zip(
        variant_df.columns, 
        ["#0072B2", "#882255", "#44AA99", "#D55E00"]):
        variant_plot.add_trace(
            go.Bar(
                x=variant_df.index, 
                y=variant_df[col], 
                name=col, 
                marker_color=color
            )
        )
    variant_plot.update_layout(
        barmode='stack', xaxis=dict(title='Ancestry', tickfont_size=14),
        yaxis=dict(title='Count', tickfont_size=14), # may want to change title font to 16
        width=1100, height=600
    )

    return variant_plot

def prune_steps(proj_samples, df_qc, out_path):
    pre_QC_total = proj_samples['IID'].count()

    sample_level = df_qc[df_qc.level == 'sample']
    funnel_df = sample_level.groupby('step', as_index=False)['pruned_count'].sum()
    funnel_df.loc[-1] = {'pruned_count': 0, 'step': 'pre_QC'}

    ordered_prune = {
        'pre_QC': 'Pre-QC',
        'callrate_prune': 'Call Rate Prune',
        'sex_prune': 'Sex Prune',
        'het_prune': 'Heterozygosity Prune',
        'related_prune': 'Duplicated Prune'
    }

    # Convert 'step' to categorical with the defined order and sort
    funnel_df['step'] = pd.Categorical(funnel_df['step'], categories=ordered_prune.keys(), ordered=True)
    funnel_df.sort_values('step', inplace = True)
    funnel_df['remaining_samples'] = pre_QC_total - funnel_df['pruned_count'].cumsum()
    funnel_df['step_name'] = funnel_df['step'].map(ordered_prune)
    funnel_df.reset_index(inplace = True, drop = True)
    # funnel_df.to_csv(f'{out_path}/funnel_plot.csv', index = False)

    funnel_path = 'funnel_plot.html' ### can remove this structure if don't need to modify save location
    funnel_plot = plot_funnel(funnel_df)
    funnel_plot.write_html(f'{out_path}/{funnel_path}')

def related_qc(rel_samples, out_path):
    ancestry_dict = {
                'AFR': 'African', 
                'SAS': 'South Asian', 
                'EAS': 'East Asian', 
                'EUR': 'European',
                'AMR': 'American', 
                'AJ': 'Ashkenazi Jewish', 
                'AAC': 'African American/Afro-Caribbean',
                'CAS': 'Central Asian', 
                'MDE': 'Middle Eastern', 
                'FIN': 'Finnish', 
                'CAH': 'Complex Admixture History'
                }

    dup_df = rel_samples[rel_samples.REL == 'duplicate'] # will include pairs
    dup_df.drop_duplicates(subset = 'IID1', inplace = True)
    dup_df['duplicated_count'] = 1
    dup_summary = dup_df.groupby('ancestry', as_index=False)['duplicated_count'].sum()

    rel_samples['related_count'] = 1
    rel_summary = rel_samples.groupby('ancestry', as_index=False)['related_count'].sum()

    related_plot = pd.merge(rel_summary, dup_summary)
    related_plot['label'] = related_plot['ancestry'].map(ancestry_dict)
    related_plot.sort_values(by = 'ancestry', inplace = True)
    related_plot.drop(columns = 'ancestry', inplace = True)

    related_plot.rename(columns = {'label': 'Ancestry Category',
                                'related_count': 'Related Sample Count',
                                'duplicated_count': 'Duplicated Sample Count'}, inplace = True)
    
    related_path = 'related_plot.csv'
    related_plot.to_csv(f'{out_path}/{related_path}', index = False)

def variant_qc(df_qc, out_path):
    # Make ordered metrics like with sample-level with cleaner labels 
    ordered_metrics = {
        'mis_removed_count': 'Case-Control Missingness Removed Variants', 
        'haplotype_removed_count': 'Haplotype Missingness Removed Variants',
        'hwe_removed_count': 'HWE Removed Variants', 
        'geno_removed_count': 'Genotype Missingness Removed Variants'
    }

    dataframes = []
    for metric in ordered_metrics.keys():
        df_metric = df_qc.query(f"metric == '{metric}'").reset_index(drop=True)
        df_metric.rename(
            columns={'pruned_count': metric}, 
            inplace=True
        )
        df_metric.drop(columns = ['metric', 'step', 'level', 'pass'], inplace = True)
        dataframes.append(df_metric)

    df_merged = dataframes[0]
    for df in dataframes[1:]:
        df_merged = pd.merge(
            df_merged, df, 
            on=['ancestry'],
            how='outer', 
            suffixes=('', '_duplicate')
        )

    df_merged.set_index('ancestry', inplace=True)
    df_merged.columns = df_merged.columns.map(ordered_metrics)
    # save with index
    # df_merged.to_csv(f'{out_path}/variant_plot.csv')

    var_plot_path = 'variant_plot.html'
    variant_plot = plot_variants(df_merged)
    variant_plot.write_html(f'{out_path}/{var_plot_path}')

def ancestry_data(proj_labels, proj_samples, ref_pca, proj_pca, out_path):
    # Ancestry Page - PCA Plots
    proj_pca['label'] = 'Predicted'
    proj_labels.rename(columns={'label': 'Predicted Ancestry', 'count': 'Counts'}, inplace = True)
    proj_samples.rename(columns={'label': 'Predicted Ancestry'}, inplace = True)

    proj_pca = proj_pca.merge(proj_samples, on = ['FID', 'IID'], how = 'inner')

    ref_pca_path = 'ref_pca_plot.csv'
    proj_pca_path = 'proj_pca_plot.csv'
    proj_labels_path = 'anc_summary.csv'

    ref_pca.to_csv(f'{out_path}/{ref_pca_path}', index = False)
    proj_pca.to_csv(f'{out_path}/{proj_pca_path}', index = False)
    proj_labels.to_csv(f'{out_path}/{proj_labels_path}', index = False)

def ancestry_model(confusion_matrix, out_path):
    # Ancestry Page - Model Performance
    tp = np.diag(confusion_matrix)
    col_sum = confusion_matrix.sum(axis=0)
    row_sum = confusion_matrix.sum(axis=1)

    class_recall = np.array(tp / row_sum)
    class_precision = np.array(tp / col_sum)

    balanced_accuracy = np.mean(class_recall)
    margin_of_error = 1.96 * np.sqrt(
        (balanced_accuracy * (1 - balanced_accuracy)) / sum(col_sum)
    )
    precision = np.mean(class_precision)
    f1 = np.mean(2 * ((class_recall * class_precision) / (class_recall + class_precision)))

    metrics = {'Balanced Accuracy': [f'{balanced_accuracy:.3f}'],
                'Margin of Error': [f'{margin_of_error:.3f}'],
                'Precision': [f'{precision:.3f}'],
                'F1 Score': [f'{f1:.3f}']}
    model_metrics = pd.DataFrame.from_dict(metrics)

    metrics_path = 'model_metrics.csv'
    cm_path = 'confusion_matrix.csv'

    model_metrics.to_csv(f'{out_path}/{metrics_path}', index = False)
    confusion_matrix.to_csv(f'{out_path}/{cm_path}', index = False)

def ancestry_breakdown(ref_pca, proj_labels, out_path):
    # Ancestry Page - Pie Charts
    df_ref_ancestry_counts = ref_pca['label'].value_counts(normalize=True).rename_axis('Ancestry Category').reset_index(name='Proportion')
    ref_counts = ref_pca['label'].value_counts().rename_axis('Ancestry Category').reset_index(name='Counts')
    ref_combo = pd.merge(df_ref_ancestry_counts, ref_counts, on='Ancestry Category')
    ref_combo.rename(columns={'Proportion': 'Ref Panel Proportion', 'Counts': 'Ref Panel Counts'}, inplace=True)

    total_count = proj_labels['count'].sum() # not the same as total master key count
    proj_labels['Predicted Proportion'] = proj_labels['count'] / total_count
    proj_labels.rename(columns = {'label': 'Ancestry Category', 'count': 'Predicted Counts'}, inplace = True)

    # ref_combo = ref_combo[['Ancestry Category', 'Ref Panel Counts']]
    ref_combo_cah = pd.DataFrame([['CAH', 'NA']], columns=['Ancestry Category', 'Ref Panel Counts'])
    ref_combo = pd.concat([ref_combo, ref_combo_cah], axis=0)
    pie_table = pd.merge(ref_combo, proj_labels, on='Ancestry Category')
    pie_table.sort_values(by = ['Predicted Counts'], ascending = False, inplace = True)

    pie_path = 'pie_table.csv'
    pie_table.to_csv(f'{out_path}/{pie_path}', index = False)

def prep_browser_files(rel: int, master_key: str, gt_output: str, temp_dir: str) -> list:
    output_file = open(gt_output)
    data = json.load(output_file)

    master = pd.read_csv(master_key)

    proj_labels = pd.DataFrame(data['ancestry_counts'])
    proj_samples = pd.DataFrame(data['ancestry_labels'])
    confusion_matrix = pd.DataFrame(data['confusion_matrix'])
    proj_pca = pd.DataFrame(data['projected_pcs'])
    df_qc = pd.DataFrame(data['QC'])
    ref_pca = pd.DataFrame(data['ref_pcs'])
    rel_samples = pd.DataFrame(data['related_samples'])

    create_master(rel, master, temp_dir)
    prune_steps(proj_samples, df_qc, temp_dir)
    related_qc(rel_samples, temp_dir)
    variant_qc(df_qc, temp_dir)
    ancestry_data(proj_labels.copy(), proj_samples, ref_pca, proj_pca, temp_dir)
    ancestry_model(confusion_matrix, temp_dir)
    ancestry_breakdown(ref_pca, proj_labels.copy(), temp_dir)

    # Get list of all output files
    final_files = [f for f in os.listdir(temp_dir) if os.path.isfile(os.path.join(temp_dir, f))]
    return final_files