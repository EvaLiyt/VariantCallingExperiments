import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import itertools

# ---- Output directory ----
output_dir = "./confusion_matrix"
os.makedirs(output_dir, exist_ok=True)

# ---- Helper to create a small confusion matrix ----
def create_cm_matrix(TP, FP, TN, FN):
    return np.array([[TP, FN],
                     [FP, TN]])

# ---- Helper to plot mini confusion matrix ----
def plot_cm(ax, mat, n=None):
    pos_vals, neg_vals = mat[0,:], mat[1,:]
    norm_pos = (pos_vals - pos_vals.min()) / max(pos_vals.max() - pos_vals.min(), 1e-6)
    norm_neg = (neg_vals - neg_vals.min()) / max(neg_vals.max() - neg_vals.min(), 1e-6)

    def soft_color(cmap, value, alpha=0.7):
        rgba = cmap(value)
        rgba_soft = [alpha*x + (1-alpha)*1 for x in rgba[:3]] + [rgba[3]]
        return rgba_soft

    for i in range(2):
        for j in range(2):
            if i == 0:  # positive row
                cmap = plt.cm.OrRd if norm_pos[j] < 0.5 else plt.cm.PuRd
                color = soft_color(cmap, norm_pos[j], alpha=0.6)
            else:  # negative row
                cmap = plt.cm.BuGn if norm_neg[j] < 0.5 else plt.cm.YlGnBu
                color = soft_color(cmap, norm_neg[j], alpha=0.6)
    
            ax.add_patch(Rectangle((j, 1-i), 1, 1, facecolor=color, edgecolor='black'))
            ax.text(j+0.5, 1-i+0.5, f"{mat[i,j]:.2f}" if mat[i,j] < 10 else f"{int(mat[i,j])}",
                    ha='center', va='center', fontsize=6)

    ax.text(-0.17, 1, "Actual", ha='right', va='center', rotation=90, fontsize=3.5) 
    ax.text(-0.05, 1.5, "Positive", ha='right', va='center', rotation=90, fontsize=3) 
    ax.text(-0.05, 0.5, "Negative", ha='right', va='center', rotation=90, fontsize=3) 
    ax.text(0.5, 2.2, "Predicted", ha='center', va='bottom', fontsize=3.5) 
    ax.text(0.5, 2, "Positive", ha='center', va='bottom', fontsize=3) 
    ax.text(1.5, 2, "Negative", ha='center', va='bottom', fontsize=3) 
    if n is not None: 
        ax.text(2, 2.3, f"n={n}", va='center', ha='right', fontsize=3.5)

    ax.set_xlim(0,2)
    ax.set_ylim(0,2)
    ax.axis('off')

# ---- Helper to plot F1 score ----
def plot_f1_from_csv(ax, f1, n=None, factor=0.6):
    rgba = plt.cm.viridis(f1)
    rgba_light = [1 - factor*(1-x) for x in rgba[:3]] + [rgba[3]]

    ax.add_patch(Rectangle((0,0),1,1, facecolor=rgba_light, edgecolor='black'))
    ax.text(0.5, 0.5, f"{f1:.2f}", ha='center', va='center', fontsize=6)
    if n is not None:
        ax.text(1.05, 1.05, f"n={n}", ha='right', va='bottom', fontsize=4)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.axis('off')

# ---- General function to plot confusion matrix scatter ----
def scatter_cm(input_file, model_label="Model", count_type="number"):
    df = pd.read_csv(input_file)
    genotypes = [c.replace("_TPR","") for c in df.columns if c.endswith("_TPR")]

    # ---- Generate custom minD/minVar combinations ----
    all_combinations = []
    for minD in range(1,9):           # minD 1-8
        for minVar in range(1,8):     # minVar 1-7
            if minD == 1 and minVar == 1:
                all_combinations.append((minD, minVar))
            elif minD > 1 and minVar < minD:
                all_combinations.append((minD, minVar))
    all_combinations = pd.DataFrame(all_combinations, columns=['minD_value','minVar_value'])

    for g in genotypes:
        if count_type == "number":
            TP_col, FP_col, TN_col, FN_col = f"{g}_TP", f"{g}_FP", f"{g}_TN", f"{g}_FN"
        else:
            TP_col, FP_col, TN_col, FN_col = f"{g}_TPR", f"{g}_FPR", f"{g}_TNR", f"{g}_FNR"

        grouped = df.groupby(['minD_value','minVar_value']).agg(
            TP=pd.NamedAgg(column=TP_col, aggfunc='mean'),
            FP=pd.NamedAgg(column=FP_col, aggfunc='mean'),
            TN=pd.NamedAgg(column=TN_col, aggfunc='mean'),
            FN=pd.NamedAgg(column=FN_col, aggfunc='mean'),
            n=pd.NamedAgg(column=TP_col, aggfunc='count')
        ).reset_index()

        # ---- Merge with custom combinations ----
        grouped = all_combinations.merge(grouped, on=['minD_value','minVar_value'], how='left')
        grouped[['TP','FP','TN','FN','n']] = grouped[['TP','FP','TN','FN','n']].fillna(0)

        fig, ax = plt.subplots(figsize=(8,7))
        ax.set_xlim(0, 9)
        ax.set_ylim(0, 9)
        ax.set_xlabel("minD value")
        ax.set_ylabel("minVar value")
        ax.set_title(f"{model_label} – Genotype {g} Confusion Matrices ({count_type})")

        width = 0.8
        height = 0.8
        for _, row in grouped.iterrows():
            ax_inset = ax.inset_axes([
                (row['minD_value'] - width/2) / 9,
                (row['minVar_value'] - height/2) / 9,
                width / 9,
                height / 9
            ])
            cm = create_cm_matrix(row['TP'], row['FP'], row['TN'], row['FN'])
            plot_cm(ax_inset, cm, n=int(row['n']))

        # Positive values (orange–purple)
        norm_pos = Normalize(vmin=0, vmax=1)
        sm_pos = ScalarMappable(norm=norm_pos, cmap=plt.cm.PuRd)
        cbar_pos = fig.colorbar(sm_pos, ax=ax, fraction=0.025, pad=0.02)
        cbar_pos.set_label('Positive scale (TP/FN)', fontsize=8)
        cbar_pos.ax.tick_params(labelsize=6)

        # Negative values (blue–yellow)
        norm_neg = Normalize(vmin=0, vmax=1)
        sm_neg = ScalarMappable(norm=norm_neg, cmap=plt.cm.YlGnBu)
        cbar_neg = fig.colorbar(sm_neg, ax=ax, fraction=0.025, pad=0.08)
        cbar_neg.set_label('Negative scale (FP/TN)', fontsize=8)
        cbar_neg.ax.tick_params(labelsize=6)
        
        plt.savefig(os.path.join(output_dir, f"cm_{count_type}_{model_label.lower().replace(' ','_')}_{g}.png"), dpi=300)
        plt.close()
        print(f"Saved CM scatter plot for genotype {g}")

        # ---- F1 plot for rate plots ----
        if count_type == "rate" and g == "Mean" and "Mean_F1Score" in df.columns:
            grouped_f1 = df.groupby(['minD_value','minVar_value']).agg(
                Mean_F1=pd.NamedAgg(column="Mean_F1Score", aggfunc='mean'),
                n=pd.NamedAgg(column="Mean_F1Score", aggfunc='count')
            ).reset_index()

            grouped_f1 = all_combinations.merge(grouped_f1, on=['minD_value','minVar_value'], how='left')
            grouped_f1[['Mean_F1','n']] = grouped_f1[['Mean_F1','n']].fillna(0)

            fig, ax = plt.subplots(figsize=(8,7))
            ax.set_xlim(0, 9)
            ax.set_ylim(0, 9)
            ax.set_xlabel("minD value")
            ax.set_ylabel("minVar value")
            ax.set_title(f"{model_label} – Genotype {g} F1 Scores ({count_type})")

            for _, row in grouped_f1.iterrows():
                ax_inset = ax.inset_axes([
                    (row['minD_value'] - width/2) / 9,
                    (row['minVar_value'] - height/2) / 9,
                    width / 9,
                    height / 9
                ])
                plot_f1_from_csv(ax_inset, row['Mean_F1'], n=int(row['n']))

            # ---- Add F1 colorbar ----
          
            norm_f1 = Normalize(vmin=0, vmax=1)
            sm_f1 = ScalarMappable(norm=norm_f1, cmap=plt.cm.viridis)
            cbar_f1 = fig.colorbar(sm_f1, ax=ax, fraction=0.025, pad=0.02)
            cbar_f1.set_label('Mean F1 Score', fontsize=8)
            cbar_f1.ax.tick_params(labelsize=6)

            
            plt.savefig(os.path.join(output_dir, f"f1_{count_type}_{model_label.lower().replace(' ','_')}_{g}.png"), dpi=300)
            plt.close()
            print(f"Saved F1 scatter plot for genotype {g}")

# ---- Run for multiple models ----
for model, path in [("Complex Model", "../realisticModel/4-analysis/summary/summary.csv"),
                    ("Calibrated Model", "../calibratedModel/4-analysis/summary/summary.csv"),
                    ("Simple Model", "../simpleModel/4-analysis/summary/summary.csv")]:
    scatter_cm(path, model_label=model, count_type="number")
    scatter_cm(path, model_label=model, count_type="rate")
