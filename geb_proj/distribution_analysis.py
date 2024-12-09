import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests



def plot_corr_distribution(astrocyte_human, astrocyte_mouse, vip_human, vip_mouse):
    labels = ["Astrocyte_Human", "Astrocyte_Mouse", "VIP_Human", "VIP_Mouse"]
    data = np.concatenate([astrocyte_human.to_numpy().reshape(-1,1), astrocyte_mouse.to_numpy().reshape(-1,1), vip_human.to_numpy().reshape(-1,1), vip_mouse.to_numpy().reshape(-1,1)], axis = 1)

    plt.hist(data, alpha = 0.7, label = labels)
    plt.legend()
    plt.xlabel("Adjusted Correlation")
    plt.ylabel("Frequency")
    plt.title("Distribution of Adjusted Correlation")
    plt.savefig("adj_corr_distribution.png")
    plt.clf()
    return


def plot_p_value_distribution(astrocyte_human, astrocyte_mouse, vip_human, vip_mouse):
    labels = ["Astrocyte_Human", "Astrocyte_Mouse", "VIP_Human", "VIP_Mouse"]
    data = np.concatenate([astrocyte_human.to_numpy().reshape(-1,1), astrocyte_mouse.to_numpy().reshape(-1,1), vip_human.to_numpy().reshape(-1,1), vip_mouse.to_numpy().reshape(-1,1)], axis = 1)

    plt.hist(data, alpha = 0.7, label = labels)
    plt.legend()
    plt.xlabel("P-value")
    plt.ylabel("Frequency")
    plt.title("Distribution of P-values")
    plt.savefig("p_value_distribution.png")
    plt.clf()
    return

def compute_rank_sum(astrocyte_human, astrocyte_mouse, vip_human, vip_mouse):
    print(len(astrocyte_human), len(astrocyte_mouse), len(vip_human), len(vip_mouse))
    mu_u = len(astrocyte_human) * len(astrocyte_mouse)/2
    n_total = len(astrocyte_human) + len(astrocyte_mouse)
    sigma_u = np.sqrt(len(astrocyte_human) * len(astrocyte_mouse) * (len(astrocyte_human) + len(astrocyte_mouse) + 1)/12)


    stat_1, p_1 = mannwhitneyu(astrocyte_human, astrocyte_mouse)#mannwhitneyu(astrocyte_human, astrocyte_mouse)
    print(((stat_1 - mu_u)/sigma_u)/(np.sqrt(n_total)))
    stat_2, p_2 = mannwhitneyu(astrocyte_human, vip_human) #mannwhitneyu(astrocyte_human, vip_human)
    print(((stat_2 - mu_u)/sigma_u)/(np.sqrt(n_total)))
    stat_3, p_3 = mannwhitneyu(astrocyte_human, vip_mouse)
    
    stat_4, p_4 = mannwhitneyu(astrocyte_mouse, vip_human)
    
    stat_5, p_5 = mannwhitneyu(astrocyte_mouse, vip_mouse)
    print(((stat_5 - mu_u)/sigma_u)/(np.sqrt(n_total)))
    stat_6, p_6 = mannwhitneyu(vip_human, vip_mouse)
    print(((stat_6 - mu_u)/sigma_u)/(np.sqrt(n_total)))
    p_vals = [p_1,p_2,p_3,p_4,p_5,p_6]
    adjusted_p = multipletests(p_vals, method='fdr_bh')[1]
    print(f"Adjusted p-values: {adjusted_p}")
    return
    


if __name__ == "__main__":
    astrocyte_human = pd.read_csv("Bridges Files/Peak_Phylolm/peakPhyloResults_astrocyte_human_new.csv")
    astrocyte_mouse = pd.read_csv("Bridges Files/Peak_Phylolm/peakPhyloResults_astrocyte_mouse_new.csv")
    vip_human = pd.read_csv("Bridges Files/Peak_Phylolm/peakPhyloResults_vip_human_new.csv")
    vip_mouse = pd.read_csv("Bridges Files/Peak_Phylolm/peakPhyloResults_vip_mouse_new.csv")

    #plot_p_value_distribution(astrocyte_human['pvalue'], astrocyte_mouse['pvalue'], vip_human['pvalue'], vip_mouse['pvalue'])
    #plot_corr_distribution(astrocyte_human['adjCorrelation'], astrocyte_mouse['adjCorrelation'], vip_human['adjCorrelation'], vip_mouse['adjCorrelation'])

    #compute_rank_sum(astrocyte_human['pvalue'], astrocyte_mouse['pvalue'], vip_human['pvalue'], vip_mouse['pvalue'])
    compute_rank_sum(astrocyte_human['adjCorrelation'], astrocyte_mouse['adjCorrelation'], vip_human['adjCorrelation'], vip_mouse['adjCorrelation'])