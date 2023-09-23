import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

#https://github.com/zhanglab-aim/CROTONdb
# Upload code to GH

MEDIUM_SIZE = 14

###########################
## Filter PAM-disrupters ##
###########################

# we should still filter out PAM-disrupting variants, as repair happens after 
# cleavage, and PAM-disrupting variants most likely won’t be cleaved;

def filter_PAM_distrupt_vars(gene_df):
    # Input: CROTONdb downloaded dataframe
    # Output: same dataframe (retains cols) but filtered

    # gene_df["Ref. Sequence"] = gene_df["Ref. and Alt. Sequences"].str[0:60]
    # gene_df["Alt. Sequence"] = gene_df["Ref. and Alt. Sequences"].str[60:]
    # gene_df["Ref. PAM"] = gene_df["Ref. and Alt. Sequences"].str[33:36]
    # gene_df["Alt. PAM"] = gene_df["Ref. and Alt. Sequences"].str[93:96]
    print("All variants", gene_df.shape[0])
    gene_df = gene_df[gene_df["Ref. and Alt. Sequences"].str[94:96] == "GG"]\
                                                        .reset_index(drop = True)
    print("Filtered variants", gene_df.shape[0])

    # gene_df.to_csv("hi.csv")
    return gene_df

##########
## CD33 ##
##########

# ------------------------------------------------------------------------------
# ? Plot histogram of max variant effects 
# ? (over the whole gene and over the given PAM site –– CD33|17)
# GCTGCAAGTGCAGGAGTCAGTGACGGTACA|GGAGGGTTTGTGCGTCCTCGTGCCCTGCAC
# gRNA4: GAGTCAGTGACGGTACAGGA

def CD33_maxvar_hist(cd33_df):
    plt.figure(figsize = (7, 8))
    plt.ylabel("Density", fontsize = MEDIUM_SIZE)
    plt.xlabel("Maximum Variant Effect (%)", fontsize = MEDIUM_SIZE)
    plt.tick_params(labelsize = MEDIUM_SIZE - 2)

    overall_patch = mpatches.Patch(color = "tab:blue", alpha = 0.6, label = "Overall")
    target_patch = mpatches.Patch(color = "tab:orange", alpha = 0.6, label = "Targeted\nPAM")
    plt.legend(handles=[overall_patch, target_patch], fontsize = MEDIUM_SIZE)
    
    plt.hist(cd33_df["Max Variant Effect (%)"], density = True, bins = 15, alpha = 0.6)
    plt.hist(cd33_df[cd33_df["PAM ID"] == "CD33|17"]["Max Variant Effect (%)"], 
            density = True, bins = 15, alpha = 0.5)
    
    plt.savefig("images/CD33_maxvar_hist.png", bbox_inches = "tight", dpi = 350)

# ------------------------------------------------------------------------------
# ? Specific variant talking about with highest variant effect
# ? Three sets of bar plots––one for reference, one for alternative allele 
# “rs1339188502,” a missense A to G mutation

def CD33_refvsalt_bars(cd33_df):
    highest_var = cd33_df[(cd33_df["Variant ID"] == "rs1339188502") \
                    & (cd33_df["PAM ID"] == "CD33|17")] 
    highest_var_dict = {"Predicted Value": [], "Allele": ["Reference"] * 3 + ["Alternate\nrs1339188502"] * 3, 
                        "Editing Outcome": ["1 bp Insertion", "1 bp Frameshift", "2 bp Frameshift"] * 2}

    highest_var_dict["Predicted Value"].extend(highest_var[["Ref. 1bp Insertion (%)", 
                    "Ref. 1bp Frameshift (%)", "Ref. 2bp Frameshift (%)", 
                    "Alt. 1bp Insertion (%)", "Alt. 1bp Frameshift (%)", 
                    "Alt. 2bp Frameshift (%)"]].values.tolist()[0])
    
    plt_df = pd.DataFrame.from_dict(highest_var_dict)

    plt.figure(figsize = (10, 8))
    ax = sns.barplot(data = plt_df, x = "Editing Outcome", y = "Predicted Value", hue = "Allele")
    # ax = [above] AND THEN # show_values(ax)
    ax.set_ylabel("Predicted Value", fontsize = MEDIUM_SIZE)
    ax.set_xlabel("Editing Oucome", fontsize = MEDIUM_SIZE)
    ax.tick_params(labelsize = MEDIUM_SIZE - 2)
    ax.legend(title = "Allele", fontsize = MEDIUM_SIZE)
    plt.setp(ax.get_legend().get_title(), fontsize= MEDIUM_SIZE) # for legend title

    plt.savefig("images/CD33_refvsalt_bars.png", bbox_inches = "tight", dpi = 350)

# ------------------------------------------------------------------------------
# ? Since already comparing gene, histogram of SNV mutations
# (nucleotide substitution frequency map for CD33)

def CD33_SNV_hist(cd33_df):
    CD33allele_cts = cd33_df.value_counts(subset = ["Ref. Allele", "Alt. Allele"]) \
                            .reset_index().rename(columns = {0:"Count"})
    allele_dict = {"A": [], "T": [], "C": [], "G": []}

    x = list(allele_dict.keys())
    for ref in x: 
        for alt in x:
            allele_ct = (CD33allele_cts[(CD33allele_cts["Ref. Allele"] == ref) \
                        & (CD33allele_cts["Alt. Allele"] == alt)]["Count"].tolist())
            if allele_ct == []: 
                allele_dict[ref].append(0)
            else: 
                allele_dict[ref].extend(allele_ct)

    plt.figure(figsize = (11, 8))
    plt.ylabel("Frequency", fontsize = MEDIUM_SIZE)
    plt.xlabel("Reference Allele", fontsize = MEDIUM_SIZE)
    plt.tick_params(labelsize = MEDIUM_SIZE - 2)

    # Plot stacked bar chart
    y1 = np.array(allele_dict["A"])
    y2 = np.array(allele_dict["T"])
    y3 = np.array(allele_dict["C"])
    y4 = np.array(allele_dict["G"])
    plt.bar(x, y1)
    plt.bar(x, y2, bottom = y1)
    plt.bar(x, y3, bottom = y1 + y2)
    plt.bar(x, y4, bottom = y1 + y2 + y3)
    legend = plt.legend(x, title = "Alternate\nAllele", fontsize = MEDIUM_SIZE)

    # legend = plt.legend(title = "Alternate Allele", fontsize = MEDIUM_SIZE)
    plt.setp(legend.get_title(), fontsize = MEDIUM_SIZE)
  
    
    plt.savefig("images/CD33_SNV_hist.png", bbox_inches = "tight", dpi = 350)

# ------------------------------------------------------------------------------
# ? They have identified other targets look at other targets 
# ? and see if they’re better in terms of variant effect

def get_altgRNAs(cd33_df): 
    # FROM PAPER (FULL LENGTH VERSIONS): 
    # gRNA1 tggggtgattatgagcaccg	
    # gRNA2 tgagcatcgtagacgccagg	
    # gRNA3 atccctggcactctagaacc	
    # gRNA4 gagtcagtgacggtacagga ** THE CHOSEN ONE -->
    # gRNA5 tgtcacatgcacagagagct	

    # FROM USING GENOME BROWSER: 
    # gRNA1 ggtcctggggccgTGGGGTGATTATGAGCA|CCGaggagtgagtagtcctggggcccaggg
    # gRNA2 aacaactgctcccTGAGCATCGTAGACGCC|AGGaggagggataatcgttcatacttcttt
    # gRNA3 gcccaaaatcctcATCCCTGGCACTCTAGA|ACCcggccactccaaaaacctgacctgctc
    # gRNA4 gctgcaagtgcagGAGTCAGTGACGGTACA|GGAgggtttgtgcgtcctcgtgccctgcac
    # gRNA5 cctgtgcctcaccTGTCACATGCACAGAGA|GCTggggagatttgtaactgtgtttggtac
    
    cd33_df["Ref. Frameshift"] = cd33_df["Ref. 1bp Frameshift (%)"] + cd33_df["Ref. 2bp Frameshift (%)"]
    cd33_df["Alt. Frameshift"] = cd33_df["Alt. 1bp Frameshift (%)"] + cd33_df["Alt. 2bp Frameshift (%)"]
    cd33_df = cd33_df[["Variant ID", "PAM ID", "Max Variant Effect (%)", 
                       "Ref. Frameshift", "Alt. Frameshift", ]]

    # PAM IDs associated with gRNA1, gRNA2, gRNA3, gRNA4, gRNA5 in that order
    PAMid_lst = ["CD33|112", "CD33|58", "CD33|79", "CD33|17", "CD33|68"]
    fs_maxvar_lst = []
    for PAMid in PAMid_lst: 
        gRNA_df =  cd33_df[cd33_df["PAM ID"] == PAMid]
        fs = gRNA_df["Ref. Frameshift"].tolist()[0]
        maxvar = gRNA_df["Max Variant Effect (%)"].max()
        fs_maxvar_lst.append((fs, maxvar))

    return fs_maxvar_lst

# ------------------------------------------------------------------------------
def get_CD33_plts():
    plt.figure(figsize = (7, 8))
    cd33_df = filter_PAM_distrupt_vars(pd.read_csv("data/CD33.csv")).drop_duplicates(subset=["PAM ID"])
    plt.hist(cd33_df["Ref. 1bp Frameshift (%)"] + cd33_df["Ref. 2bp Frameshift (%)"], 
             density = True, bins = 15, alpha = 0.6)
    plt.ylabel("Density", fontsize = MEDIUM_SIZE)
    plt.xlabel("Frameshift Frequency (%)", fontsize = MEDIUM_SIZE)
    plt.tick_params(labelsize = MEDIUM_SIZE - 2)
    plt.savefig("images/CD33_fs_dist.png", bbox_inches = "tight", dpi = 350)

    # cd33_df = filter_PAM_distrupt_vars(pd.read_csv("data/CD33.csv"))
    # print(plt.rcParams['font.family'])
    # CD33_refvsalt_bars(cd33_df)
    
    # CD33_maxvar_hist(cd33_df) 
    # CD33_SNV_hist(cd33_df)
    # CD33_refvsalt_bars(cd33_df)
    # print(get_altgRNAs(pd.read_csv("data/CD33.csv")))
    # print(get_altgRNAs(cd33_df))
get_CD33_plts()
##########
## TTR ##
##########

# ? Plot TTR version of PDCD1 figure with frameshift frequency
# Note: summing 1 bp frameshift and 2 bp frameshift

def get_overlap(x1, x2, y1, y2):
    # INPUT: ranges x1, x2 and y1, y2 where x1 < x2 and y1 < y2
    # First, check if there is overlap between ranges
    # Then, compute and return overlap if it exists, ow return 0
    if max(x1, y1) <= min(x2, y2): 
        return min(x2, y2) - max(x1, y1) + 1 
    else:
        return 0

def TTR_fs_genewide(ttr_df): 
    # NOTE: prev version: https://github.com/zj-zhang/CROTON-dev/blob/master/src/variant/plt_variant.py

    ttr_df[["Gene", "PAM Index"]] = ttr_df["PAM ID"].str.split("|", expand=True)
    ttr_df["Ref. Frameshift"] = ttr_df["Ref. 1bp Frameshift (%)"] + ttr_df["Ref. 2bp Frameshift (%)"]
    ttr_df["Alt. Frameshift"] = ttr_df["Alt. 1bp Frameshift (%)"] + ttr_df["Alt. 2bp Frameshift (%)"]
    ttr_df[["strand", ":", "start", "-", "end"]] = ttr_df["PAM Range"].str.split(expand = True)
    plt_df = ttr_df[["PAM Index", "Variant ID", "Variant Position", "strand", "start", 
        "end", "Ref. Frameshift", "Alt. Frameshift", "Ref. and Alt. Sequences"]].reset_index(drop = True)
    
    exons_lst = []

    # Match to exons
    for i in range(len(plt_df)): # len(plt_df)
        start = int(plt_df.loc[i, "start"]) 
        end = int(plt_df.loc[i, "end"])

        # based on TTR_exons.csv --> see file
        # Based on: exons_df = pd.read_csv("data/TTR_exons.csv")
        e1_start, e1_end = 31591903, 31591971 # CDS
        e2_start, e2_end = 31592896, 31593026 # exon, CDS
        e3_start, e3_end = 31595120, 31595255 # exon, CDS
        e4_start, e4_end = 31598568, 31598821 # exon
        cutoff = 30

        # print(start_coord, end_coord)
        if get_overlap(start, end, e1_start, e1_end) >= cutoff:
            exons_lst.append("E1")
        elif get_overlap(start, end, e2_start, e2_end) >= cutoff:
            exons_lst.append("E2")
        elif get_overlap(start, end, e3_start, e3_end) >= cutoff: 
            exons_lst.append("E3")
        elif get_overlap(start, end, e4_start, e4_end) >= cutoff: 
            exons_lst.append("E4")
        else:
            exons_lst.append(np.nan)
    
    plt_df["Exonic Region"] = exons_lst
    plt_df = plt_df.dropna()
    plt_df["PAM Index"] = plt_df["PAM Index"].astype(int)
    plt_df = plt_df.sort_values(by = "PAM Index")

    # Create and stylize plot
    plt.figure(figsize = (16, 6))
    sns.set(style = "whitegrid")
    sns.boxplot(x = "PAM Index", y = "Ref. Frameshift", data = plt_df, 
                    showfliers = False, linewidth = 3, color = 'black')
    ax = sns.stripplot(x = "PAM Index", y = "Alt. Frameshift", hue = "Exonic Region", 
                    data = plt_df, marker = 'o', linewidth = 0.5, size = 6) 
    plt.show()
    # for spine in ["right", "top", "left", "bottom"]:
    #     ax.spines[spine].set_visible(False)
    # plt.ylabel("Frameshift Frequency", fontsize = MEDIUM_SIZE)
    # plt.setp(ax.get_xticklabels(), visible = False)
    # plt.tick_params(labelsize = MEDIUM_SIZE - 2)
    # legend = plt.legend(title = "Exonic Region", bbox_to_anchor=(1.01, 1), 
    #          fontsize = MEDIUM_SIZE, ncol = 2, columnspacing = 0.6, 
    #          borderaxespad = 0., handletextpad = 0.1)
    # plt.setp(legend.get_title(), fontsize = MEDIUM_SIZE)
    # plt.savefig('images/TTR_fs_genewide.png', bbox_inches='tight', dpi = 350)

    # plt.close()

# ------------------------------------------------------------------------------
# ? Zoom into predisposed variant
# rs76992529 variant effects

def TTR_maxvar_hist(trr_df):
    plt.figure(figsize = (7, 8))
    plt.ylabel("Density", fontsize = MEDIUM_SIZE)
    plt.xlabel("Maximum Variant Effect (%)", fontsize = MEDIUM_SIZE)
    plt.tick_params(labelsize = MEDIUM_SIZE - 2)

    overall_patch = mpatches.Patch(color = "tab:blue", alpha = 0.6, label = "Overall")
    var_patch = mpatches.Patch(color = "tab:orange", alpha = 0.6, label = "Variant\nrs76992529")
    plt.legend(handles=[overall_patch, var_patch], fontsize = MEDIUM_SIZE)
    
    plt.hist(trr_df["Max Variant Effect (%)"], density = True, bins = 30, alpha = 0.6)
    plt.hist(trr_df[trr_df["Variant ID"] == "rs76992529"]["Max Variant Effect (%)"], 
            density = True, bins = 2, alpha = 0.5) 

    plt.savefig('images/TTR_maxvar_hist.png', bbox_inches='tight', dpi=350)

# ------------------------------------------------------------------------------
# ? Scatter plot of max variant effect vs frameshift
# Visualizing the same data from fig1D for TTR, and we can see if there is 
# a Pareto front (https://en.wikipedia.org/wiki/Pareto_front)

def TTR_maxvar_fs(ttr_df): 
    ttr_df["Ref. Frameshift"] = ttr_df["Ref. 1bp Frameshift (%)"] + ttr_df["Ref. 2bp Frameshift (%)"]
    ttr_df["Alt. Frameshift"] = ttr_df["Alt. 1bp Frameshift (%)"] + ttr_df["Alt. 2bp Frameshift (%)"]
    plt.scatter(ttr_df["Ref. Frameshift"], ttr_df["Max Variant Effect (%)"])
    
    plt.xlabel("Reference Frameshift (%)", fontsize = MEDIUM_SIZE)
    plt.ylabel("Maximum Variant Effect (%)", fontsize = MEDIUM_SIZE)
    plt.tick_params(labelsize = MEDIUM_SIZE - 2)
    
    plt.show()

def better_TTR_candidates(ttr_df):
    "TTR|24"
    
    print(ttr_df) 

# ------------------------------------------------------------------------------
def get_TTR_plts():
    ttr_df = filter_PAM_distrupt_vars(pd.read_csv("data/TTR.csv"))
    # TTR_fs_genewide(ttr_df)
    # TTR_maxvar_hist(ttr_df)
    # TTR_maxvar_fs(ttr_df)
    better_TTR_candidates(ttr_df)

# get_TTR_plts()
# ------------------------------------------------------------------------------
# get_CD33_plts()
# get_TTR_plts()

# def show_values(axs, orient="v", space=.01):
    # def _single(ax):
    #     if orient == "v":
    #         for p in ax.patches:
    #             _x = p.get_x() + p.get_width() / 2
    #             _y = p.get_y() + p.get_height() + (p.get_height()*0.01)
    #             value = '{:.1f}'.format(p.get_height())
    #             ax.text(_x, _y, value, ha="center") 
    #     elif orient == "h":
    #         for p in ax.patches:
    #             _x = p.get_x() + p.get_width() + float(space)
    #             _y = p.get_y() + p.get_height() - (p.get_height()*0.5)
    #             value = '{:.1f}'.format(p.get_width())
    #             ax.text(_x, _y, value, ha="left")

    # if isinstance(axs, np.ndarray):
    #     for idx, ax in np.ndenumerate(axs):
    #         _single(ax)
    # else:
    #     _single(axs)