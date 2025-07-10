import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from scipy.stats import norm, skew, skewtest
from sklearn.model_selection import train_test_split
import sys
import math


# Input format: python normalize_data.py input_file output_file
def pdf(x):
    mean = np.mean(x)
    std = np.std(x)
    y_out = 1/(std * np.sqrt(2 * np.pi)) * np.exp( - (x - mean)**2 / (2 * std**2))
    return y_out


def plot_bell_curve(data):
    print(max(data))
    print(min(data))
    x = np.array(data)
    y = norm.pdf(x, 0, 1)
    # # To fill in values under the bell-curve
    x_fill = np.array(data)
    y_fill = norm.pdf(x_fill, 0, 1)

    # Plotting the bell-shaped curve
    plt.style.use('fivethirtyeight')
    plt.figure(figsize=(6, 6))
    # plt.plot(x, y, color='black')
    plt.scatter(x, y, marker='o',
                s=25, color='red')
    plt.fill_between(x_fill, y_fill, 0,
                     alpha=0.2, color='blue')
    plt.show()


def min_max(df):
    immunogenicity_name = "norm immunogenicity"
    antigenicity_name = "norm antigenicity"
    allergenicity_name = "norm allergenicity"
    # binding_affinity_name = "norm binding affinity"
    immunogenicity_max = df[immunogenicity_name].max()
    immunogenicity_min = df[immunogenicity_name].min()
    antigenicity_max = df[antigenicity_name].max()
    antigenicity_min = df[antigenicity_name].min()
    allergenicity_max = df[allergenicity_name].max()
    allergenicity_min = df[allergenicity_name].min()
    # ba_max = df[binding_affinity_name].max()
    # ba_min = df[binding_affinity_name].min()
    # norm_immunogenicities = []
    # norm_antigenicities = []
    # norm_bas = []
    for i in range(df.shape[0]):
        immunogenicity_value = df[immunogenicity_name][i]
        antigenicity_value = df[antigenicity_name][i]
        allergenicity_value = df[allergenicity_name][i]
        # ba_value = df[binding_affinity_name][i]
        new_immunogenicity_value = (immunogenicity_value - immunogenicity_min)/(immunogenicity_max - immunogenicity_min) # for CD8 use this
        # new_immunogenicity_value = (immunogenicity_value - immunogenicity_max) / (immunogenicity_min - immunogenicity_max)
        new_antigenicity_value = (antigenicity_value - antigenicity_min)/(antigenicity_max - antigenicity_min)
        new_allergenicity_value = (allergenicity_value - allergenicity_max)/(allergenicity_min - allergenicity_max)
        # new_ba_value = (ba_value - ba_min)/(ba_max - ba_min)
        # new_ba_value = (ba_value - ba_max) / (ba_min - ba_max)
        df.at[i, "norm immunogenicity"] = new_immunogenicity_value
        df.at[i, "norm antigenicity"] = new_antigenicity_value
        df.at[i, "norm allergenicity"] = new_allergenicity_value
        # df.at[i, "norm binding affinity"] = new_ba_value
    return df


def log_fifty(df):
    new_column = df["Binding Affinity"].map(lambda a: 1 - math.log(a, 50000))
    df["log(Binding Affinity)"] = new_column
    return df


def z_score_norm(df):
    count = 0
    n = df.shape[0]
    limit = 5000
    one_sd = int(limit * 0.341)
    neg_one_sd = int(limit * 0.341)
    two_sd = int(limit * 0.136)
    neg_two_sd = int(limit * 0.136)
    three_sd = int(limit * 0.021)
    neg_three_sd = int(limit * 0.021)

    # four_sd = int(n * 0.999)
    immunogenicity_name = "Immunogenicity"
    antigenicity_name = "Antigenicity"
    binding_affinity_name = "Binding Affinity"
    immunogenicity_mean = df[immunogenicity_name].mean()
    immunogenicity_std = df[immunogenicity_name].std()
    antigenicity_mean = df[antigenicity_name].mean()
    antigenicity_std = df[antigenicity_name].std()
    ba_mean = df[binding_affinity_name].mean()
    # ba_mean = 24682.836335282653
    # ba_mean = 5000
    ba_std = df[binding_affinity_name].std()

    norm_immunogenicities = []
    norm_antigenicities = []
    norm_bas = []
    to_drop = []

    one_sd_list = []
    neg_one_sd_list = []
    two_sd_list = []
    neg_two_sd_list = []
    three_sd_list = []
    neg_three_sd_list = []

    for i in range(n):
    # while i <= n:
        if count < limit:
            # print(i)
            immunogenicity_value = df[immunogenicity_name][i]
            antigenicity_value = df[antigenicity_name][i]
            ba_value = df[binding_affinity_name][i]
            new_immunogenicity_value = (immunogenicity_value - immunogenicity_mean)/immunogenicity_std
            new_antigenicity_value = (antigenicity_value - antigenicity_mean)/antigenicity_std
            new_ba_value = (ba_value - ba_mean)/ba_std

            # norm_immunogenicities.append(new_immunogenicity_value)
            # norm_antigenicities.append(new_antigenicity_value)
            # norm_bas.append(new_ba_value)
            #
            # df.at[i, "norm immunogenicity"] = new_immunogenicity_value
            # df.at[i, "norm antigenicity"] = new_antigenicity_value
            # df.at[i, "norm binding affinity"] = new_ba_value

            if 0 <= new_ba_value < 1:
                if len(one_sd_list) < one_sd:
                    # epitopes.append(df["Epitope"][i])
                    # results.append(df["Result"][i])
                    # alleles.append(df["Allele"][i])
                    # tested.append(df["Tested"][i])
                    # responded.append(df["Responded"][i])
                    df.at[i, "norm immunogenicity"] = new_immunogenicity_value
                    df.at[i, "norm antigenicity"] = new_antigenicity_value
                    df.at[i, "norm binding affinity"] = new_ba_value
                    # immunogenicities.append(immunogenicity_value)
                    # antigenicities.append(antigenicity_value)
                    # bas.append(ba_value)
                    # obas.append(df["Binding Affinity"][i])
                    norm_immunogenicities.append(new_immunogenicity_value)
                    norm_antigenicities.append(new_antigenicity_value)
                    norm_bas.append(new_ba_value)
                    # binding_scores.append(df["Binding Score"][i])
                    # potentials.append(df["Potential"][i])
                    one_sd_list.append(i)
                    count += 1
            elif -1 < new_ba_value < 0:
                if len(neg_one_sd_list) < neg_one_sd:
                    # epitopes.append(df["Epitope"][i])
                    # results.append(df["Result"][i])
                    # alleles.append(df["Allele"][i])
                    # tested.append(df["Tested"][i])
                    # responded.append(df["Responded"][i])
                    df.at[i, "norm immunogenicity"] = new_immunogenicity_value
                    df.at[i, "norm antigenicity"] = new_antigenicity_value
                    df.at[i, "norm binding affinity"] = new_ba_value
                    # immunogenicities.append(immunogenicity_value)
                    # antigenicities.append(antigenicity_value)
                    # bas.append(ba_value)
                    # obas.append(df["Binding Affinity"][i])
                    norm_immunogenicities.append(new_immunogenicity_value)
                    norm_antigenicities.append(new_antigenicity_value)
                    norm_bas.append(new_ba_value)
                    # binding_scores.append(df["Binding Score"][i])
                    # potentials.append(df["Potential"][i])
                    neg_one_sd_list.append(i)
                    count += 1
            # elif 1 <= new_ba_value < 2:
            #     if len(two_sd_list) < two_sd:
            #         # epitopes.append(df["Epitope"][i])
            #         # results.append(df["Result"][i])
            #         # alleles.append(df["Allele"][i])
            #         # tested.append(df["Tested"][i])
            #         # responded.append(df["Responded"][i])
            #         df.at[i, "norm immunogenicity"] = new_immunogenicity_value
            #         df.at[i, "norm antigenicity"] = new_antigenicity_value
            #         df.at[i, "norm binding affinity"] = new_ba_value
            #         # immunogenicities.append(immunogenicity_value)
            #         # antigenicities.append(antigenicity_value)
            #         # bas.append(ba_value)
            #         # obas.append(df["Binding Affinity"][i])
            #         norm_immunogenicities.append(new_immunogenicity_value)
            #         norm_antigenicities.append(new_antigenicity_value)
            #         norm_bas.append(new_ba_value)
            #         # binding_scores.append(df["Binding Score"][i])
            #         # potentials.append(df["Potential"][i])
            #         two_sd_list.append(i)
            #         count += 1
            # elif -2 < new_ba_value <= -1:
            #     if len(neg_two_sd_list) < neg_two_sd:
            #         # epitopes.append(df["Epitope"][i])
            #         # results.append(df["Result"][i])
            #         # alleles.append(df["Allele"][i])
            #         # tested.append(df["Tested"][i])
            #         # responded.append(df["Responded"][i])
            #         df.at[i, "norm immunogenicity"] = new_immunogenicity_value
            #         df.at[i, "norm antigenicity"] = new_antigenicity_value
            #         df.at[i, "norm binding affinity"] = new_ba_value
            #         # immunogenicities.append(immunogenicity_value)
            #         # antigenicities.append(antigenicity_value)
            #         # bas.append(ba_value)
            #         # obas.append(df["Binding Affinity"][i])
            #         norm_immunogenicities.append(new_immunogenicity_value)
            #         norm_antigenicities.append(new_antigenicity_value)
            #         norm_bas.append(new_ba_value)
            #         # binding_scores.append(df["Binding Score"][i])
            #         # potentials.append(df["Potential"][i])
            #         neg_two_sd_list.append(i)
            #         count += 1
            # elif 2 <= new_ba_value < 3:
            #     if len(three_sd_list) < three_sd:
            #         # epitopes.append(df["Epitope"][i])
            #         # results.append(df["Result"][i])
            #         # alleles.append(df["Allele"][i])
            #         # tested.append(df["Tested"][i])
            #         # responded.append(df["Responded"][i])
            #         df.at[i, "norm immunogenicity"] = new_immunogenicity_value
            #         df.at[i, "norm antigenicity"] = new_antigenicity_value
            #         df.at[i, "norm binding affinity"] = new_ba_value
            #         # immunogenicities.append(immunogenicity_value)
            #         # antigenicities.append(antigenicity_value)
            #         # bas.append(ba_value)
            #         # obas.append(df["Binding Affinity"][i])
            #         norm_immunogenicities.append(new_immunogenicity_value)
            #         norm_antigenicities.append(new_antigenicity_value)
            #         norm_bas.append(new_ba_value)
            #         # binding_scores.append(df["Binding Score"][i])
            #         # potentials.append(df["Potential"][i])
            #         three_sd_list.append(i)
            #         count += 1
            # elif -3 < new_ba_value <= -2:
            #     if len(neg_three_sd_list) < neg_three_sd:
            #         # epitopes.append(df["Epitope"][i])
            #         # results.append(df["Result"][i])
            #         # alleles.append(df["Allele"][i])
            #         # tested.append(df["Tested"][i])
            #         # responded.append(df["Responded"][i])
            #         df.at[i, "norm immunogenicity"] = new_immunogenicity_value
            #         df.at[i, "norm antigenicity"] = new_antigenicity_value
            #         df.at[i, "norm binding affinity"] = new_ba_value
            #         # immunogenicities.append(immunogenicity_value)
            #         # antigenicities.append(antigenicity_value)
            #         # bas.append(ba_value)
            #         # obas.append(df["Binding Affinity"][i])
            #         norm_immunogenicities.append(new_immunogenicity_value)
            #         norm_antigenicities.append(new_antigenicity_value)
            #         norm_bas.append(new_ba_value)
            #         # binding_scores.append(df["Binding Score"][i])
            #         # potentials.append(df["Potential"][i])
            #         neg_three_sd_list.append(i)
            #         count += 1
            else:
                to_drop.append(i)
    print(len(one_sd_list))
    print(len(neg_one_sd_list))
    # print(len(two_sd_list))
    # print(len(neg_two_sd_list))
    # print(len(three_sd_list))
    # print(len(neg_three_sd_list))
    # print()

    # print(len(norm_immunogenicities))
    # print(len(norm_antigenicities))
    # print(len(norm_bas))

    # df.drop(to_drop, inplace=True)
    df = df[df["norm binding affinity"].notna()].reset_index(drop=True)
    # df = df[(df["Tested"] != 0) & (df["Responded"] != 0)]
    print(len(df[df["Result"] == "Positive"]) + len(df[df["Result"] == "Positive-High"]) + len(df[df["Result"] == "Positive-Intermediate"]) + len(df[df["Result"] == "Positive-Low"]))
    print(len(df[df["Result"] == "Negative"]))

    norm_immunogenicities.sort()
    norm_antigenicities.sort()
    norm_bas.sort()

    # print(norm_immunogenicities)
    # print(norm_antigenicities)
    # print(norm_bas)

    plot_bell_curve(norm_immunogenicities)
    plot_bell_curve(norm_antigenicities)
    plot_bell_curve(norm_bas)

    return df


def z_score_norm_test(df, immunogenicity_name, antigenicity_name, allergenicity_name):
    # immunogenicity_name = "immunogenicity"
    # antigenicity_name = "antigenicity"
    # binding_affinity_name = "binding affinity (nM)"
    # immunogenicity_name = "Immunogenicity"
    # antigenicity_name = "Antigenicity"
    # binding_affinity_name = "Binding Affinity"
    immunogenicity_mean = df[immunogenicity_name].mean()
    immunogenicity_std = df[immunogenicity_name].std()
    antigenicity_mean = df[antigenicity_name].mean()
    antigenicity_std = df[antigenicity_name].std()
    # ba_mean = df[binding_affinity_name].mean()
    # ba_std = df[binding_affinity_name].std()
    allergenicity_mean = df[allergenicity_name].mean()
    allergenicity_std = df[allergenicity_name].std()
    for i in range(df.shape[0]):
        immunogenicity_value = df[immunogenicity_name][i]
        antigenicity_value = df[antigenicity_name][i]
        # ba_value = df[binding_affinity_name][i]
        allergenicity_value = df[allergenicity_name][i]
        new_immunogenicity_value = (immunogenicity_value - immunogenicity_mean) / immunogenicity_std
        new_antigenicity_value = (antigenicity_value - antigenicity_mean) / antigenicity_std
        # new_ba_value = (ba_value - ba_mean) / ba_std
        new_allergenicity_value = (allergenicity_value - allergenicity_mean) / allergenicity_std
        df.at[i, "norm immunogenicity"] = new_immunogenicity_value
        df.at[i, "norm antigenicity"] = new_antigenicity_value
        df.at[i, "norm allergenicity"] = new_allergenicity_value
        # df.at[i, "norm binding affinity"] = new_ba_value
    return df


infile = sys.argv[1]
outfile = sys.argv[2]
# df_prior = pd.read_excel("training_set_wvar4_nepmhci.xlsx")
df_prior = pd.read_excel(infile)
# df_prior = log_fifty(df_prior)
normalized_df = z_score_norm_test(df_prior, "immunogenicity", "antigenicity", "allergenicity")
normalized_df = min_max(normalized_df)
# normalized_df.to_excel("normalized_training_set_cd8.xlsx", index=False, header=True)
normalized_df.to_excel(outfile, index=False, header=True)
# train, test = train_test_split(normalized_df, test_size=0.2)
# train.to_excel("normalized_training_set_wvar4.xlsx", header=True, index=False)
# test.to_excel("normalized_test_set_wvar4.xlsx", header=True, index=False)
# #
# cancer = "CRC"

# cancer_mutations_dict = {"CRC": ["R38H", "R88Q", "G106V", "C420R", "E453Q", "E542K", "E545K", "R1023Q", "M1043I", "H1047R"],
#                          "Meningioma": ["E110K", "I391M", "R108H", "G914R", "N345K", "E453K", "Y165H", "H1047R", "E545K"],
#                          "BC": ["E542K", "E542V", "E545K", "Q546E", "Q546R", "H1047L", "H1047R", "N345K", "E726K", "C420R", "G118D", "E453K", "Q546K", "G1049R", "M1043I", "K111E", "E81K", "E545A", "E545G", "N1044K", "S405P"],
#                          "Endometrial": ["E542K", "E542Q", "E545K", "E545G", "G1007R", "Y1021H", "Y1021C", "A1035V", "M1043I", "H1047Y", "H1047R", "G1050D", "T1052K", "H1065L"],
#                          "Glioblastoma": ["R88E", "E542K", "E545A", "T1025N", "Y1021N", "R88Q", "P298T", "R310C", "T1031G", "V344G", "E453K", "E545K", "Y1021C", "M1043I", "N1044S", "H1047Y", "G1049S"]}
# cancer_mutations_dict = {"CRC": ["R38H", "R88Q", "G106V", "C420R", "E453Q", "E542K", "E545K", "R1023Q", "M1043I", "H1047R"]}
# point_mutants = cancer_mutations_dict[cancer]
# final_out = f"normalized_{cancer}_mutants_mhci.xlsx"
# with pd.ExcelWriter(final_out, engine='openpyxl') as w:
#     for pm in point_mutants:
#         df_var = pd.read_excel(f"/Users/mvsamudrala/CancerVaccine/Peptides/epitopes/Automated_Collection/all_variables_mhci_{cancer}.xlsx", sheet_name=pm)
#         # df_var = log_fifty(df_var)
#         df_var = z_score_norm_test(df_var, "immunogenicity", "antigenicity", "binding affinity (nM)", "allergenicity")
#         df_var = min_max(df_var)
#         df_var.to_excel(w, header=True, index=False, sheet_name=pm)


# df_prior = pd.read_excel("training_set_wvar_mhcii3_fixed.xlsx")
# # df_prior = log_fifty(df_prior)
# normalized_df = z_score_norm_test(df_prior, "Immunogenicity", "Antigenicity", "Binding Affinity", "Allergenicity")
# normalized_df = min_max(normalized_df)
# normalized_df.to_excel("tweaked_training_set_mhcii3.xlsx", index=False, header=True)
# train, test = train_test_split(normalized_df, test_size=0.2)
# train.to_excel("normalized_training_set_wvar_mhcii3.xlsx", header=True, index=False)
# test.to_excel("normalized_test_set_wvar_mhcii3.xlsx", header=True, index=False)

# cancer = "Endometrial"
#
# cancer_mutations_dict = {"CRC": ["R38H", "R88Q", "G106V", "C420R", "E453Q", "E542K", "E545K", "R1023Q", "M1043I", "H1047R"],
#                          "Meningioma": ["E110K", "I391M", "R108H", "G914R", "N345K", "E453K", "Y165H", "H1047R", "E545K"],
#                          "BC": ["E542K", "E542V", "E545K", "Q546E", "Q546R", "H1047L", "H1047R", "N345K", "E726K", "C420R", "G118D", "E453K", "Q546K", "G1049R", "M1043I", "K111E", "E81K", "E545A", "E545G", "N1044K", "S405P"],
#                          "Endometrial": ["E542K", "E542Q", "E545K", "E545G", "G1007R", "Y1021H", "Y1021C", "A1035V", "M1043I", "H1047Y", "H1047R", "G1050D", "T1052K", "H1065L"],
#                          "Glioblastoma": ["R88E", "E542K", "E545A", "T1025N", "Y1021N", "R88Q", "P298T", "R310C", "T1031G", "V344G", "E453K", "E545K", "Y1021C", "M1043I", "N1044S", "H1047Y", "G1049S"]}
# point_mutants = cancer_mutations_dict[cancer]
# final_out = f"normalized_{cancer}_mutants_mhcii.xlsx"
# with pd.ExcelWriter(final_out, engine='openpyxl') as w:
#     for pm in point_mutants:
#         df_var = pd.read_excel(f"/Users/mvsamudrala/CancerVaccine/Peptides/epitopes/Automated_Collection/all_variables_mhcii_{cancer}.xlsx", sheet_name=pm)
#         # df_var = log_fifty(df_var)
#         df_var = z_score_norm_test(df_var, "immunogenicity", "antigenicity", "binding affinity (nM)", "allergenicity")
#         df_var = min_max(df_var)
#         df_var.to_excel(w, header=True, index=False, sheet_name=pm)


