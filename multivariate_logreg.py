import pandas as pd
import numpy as np
from numpy import mean
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import sys
import openpyxl


# Input format: python multivariate_logreg.py [I or II] input_file output_file
# cd4 = pd.read_excel("refactored_trainingset_cd4.csv")
# cd8 = pd.read_excel("refactored_trainingset_cd8.csv")
# # df_test["predicted potential"] = [None for i in range(df_test.shape[0])]
#
# # X = df[["norm immunogenicity", "norm antigenicity", "norm allergenicity", "norm binding affinity"]].values
# # # X = df[["norm immunogenicity", "norm antigenicity", "norm binding affinity", "norm allergenicity"]].values
# # # X = df[["norm immunogenicity", "norm antigenicity", "norm binding affinity"]].values
# # # X = df[["Immunogenicity", "Binding Affinity"]].values
# # y = df[["Potential"]]
#
# cd8["logBindingAffinity"] = np.log(cd8["Binding Affinity"])
# cd4["logBindingAffinity"] = np.log(cd4["Binding Affinity"])

mhc = sys.argv[1]

if mhc == "I":
    training_df = pd.read_csv("refactored_trainingset_cd8.csv")
    training_df["logBindingAffinity"] = np.log(training_df["Binding Affinity"])
    X = training_df[["norm immunogenicity", "norm antigenicity", "norm allergenicity", "logBindingAffinity"]].values
    y = training_df[["Result"]]
else:
    training_df = pd.read_csv("refactored_trainingset_cd4.csv")
    training_df["logBindingAffinity"] = np.log(training_df["Binding Affinity"])
    X = training_df[["norm immunogenicity", "norm antigenicity", "norm allergenicity", "logBindingAffinity"]].values
    y = training_df[["Result"]]
logr = LogisticRegression()
logr.fit(X, y)

infile = sys.argv[2]
outfile = sys.argv[3]
df_normalized_epitopes = pd.read_excel(infile)

for i in range(df_normalized_epitopes.shape[0]):
    df_normalized_epitopes["log binding affinity"] = np.log(df_normalized_epitopes["binding affinity (nM)"])
    df_x = df_normalized_epitopes[
        ["norm immunogenicity", "norm antigenicity", "norm allergenicity", "log binding affinity"]].values
    predicted_potential = logr.predict_proba(df_x)[:, 1]
    df_normalized_epitopes["potential"] = predicted_potential
df_normalized_epitopes_ranked = df_normalized_epitopes.sort_values(by=["potential"], ascending=False)
df_normalized_epitopes_ranked.to_excel(outfile, index=False, header=True)

# return df_normalized_epitopes_ranked

# cd8.to_excel("linreg_cd8_training.xlsx", header=True, index=False)
# cd4.to_excel("linreg_cd4_training.xlsx", header=True, index=False)

# cd8_train, cd8_test = train_test_split(cd8, test_size=0.3)
# cd8_x = cd8_train[["norm immunogenicity", "norm antigenicity", "norm allergenicity", "logBindingAffinity"]].values
# cd8_y = cd8_train[["Potential"]]
# cd4_train, cd4_test = train_test_split(cd4, test_size=0.3)
# cd4_x = cd4_train[["norm immunogenicity", "norm antigenicity", "norm allergenicity", "logBindingAffinity"]].values
# cd4_y = cd4_train[["Potential"]]


# kf = KFold(n_splits=5, random_state=1, shuffle=True)
# regr = linear_model.LinearRegression()
# regr.fit(cd8_x, cd8_y)
# print(mean(cross_val_score(regr, cd8_x, cd8_y, scoring='r2', cv=kf, n_jobs=-1)))
# # print("Immunogenicity Coefficient, Antigenicity Coefficient, Binding Affinity Coefficient")
# # print("Immunogenicity Coefficient, Binding Affinity Coefficient")
# # print(regr.coef_)
#
# kf = KFold(n_splits=5, random_state=1, shuffle=True)
# regr = linear_model.LinearRegression()
# regr.fit(cd4_x, cd4_y)
# print(mean(cross_val_score(regr, cd4_x, cd4_y, scoring='r2', cv=kf, n_jobs=-1)))

# df_test = pd.read_excel("tweaked_training_set2.xlsx")
# for i in range(df_test.shape[0]):
#     predicted_potential = regr.predict([[df_test["norm immunogenicity"][i], df_test["norm antigenicity"][i], df_test["norm binding affinity"][i]]])
#     df_test.at[i, "predicted potential"] = predicted_potential
#
# df_test.to_excel("predicted_output_test_set.xlsx", header=True, index=False)


cancer = "CRC"
# #
# cancer_mutations_dict = {"CRC": ["R38H", "R88Q", "G106V", "C420R", "E453Q", "E542K", "E545K", "R1023Q", "M1043I", "H1047R"],
#                          "Meningioma": ["E110K", "I391M", "R108H", "G914R", "N345K", "E453K", "Y165H", "H1047R", "E545K"],
#                          "BC": ["E542K", "E542V", "E545K", "Q546E", "Q546R", "H1047L", "H1047R", "N345K", "E726K", "C420R", "G118D", "E453K", "Q546K", "G1049R", "M1043I", "K111E", "E81K", "E545A", "E545G", "N1044K", "S405P"],
#                          "Endometrial": ["E542K", "E542Q", "E545K", "E545G", "G1007R", "Y1021H", "Y1021C", "A1035V", "M1043I", "H1047Y", "H1047R", "G1050D", "T1052K", "H1065L"],
#                          "Glioblastoma": ["R88E", "E542K", "E545A", "T1025N", "Y1021N", "R88Q", "P298T", "R310C", "T1031G", "V344G", "E453K", "E545K", "Y1021C", "M1043I", "N1044S", "H1047Y", "G1049S"]}
# cancer_mutations_dict = {"CRC": ["R38H", "R88Q", "G106V", "C420R", "E453Q", "E542K", "E545K", "R1023Q", "M1043I", "H1047R"]}
# point_mutants = cancer_mutations_dict[cancer]
# # final_out = f"predicted_output_{cancer}_mutants_training_mhcii.xlsx"
# final_out = f"predicted_output_{cancer}_mutants_ctraining.xlsx"
# with pd.ExcelWriter(final_out, engine='openpyxl') as w:
#     for pm in point_mutants:
#         # df_test = pd.read_excel(f"normalized_{cancer}_mutants_mhcii.xlsx", sheet_name=pm)
#         df_test = pd.read_excel(f"normalized_{cancer}_mutants_mhci.xlsx", sheet_name=pm)
#         for i in range(df_test.shape[0]):
#             predicted_potential = regr.predict([[df_test["norm immunogenicity"][i], df_test["norm antigenicity"][i], df_test["norm allergenicity"][i], df_test["norm binding affinity"][i]]])
#             # predicted_potential = regr.predict([[df_test["norm immunogenicity"][i], df_test["norm antigenicity"][i], df_test["norm binding affinity"][i]]])
#             df_test.at[i, "predicted potential"] = predicted_potential
#         df_test.to_excel(w, header=True, index=False, sheet_name=pm)

# Accidentally ran neoepitope training set on deepimmuno set
# 3 = iedb regular epitopes, 4 = cedar cancer neoepitopes, 5 = deepimmuno training set
# tweaked = iedb regular, tweaked2 = cedar, tweaked3 = deepimmuno

# cleaned_training_set_mhcii: n=600
# cleaned_training_set_mhcii2: n=685
# cleaned_training_set_mhcii3: n=4500
