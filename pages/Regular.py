import sys
import subprocess
import pandas as pd
import numpy as np
import os
from pathlib import Path
import platform
import openpyxl
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
from time import sleep
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver import ChromeOptions
import chromedriver_autoinstaller
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from functools import partial
import streamlit as st


st.markdown("<h1 style='text-align: center;'>Regular Cancer Vaccine Design</h1>", unsafe_allow_html=True)
form1 = st.form(key="Options")
form1.markdown("<p style='text-align: center;'>Please input either the name for your target gene or a file for the input gene in .fasta format.</p>", unsafe_allow_html=True)
col1, col2 = form1.columns(2)
gene_name = col1.text_input("Gene Name")
gene_file = col2.file_uploader("Input Gene File", type=['fasta'])
mhc_classes = form1.multiselect("Choose your MHC class(es)", ["Class I", "Class II"])
cancer_mutations = form1.text_area("Input your cancer subtypes and corresponding point mutations", "Cancer: PM1,PM2,PM3,etc.\nExample input:\nGBM: R38H,G106V,E545K\nCRC: H1047R,Y1021K", height=150)
epitope_properties = form1.multiselect("Choose what epitope characteristics you want to collect", ["Immunogenicity", "Antigenicity", "Allergenicity", "Aliphatic Index", "GRAVY Score", "Isoelectric Point", "Half-Life", "Instability Index", "Toxicity", "IFN-γ Release"])
scoring_function = form1.checkbox("Use machine learning-driven epitope scoring function")
filtering = form1.checkbox("Filter epitopes using collected properties")
population_coverage = form1.checkbox("Calculate population coverage of final vaccine composition")

form1.form_submit_button("Submit Params")