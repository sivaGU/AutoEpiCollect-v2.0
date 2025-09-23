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
from selenium.webdriver.chrome.options import Options
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from functools import partial
import streamlit as st


st.set_page_config(page_title="Multipage App")
st.markdown("<h1 style='text-align: center;'>AutoEpiCollect v2.0</h1>", unsafe_allow_html=True)
st.image("personalized_vac_ga.png")
st.markdown("<p style='text-align: center;'>Welcome to AutoEpiCollect v2.0. This updated GUI interface will allow you to design your own personalized cancer vaccine. There are two different kinds vaccine designs that we offer:</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: left;'>1. <strong>Regular Cancer Vaccine Design</strong>: With this type of design, you will be able to determine top potentially immunogenic peptide epitopes for a cancer vaccine targeting multiple mutations in a singular target gene. Download the software for this cancer vaccine design from this citation:</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: left;'>1. <strong>Samudrala, M.; Dhaveji, S.; Savsani, K.; Dakshanamurthy, S. AutoEpiCollect, a Novel Machine Learning-Based GUI Software for Vaccine Design: Application to Pan-Cancer Vaccine Design Targeting PIK3CA Neoantigens. Bioengineering 2024, 11, 322. https://doi.org/10.3390/bioengineering11040322</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>For extensive documentation about AutoEpiCollect, please visit: https://autoepicollect.readthedocs.io/en/latest/</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: left;'>2. <strong>Personalized Cancer Vaccine Design</strong>: Our new release features an updated personalized cancer vaccine design that we have developed to target multiple genes using RNAseq data.</p>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>Please click on the '>' in the top left corner to access the sidebar and navigate to the personalized cancer vaccine design software.</p>", unsafe_allow_html=True)

st.sidebar.success("Select a page above")

