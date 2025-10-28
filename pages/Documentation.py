import sys
import shutil
import subprocess
import pandas as pd
import numpy as np
import os
import re
import io
from pathlib import Path
import zipfile
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
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from webdriver_manager.core.os_manager import ChromeType
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from functools import partial
import urllib.parse, urllib.request
import streamlit as st
import json
import time
from typing import List, Dict, Any
import requests


st.markdown(
        """
        ## Personalized Cancer Vaccine Design — Documentation

        Welcome to the documentation tab! Below you’ll find a detailed, step-by-step guide to what each section of the 
        Personalized Cancer Vaccine workflow does and how to use it.

        ## ① Generate Aggregated MHC Files

        ### 1. Upload your .CSV output file from Partek Flow.

        Reads your Partek Flow CSV of variants and returns a list of genes with point mutations, sorted by most mutations.
        The .csv file can only be from Partek Flow. Partek analyzes desired tissue samples to generate a list of potential mutations.

        **How it works:**
        1. Reads the uploaded CSV.
        2. Filters rows where **Gene section** contains “Exon”.
        3. Keeps only rows where **AA change** is not null.
        4. Excludes any nonsense mutations (those containing `*`).
        5. Groups by **Gene ID**, collects all AA changes into a list.
        6. Sorts genes by number of mutations (descending).
        7. Returns the list.

        ### 2. Select your vaccine design module

        There are two different modules "pathways" that AutoEpiCollect 2.0 can use to select epitopes for the cancer vaccine.

        1. The first is by verifying the carcinogenic properties of the filtered genes from Partek Flow using the Human 
        Protein Atlas Database. This is considered a more generalized workflow that keeps genes if they have any 
        potential carcinogenic effect. 
        2. The second module is more cancer-specific. If you select module 2, you will be able to enter a specific cancer
        subtype you are interested in treating. AutoEpiCollect 2.0 will then enter this cancer subtype into the COSMIC 
        Database and select genes that overlap with the top 20 most commonly mutated genes for the desired type of cancer. 

        Once you submit your Partek Flow .csv file and it is successfully read, another box will appear depending on the module you have selected.
        - If you selected module 1, you will be able to choose the number of genes that you want AutoEpiCollect 2.0 to use during the epitope selection process.
        - If you selected module 2, you will be able to enter your cancer subtype of interest.

        ### 3. Choose MHC class(es) and additional epitope characteristics

        The final options in the first phase allow you to choose which MHC classes and additional epitope 
        characteristics you want AutoEpiCollect 2.0 to collect. 

        1. MHC Classes: You can either choose from MHC Class I, Class II, or both. "Class I" means that AutoEpiCollect 2.0 will 
        collect and analyze epitopes that bind MHC Class I molecules (these present peptides to CD8+ T cells. Likewise, "Class II"
        means that collected epitopes will bind MHC Class II molecules (these present peptides to CD4+ T cells). 
        2. Additional Epitope Characteristics: AutoEpiCollect 2.0 will automatically obtain three epitope characteristics, since these
        are used for ranking the epitopes via machine learning (logistic regression): binding affinity, immunogenicity, 
        allergenicity. Users have the options to select additional characteristics that AutoEpiCollect will obtain for
        each epitope.
           - Aliphatic Index: The relative volume occupied by aliphatic (non-aromatic) side chains—namely alanine, 
           valine, isoleucine and leucine—in a peptide or protein. A higher aliphatic index often correlates with 
           greater thermal and structural stability. Stable peptides are less likely to degrade before they’re 
           taken up by antigen-presenting cells, improving vaccine shelf life and in-vivo persistence. Tool used for collection: **ProtParam**
           - GRAVY Score: The arithmetic mean of hydropathy values for all amino acids in the sequence. Positive values 
           are more hydrophobic; negative are more hydrophilic.Hydrophobicity affects solubility and how peptides 
           interact with both MHC binding grooves and the aqueous environment. Moderately hydrophilic peptides tend 
           to be more soluble in vaccine formulations and better processed by antigen-presenting cells. Tool used for collection: **ProtParam**
           - Isoelectric Point: The pH at which the peptide carries no net electrical charge. Peptides whose pI is far 
           from physiological pH (7.4) may aggregate or precipitate in solution. Choosing epitopes with pI near 
           physiological pH improves solubility and bioavailability once injected. Tool used for collection: **ProtParam**
           - Half-Life: What it is: An in-silico prediction of how long the peptide remains intact in a given 
           environment before being degraded by proteases. Longer half-life means 
           the peptide will persist longer in circulation or inside antigen-presenting cells, increasing the window 
           during which it can be loaded onto MHC molecules and recognized by T cells. Tool used for collection: **ProtParam**
           - Instability Index: A score (below 40 indicates stable; above 40 indicates unstable) predicting whether the peptide will 
           remain folded and intact under in-vitro conditions. Unstable peptides are prone to rapid degradation or 
           unfolding, which can lead to poor MHC binding or premature clearance—both of which reduce vaccine efficacy. Tool used for collection: **ProtParam**
           - Toxicity: A prediction of whether the peptide sequence is likely to be toxic or harmful to human cells. You 
           want to avoid peptides that could damage host tissues or provoke unwanted side-effects. Non-toxic epitopes 
           ensure the vaccine is safe. Tool used for collection: **ToxinPred**
           - IFN-γ Release: A prediction of a peptide’s ability to induce Interferon-gamma secretion by T cells 
           (usually CD4⁺ Th1 or CD8⁺ CTLs). IFN-γ is a key cytokine for anti-tumor immunity—it activates macrophages 
           and enhances cytotoxic T cell responses. Peptides that drive strong IFN-γ release are more likely to elicit 
           effective, tumor-killing immune responses. Tool used for collection: **IFNEpitope**

        ### 4. Resultant epitopes and FASTA file for antigenicity

        After clicking the submission button, AutoEpiCollect 2.0 will collect epitope and epitope characteristics for each 
        gene and point mutation of interest (using the parameters chosen). When the initial run is complete, you will
        be able to download a ZIP file that contains the following files depending on which MHC Class option you selected.
        - 'aggregated_mhc_I.xlsx' + 'all_peptides_anti_mhci.fasta'
        - 'aggregated_mhc_II.xlsx` + 'all_peptides_anti_mhcii.fasta'

        The aggregated .xlsx files contain a master spreadsheet of all the epitopes collected with each epitope's 
        desired characteristics, as well as the point mutation and gene it originated from. The .fasta file contains all 
        unique epitopes in FASTA format. **This is important as it is file you will need to submit to the VaxiJen tool to 
        collect antigenicity values. Antigenicity values are required in order to use AutoEpiCollect 2.0's ranking function. 
        Since VaxiJen currently has a CAPTCHA system in place, AutoEpiCollect cannot 
        automate webscraping from this website. Follow the steps below to obtain the necessary files from VaxiJen for 
        submission in the next tab.**

        1. Visit https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html. This is the VaxiJen website where you will
        obtain antigenicity values
        2. Click "Choose File", then upload the 'all_peptides_anti_mhc(i or ii).fasta' file. You will need to do this 
        process once for each MHC Class you are running.
        3. Select 'Tumour' for the "Target Organism".
        4. Click "Submit" (If there is an error saying "Web server is returning an unknown error", please wait a few hours and try again later. It is likely an issue on VaxiJen's end.)
        5. Copy and paste the enter output into a new .txt file. **Please include the MHC Class of the peptides in the 
        name of the .txt file in the format "mhc_(I or II)". For example, "vaxijen_values_mhc_I.txt" or "vaxijen_values_mhc_II.txt".**
        Text below shows a sample output from VaxiJen copy and pasted into a new .txt file:
        ```txt
        Your Sequence:

        >peptide_1

        GLAGLLGLI

        Overall Prediction for the Protective Antigen = -1.2806 ( Probable NON-ANTIGEN ).

        Your Sequence:

        >peptide_2

        GLITCLICGV

        Overall Prediction for the Protective Antigen = 0.5635 ( Probable ANTIGEN ).

        Your Sequence:

        >peptide_3

        GLLGLITCL

        Overall Prediction for the Protective Antigen = -0.4888 ( Probable NON-ANTIGEN ).
        ```

        Now, gather you're aggregated .xlsx files with epitope data and the .txt files with copied antigenicity output 
        data and proceed to the next tab (Process Aggregated Files).

        ## ② Process Aggregated Files

        Once you have your aggregated `.xlsx` files and antigenicity `.txt` outputs from VaxiJen, use this second tab to
        merge, rank, filter, and—optionally—calculate population coverage.

        ### 1. Upload antigenicity `.txt` files  
        - Click "Choose Files" under "Upload new antigenicity values in TXT".  
        - You may supply one file per MHC class (e.g. 'vaxijen_values_mhc_I.txt', 'vaxijen_values_mhc_II.txt').  
        - Internally, each file is parsed by a function which uses regex to extract:
          1. `peptide_#` sequence  
          2. antigenicity score  
        - These are converted into DataFrames for merging.

        ### 2. Upload aggregated .xlsx epitope files  
        - Use "Upload aggregated all_variables_mhc*.xlsx" to select the 'aggregated_mhc_I.xlsx' and/or  
          'aggregated_mhc_II.xlsx' you downloaded in tab 1.  
        - AutoEpiCollect 2.0 reads each sheet into a DataFrame and concatenates by MHC class.

        ### 3. Merge antigenicity with your epitope data  
        - If you provided matching .txt and .xlsx for a class, the antigenicity scores are joined on 'peptide'.  
        - Any peptides lacking a score will show 'NaN'—you can still rank/filter, but we recommend including both files.

        ### 4. Select which extra epitope characteristic columns already exist  
        - The multiselect "Which additional characteristics exist" tells the app which flags are already in your  
          .xlsx files.  
        - By default your tables already include:
          - Binding Affinity  
          - Immunogenicity
          - Allergenicity
        - Any boxes you check here (e.g. Toxicity, Instability Index) will be used in scoring/filtering.

        ### 5. Choose processing options  
        - Use machine learning–driven epitope scoring  
          - Runs normalization + logistic-regression scoring for each peptide.  
        - Filter epitopes 
          - Applies the filtration function to select only top candidates based on **half-life, instability, toxicity, and IFN-γ release parameters.**
        - **Calculate population coverage**  
          - Computes both "regular" and "optimized" coverage sheets for epitopes.

        ### 6. Download your final ZIP  
        After you click the submit button and AutoEpiCollect 2.0 is finished running, you will be able to collect your final
        .xlsx files with all epitopes, characteristics, and scores. **Of note, the column labeled "potential" contains numbers from
        0 to 1 indicating how well an epitope will perform in a cancer vaccine based on our machine-learning ranking model. Higher
        numbers indicate better performance.**
        - Contains:
          - 'final_aggregated_mhc_I.xlsx' and/or 'final_aggregated_mhc_II.xlsx' 
          - (if you checked population coverage) 'PopCov_Regular' & 'PopCov_Optimized' sheets  
        - Use these to finalize vaccine candidate selection and reporting.

        **Tip:**  
        - Always verify that your antigenicity .txt filenames include the MHC class.  
        - If you add new epitope features in the aggregated '.xlsx', remember to check them in step 4 so they’re incorporated into scoring/filtering.
        """
    )
