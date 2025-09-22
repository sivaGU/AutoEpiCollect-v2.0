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
import streamlit as st


def make_driver():
    opts = Options()
    # required for headless in containers
    opts.add_argument("--headless=new")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    # stability niceties
    opts.add_argument("--disable-gpu")
    opts.add_argument("--window-size=1920,1080")

    # point explicitly to the system Chromium we install via packages.txt
    chrome_bin = shutil.which("chromium") or shutil.which("google-chrome")
    if chrome_bin:
        opts.binary_location = chrome_bin

    # get a matching chromedriver
    driver_path = ChromeDriverManager(chrome_type=ChromeType.CHROMIUM).install()
    return webdriver.Chrome(service=Service(driver_path), options=opts)


# Function to obtain the fasta-formatted gene sequence of interest from UniProt
def get_gene_sequence(target_gene):
    #  Instantiates the web scraping tool and accesses UniProt website
    driver = make_driver()
    driver.get('https://www.uniprot.org/')

    sleep(1)

    # Anytime By.XPATH is used, it means using the inspect element button to find the html pointer of a certain button
    # or text box or text
    search_box = driver.find_element(By.XPATH,
                                     '/html/body/div[1]/div/div/main/div/div[1]/div/section/form/div[2]/input')
    search_box.send_keys(target_gene)

    sleep(1.5)

    submit_button = driver.find_element(By.XPATH, '//*[@type="submit"]')
    submit_button.click()

    sleep(3.5)

    try:
        choose_table = driver.find_element(By.XPATH, '/html/body/form/div/span/label[2]/input')
    except NoSuchElementException:
        print("No format selector")
    else:
        choose_table.click()

    sleep(0.5)

    view_results_button = driver.find_element(By.XPATH, '/html/body/form/div/section/button')
    view_results_button.click()

    sleep(1)

    filter_human = driver.find_element(By.XPATH, '/html/body/div[1]/div/div/div/aside/div/ul/li[2]/div/ul/li[1]/a')
    filter_human.click()

    sleep(1)

    gene_button = driver.find_element(By.XPATH,
                                      '/html/body/div[1]/div/div/div/main/div[4]/table/tbody/tr[1]/td[2]/span/a')
    gene_button.click()

    sleep(5)

    download_button = driver.find_element(By.XPATH, '/html/body/div[1]/div/div/div/main/div/div[2]/div/button[1]')
    driver.execute_script("arguments[0].scrollIntoView();", download_button)
    driver.execute_script("arguments[0].click();", download_button)

    sleep(0.5)

    fasta_button = driver.find_element(By.XPATH, '/html/body/aside/div[2]/fieldset[2]/label/select/option[2]')
    fasta_button.click()

    download_button2 = driver.find_element(By.XPATH, '/html/body/aside/div[2]/section[1]/a')
    download_button2.click()

    sleep(2)

    child_tab = driver.window_handles[1]
    driver.switch_to.window(child_tab)

    result_fasta = driver.find_element(By.XPATH, '/html/body/pre').text
    with open(f"{target_gene}.fasta", "w") as g:
        g.write(result_fasta)

    for tab in driver.window_handles:
        driver.switch_to.window(tab)
        # driver.close()
        driver.quit()


# Function to create whole mutant fasta genes from each point mutation of interest
def make_mutant_genes(mutant_list, gene_seq, parent_dir):
    for m in mutant_list:
        file_out = parent_dir / "mutant_gene_fastas" / f"{m}.fasta"
        os.makedirs(os.path.dirname(file_out), exist_ok=True)
        with open(file_out, "w") as f:
            for seq_record in SeqIO.parse(open(gene_seq, mode='r'), 'fasta'):
                wild_type = seq_record.seq
                new = m[-1]
                loc = m[:-1]
                loc = loc[1:]
                loc = int(loc) - 1
                new_seq = wild_type[:loc] + new + wild_type[loc + 1:]
                mutant_seq = SeqRecord(seq=new_seq, id=seq_record.id, description=seq_record.description)
                r = SeqIO.write(mutant_seq, f, 'fasta')
                if r != 1:
                    print('Error while writing sequence:  ' + seq_record.id)


# Function that uses the IEDB Tools API to obtain potential epitopes from the mutant gene fasta files, as well as
# the binding affinity of the epitopes to MHC I and II molecules
def get_epitopes_ba(mutant_list, mhc, parent_dir, log_out, log_lin):
    if mhc == "I":
        epitopes_dict = {}
        for m in mutant_list:
            target_sequence = ""
            epitope_lengths = ""
            fasta_file = parent_dir / "mutant_gene_fastas" / f"{m}.fasta"
            for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                sequence = str(seq_record.seq)
                loc = m[:-1]
                loc = loc[1:]
                loc = int(loc) - 1
                # target_sequence = sequence[loc - 19:loc] + sequence[loc:loc + 20]
                target_sequence = sequence[max(0, loc - 19): min(len(sequence), loc + 20)]
                print(target_sequence)
                with open("MHCI_HLA_input.txt", "r") as h:
                    alleles = h.read() * 2
                    alleles = alleles[:-1]
                epitope_lengths = "9," * 27 + "10," * 26 + "10"

            mhc_i = subprocess.run(["curl", "--data",
                                    f'method=netmhcpan_ba&sequence_text={target_sequence}&allele={alleles}&length={epitope_lengths}',
                                    "http://tools-cluster-interface.iedb.org/tools_api/mhci/"], capture_output=True,
                                   text=True)
            output = mhc_i.stdout
            table = StringIO(output)
            print(table)
            df = pd.read_table(table, sep=r"\s+")
            df = df.rename(columns={"ic50": "binding affinity (nM)"})
            column_titles = ["allele", "seq_num", "start", "end", "length", "peptide", "core", "icore",
                             "percentile_rank", "binding affinity (nM)"]
            df = df.reindex(columns=column_titles)
            for i in range(df.shape[0]):
                affinity = float(df["binding affinity (nM)"][i])
                if affinity <= 50:
                    df.at[i, "binding score"] = "STRONG"
                elif 50 < affinity <= 500:
                    df.at[i, "binding score"] = "NORMAL"
                elif 500 < affinity <= 5000:
                    df.at[i, "binding score"] = "WEAK"
                else:
                    df.at[i, "binding score"] = "N/A"
            print(df)
            epitopes_dict[m] = df
            print(f"{m} done")
            log_lin.append(f"{m} done")
            log_out.text("\n".join(log_lin))
    else:
        epitopes_dict = {}
        for m in mutant_list:
            target_sequence = ""
            epitope_lengths = ""
            fasta_file = parent_dir / "mutant_gene_fastas" / f"{m}.fasta"
            for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                sequence = str(seq_record.seq)
                loc = m[:-1]
                loc = loc[1:]
                loc = int(loc) - 1
                # target_sequence = sequence[loc - 49:loc] + sequence[loc:loc + 50]
                target_sequence = sequence[max(0, loc - 49): min(len(sequence), loc + 50)]
                print(target_sequence)
                with open("MHCII_HLA_input.txt", "r") as h:
                    alleles = h.read()
                    alleles = alleles[:-1]
                epitope_lengths = "15," * 26 + "15"

            mhc_ii = subprocess.run(["curl", "--data",
                                     f'method=netmhciipan&sequence_text={target_sequence}&allele={alleles}&length={epitope_lengths}',
                                     "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"],
                                    capture_output=True, text=True)
            output = mhc_ii.stdout
            table = StringIO(output)
            df = pd.read_table(table, sep=r"\s+")
            df = df.rename(columns={"ic50": "binding affinity (nM)"})
            column_titles = ["allele", "seq_num", "start", "end", "length", "core_peptide", "peptide", "rank",
                             "binding affinity (nM)"]
            df = df.reindex(columns=column_titles)
            for i in range(df.shape[0]):
                affinity = float(df["binding affinity (nM)"][i])
                if affinity <= 50:
                    df.at[i, "binding score"] = "STRONG"
                elif 50 < affinity <= 500:
                    df.at[i, "binding score"] = "NORMAL"
                elif 500 < affinity <= 5000:
                    df.at[i, "binding score"] = "WEAK"
                else:
                    df.at[i, "binding score"] = "N/A"
            print(df)
            epitopes_dict[m] = df
            print(f"{m} done")
            log_lin.append(f"{m} done")
            log_out.text("\n".join(log_lin))
    return epitopes_dict


# Function that only keeps the mutant epitopes from the list of all potential epitopes
def get_mutant_epitopes(mutant_list, mhc, all_epitopes_dict, parent_dir):
    if mhc == "I":
        epitopes_dict = {}
        for m in mutant_list:
            fasta_file = parent_dir / "mutant_gene_fastas" / f"{m}.fasta"
            df = all_epitopes_dict[m].copy()
            for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                sequence = seq_record.seq
                loc = m[:-1]
                loc = loc[1:]
                loc = int(loc) - 1
                start_index = max(0, loc - 9)
                end_index = min(len(sequence), loc + 10)
                target_epitope = sequence[start_index:loc] + sequence[loc:end_index]
                print(m + ": " + target_epitope)
                index = 0
                bad_epitopes_indexes = []
                loc += 1
                for pep in df["peptide"]:
                    # print(pep)
                    current_start_index = sequence.find(pep, loc - 10) + 1
                    # print(start_index)
                    current_end_index = current_start_index + len(pep) - 1
                    # print(end_index)
                    if pep not in target_epitope or not current_start_index <= loc <= current_end_index:
                        # if pep not in target_epitope:
                        bad_epitopes_indexes.append(index)
                    else:
                        print(f"{pep}, {current_start_index}, {current_end_index}")
                    index += 1
                # print(bad_epitopes_indexes)
                df_dropped = df.drop(bad_epitopes_indexes).reset_index(drop=True)
                epitopes_dict[m] = df_dropped.reset_index(drop=True)
    else:
        epitopes_dict = {}
        for m in mutant_list:
            fasta_file = parent_dir / "mutant_gene_fastas" / f"{m}.fasta"
            df = all_epitopes_dict[m].copy()
            # print(f"hello: {df['peptide']}")
            for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                sequence = seq_record.seq
                loc = m[:-1]
                loc = loc[1:]
                loc = int(loc) - 1
                # print(loc)
                start_index = max(0, loc - 14)
                end_index = min(len(sequence), loc + 15)
                # print(start_index)
                # print(end_index)
                target_epitope = sequence[start_index:loc] + sequence[loc:end_index]
                print(m + ": " + target_epitope)
                bad_epitopes_indexes = []
                loc += 1
                for i in range(df.shape[0]):
                    # start_index = int(df["start"][i])
                    # end_index = int(df["end"][i])
                    pep = df["peptide"][i]
                    # print(pep)
                    # print(loc)
                    current_start_index = sequence.find(pep, loc - 15) + 1
                    # print(current_start_index)
                    current_end_index = current_start_index + len(pep) - 1
                    # print(current_end_index)
                    if pep not in target_epitope or not current_start_index <= loc <= current_end_index:
                        # if pep not in target_epitope:
                        bad_epitopes_indexes.append(i)
                    else:
                        print(f"{pep}, {current_start_index}, {current_end_index}")
                # print(bad_epitopes_indexes)
                df_dropped = df.drop(bad_epitopes_indexes).reset_index(drop=True)
                epitopes_dict[m] = df_dropped
    return epitopes_dict


# Function that makes .txt and .fasta files of peptides created from each point mutation of interest, test
def get_peptides(point_mutants, mut_epitopes_dict, mhc, parent_dir):
    for m in point_mutants:
        file_out = f"{parent_dir}/Sequences/{m}peptides_{mhc}.txt"
        os.makedirs(os.path.dirname(file_out), exist_ok=True)
        df = mut_epitopes_dict[m]
        peptide_set = set()
        for pep in df.loc[:, "peptide"]:
            peptide_set.add(pep)
        with open(file_out, "w") as f_out:
            for pep in peptide_set:
                f_out.write(pep + "\n")
        file_out = parent_dir / "Sequences" / f"{m}peptides_{mhc}.fasta"
        os.makedirs(os.path.dirname(file_out), exist_ok=True)
        with open(file_out, "w") as f_out:
            count = 0
            for pep in peptide_set:
                heading = f">Epitope{count}"
                f_out.write(heading + "\n" + pep + "\n")
                count += 1


# Function that obtains the immunogenicity scores of MHC I peptides using an already downloaded library from the IEDB
def get_local_immunogenicity_mhci(immunogenicity_file, peptide_file, current_df):
    completed_run = subprocess.run(["python", immunogenicity_file, peptide_file], capture_output=True, text=True)
    output = completed_run.stdout
    output = output[output.find("peptide"):]
    table = StringIO(output)
    df = pd.read_table(table, sep=",")

    for i in range(df.shape[0]):
        result = df["score"][i]
        e = df["peptide"][i]
        current_df.loc[current_df["peptide"] == e, "immunogenicity"] = float(result)
    return current_df


# Function that obtain the immunogenicity scores of MHC II peptides by webscraping the IEDB website
# This function is a work in progress as this website is finnicky and slow right now
def get_immunogenicity_mhcii(peptide_list, p, current_df):
    driver = make_driver()
    driver.get('http://tools.iedb.org/CD4episcore/')

    elem = WebDriverWait(driver, 60).until(
        ec.presence_of_element_located((By.XPATH, '/html/body/div[3]/form/table/tbody/tr[3]/td[2]/textarea')))

    searchbox = driver.find_element(By.XPATH, '/html/body/div[3]/form/table/tbody/tr[3]/td[2]/textarea')
    print(p)
    searchbox.send_keys(p)

    sleep(1)

    threshold_button = driver.find_element(By.XPATH, '/html/body/div[3]/form/table/tbody/tr[9]/td[2]/select/option[10]')
    threshold_button.click()

    sleep(2)

    submit_button = driver.find_elements(By.XPATH, '/html/body/div[3]/form/table/tbody/tr[12]/th/div/input[1]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 120).until(
        ec.presence_of_element_located((By.XPATH, '/html/body/div[3]/div[1]/h2')))

    for x in range(len(peptide_list)):
        result = driver.find_element(By.XPATH, f'/html/body/div[3]/div[1]/div[3]/table/tbody/tr[{x + 1}]/td[6]').text
        e = driver.find_element(By.XPATH, f'/html/body/div[3]/div[1]/div[3]/table/tbody/tr[{x + 1}]/td[3]').text
        current_df.loc[current_df["peptide"] == e, "immunogenicity"] = float(result)

    # driver.close()
    driver.quit()
    return current_df


# Function that obtains the antigenicity scores of both MHC I and II peptides.
# This function is a work in progress as VaxiJen v2.0 has a human verification loop that is blocking automation
def get_antigenicity(peptide_list, peptide_fasta, current_df):
    driver = make_driver()
    driver.get('http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html')

    file_box = driver.find_element(By.XPATH, '//input[@type="FILE"]')
    file_box.send_keys(str(peptide_fasta))

    organism_button = driver.find_element(By.XPATH,
                                          '/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[2]/p/select/option[3]')
    organism_button.click()

    submit_button = driver.find_elements(By.XPATH,
                                         '/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[3]/td[2]/input[1]')
    submit_button[0].click()

    sleep(0.5)

    for x in range(len(peptide_list)):
        result = driver.find_element(By.XPATH,
                                     f'/html/body/div/table/tbody/tr[4]/td[3]/table/tbody/tr/td/b[{3 * (x + 1)}]').text
        e = driver.find_element(By.XPATH,
                                f'/html/body/div/table/tbody/tr[4]/td[3]/table/tbody/tr/td/font[{2 * (x + 1)}]').text
        current_df.loc[current_df["peptide"] == e, "antigenicity"] = float(result)

    # driver.close()
    driver.quit()
    return current_df


# Function that obtains the allergenicity of MHC class I epitopes using web scraping
def get_allergenicity_algpred(peptide_list, pf, current_df):
    driver = make_driver()
    driver.get('https://webs.iiitd.edu.in/raghava/algpred2/batch.html')

    text_box = driver.find_element(By.XPATH,
                                   '/html/body/header/div[3]/section/form/table/tbody/tr/td/font/p/font[1]/textarea')
    text_box.send_keys(pf)

    threshold_value = driver.find_element(By.XPATH,
                                          '/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[1]/b/select/option[10]')
    threshold_value.click()

    sleep(0.5)

    submit_button = driver.find_elements(By.XPATH,
                                         '/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[2]/input[2]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 300).until(
        ec.presence_of_element_located((By.XPATH, '/html/body/header/div[3]/main/h1/strong/font/b')))

    for x in range(len(peptide_list)):
        result = driver.find_element(By.XPATH,
                                     f'/html/body/header/div[3]/main/div/table[2]/tbody/tr[{x + 1}]/td[5]').text
        current_df.loc[current_df["peptide"] == peptide_list[x], "allergenicity"] = float(result)

    # driver.close()
    driver.quit()
    return current_df


# Function that obtains the allergenicity of MHC Class II epitopes using webscraping
def get_allergenicity_netallergen(peptide_list, pf, current_df, peptide_fasta):
    input_fasta = pf + ">Epitope119" + "\n" + "TIETLMLLALIAAAA"

    driver = make_driver()
    driver.get('https://services.healthtech.dtu.dk/services/NetAllergen-1.0/')

    try:
        accept_cookies = driver.find_element(By.XPATH, '/html/body/div[5]/div[4]/div[2]')
    except NoSuchElementException:
        print("No cookies")
    else:
        accept_cookies.click()

    sleep(3)

    # text_box = driver.find_element(By.XPATH, '/html/body/main/div/div[3]/div/div[2]/div[1]/form/table/tbody/tr[1]/td/textarea')
    # text_box = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[1]/form/table/tbody/tr[1]/td/textarea')
    # driver.execute_script("arguments[0].scrollIntoView();", text_box)
    # text_box.send_keys(input_fasta)

    file_box = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[1]/form/table/tbody/tr[1]/td/input')
    driver.execute_script("arguments[0].scrollIntoView();", file_box)
    file_box.send_keys(str(peptide_fasta))

    sleep(1)

    # blast_button = driver.find_element(By.XPATH, '/html/body/main/div/div[3]/div/div[2]/div[1]/form/table/tbody/tr[4]/td/input')
    # blast_button = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[1]/form/table/tbody/tr[4]/td/input')
    # driver.execute_script("arguments[0].click();", blast_button)

    # submit_button = driver.find_element(By.XPATH, '/html/body/main/div/div[3]/div/div[2]/div[1]/form/table/tbody/tr[5]/td/input[1]')
    submit_button = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[1]/form/table/tbody/tr[5]/td/input[1]')
    driver.execute_script("arguments[0].click();", submit_button)

    elem = WebDriverWait(driver, 10000).until(
        ec.presence_of_element_located((By.XPATH, '/html/body/font/table/tbody/tr/td[3]/h2')))

    output = driver.find_element(By.XPATH, '/html/body/pre').text
    output = output[:output.find("Explain")]
    table = StringIO(output)
    df = pd.read_table(table, sep=r"\s+")

    count = 0
    for pep in peptide_list:
        result = df["Score_60F"][count]
        current_df.loc[current_df["peptide"] == pep, "allergenicity"] = float(result)
        count += 1
    driver.quit()
    return current_df


# Function that access ProtParam and obtains the half-life, instability, aliphatic index, isoelectric point, and
# GRAVY score of peptides. Only half-life and instability are used for filtering epitopes.
def get_protparam(peptide_list, h, ins, ali, iso, g, current_df):
    driver = make_driver()
    for e in peptide_list:
        driver.get('https://web.expasy.org/protparam/')

        # searchbox = driver.find_element(By.XPATH, '//*[@id="sib_body"]/form/textarea')
        searchbox = driver.find_element(By.XPATH, '/html/body/main/div/form/textarea')
        searchbox.send_keys(e)

        # submitButton = driver.find_elements(By.XPATH, '/html/body/div[2]/div[2]/form/p[1]/input[2]')
        submitButton = driver.find_element(By.XPATH, '/html/body/main/div/form/input[3]')
        submitButton.click()

        sleep(2)

        # results = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/pre[2]').text
        results = driver.find_element(By.XPATH, '/html/body/main/div/pre[2]').text
        pi = results[
             results.find('Theoretical pI: ') + 16:results.find('\n', results.find('Theoretical pI: ') + 16, -1)]
        half_life = results[results.find('The estimated half-life is: ') + 28:results.find('hours', results.find(
            'The estimated half-life is: ') + 28, -1) - 1]
        instability = results[results.find('The instability index (II) is computed to be ') + 45:results.find('\n',
                                                                                                              results.find(
                                                                                                                  'The instability index (II) is computed to be ') + 45,
                                                                                                              -1)]
        alipathy = results[
                   results.find('Aliphatic index: ') + 17:results.find('\n', results.find('Aliphatic index: ') + 17,
                                                                       -1)]
        gravy = results[results.find('Grand average of hydropathicity (GRAVY):') + 40:]

        if h:
            if ">" not in half_life:
                current_df.loc[current_df["peptide"] == e, "half-life"] = float(half_life)
            else:
                current_df.loc[current_df["peptide"] == e, "half-life"] = float(half_life[half_life.find(">") + 1:])
        if ins:
            current_df.loc[current_df["peptide"] == e, "instability"] = float(instability)
        if iso:
            current_df.loc[current_df["peptide"] == e, "isoelectric point"] = float(pi)
        if ali:
            current_df.loc[current_df["peptide"] == e, "aliphatic index"] = float(alipathy)
        if g:
            current_df.loc[current_df["peptide"] == e, "GRAVY score"] = float(gravy)
    driver.quit()
    return current_df


# Obtains the toxicity of peptides using webscraping
def get_toxicity(peptide_list, pf, current_df):
    driver = make_driver()
    driver.get('https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php')

    submission_box = driver.find_element(By.XPATH, '//*[@id="input_box"]')
    submission_box.send_keys(pf)

    submit_button = driver.find_elements(By.XPATH,
                                         '/html/body/table[2]/tbody/tr/td/form/fieldset/table[2]/tbody/tr[3]/td/input[2]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 60).until(
        ec.presence_of_element_located((By.XPATH, '/html/body/div[2]/table/thead/tr[1]/td[1]')))

    for x in range(len(peptide_list)):
        result = driver.find_element(By.XPATH, f'/html/body/div[2]/table/tbody/tr[{x + 1}]/td[4]').text
        e = driver.find_element(By.XPATH, f'/html/body/div[2]/table/tbody/tr[{x + 1}]/td[2]/a').text
        current_df.loc[current_df["peptide"] == e, "toxicity"] = result

    # driver.close()
    driver.quit()
    return current_df


# Obtains the IFN-gamma release of epitopes. IFN-gamma only used for filtering MHC Class II epitopes.
def get_ifn(peptide_list, pf, current_df):
    driver = make_driver()
    driver.get('https://webs.iiitd.edu.in/raghava/ifnepitope/predict.php')

    submission_box = driver.find_element(By.XPATH, '//*[@name="sequence"]')
    submission_box.send_keys(pf)

    submit_button = driver.find_elements(By.XPATH, '//*[@type="submit"]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 300).until(ec.presence_of_element_located((By.XPATH, '//*[@class="sorting_asc"]')))

    entry_number = driver.find_element(By.XPATH,
                                       '/html/body/div/div/div[3]/div/div/div/div/div[1]/label/select/option[4]')
    entry_number.click()

    for x in range(len(peptide_list)):
        result = driver.find_element(By.XPATH,
                                     f'/html/body/div/div/div[3]/div/div/div/div/table/tbody/tr[{x + 1}]/td[5]').text
        e = driver.find_element(By.XPATH,
                                f'/html/body/div/div/div[3]/div/div/div/div/table/tbody/tr[{x + 1}]/td[3]/a').text
        current_df.loc[current_df["peptide"] == e, "IFNg"] = result

    driver.quit()
    return current_df


# Applies a z-score and min-max normalization on immunogenicity, antigenicity, and allergenicity data.
# Binding affinity data has a log transformation applied to it
def normalize_data(df_collected_epitopes, mhc):
    immunogenicity_name = "immunogenicity"
    antigenicity_name = "antigenicity"
    allergenicity_name = "allergenicity"
    immunogenicity_mean = df_collected_epitopes[immunogenicity_name].mean()
    immunogenicity_std = df_collected_epitopes[immunogenicity_name].std()
    antigenicity_mean = df_collected_epitopes[antigenicity_name].mean()
    antigenicity_std = df_collected_epitopes[antigenicity_name].std()
    allergenicity_mean = df_collected_epitopes[allergenicity_name].mean()
    allergenicity_std = df_collected_epitopes[allergenicity_name].std()
    for i in range(df_collected_epitopes.shape[0]):
        immunogenicity_value = df_collected_epitopes[immunogenicity_name][i]
        antigenicity_value = df_collected_epitopes[antigenicity_name][i]
        allergenicity_value = df_collected_epitopes[allergenicity_name][i]
        new_immunogenicity_value = (immunogenicity_value - immunogenicity_mean) / immunogenicity_std
        new_antigenicity_value = (antigenicity_value - antigenicity_mean) / antigenicity_std
        new_allergenicity_value = (allergenicity_value - allergenicity_mean) / allergenicity_std
        df_collected_epitopes.at[i, "norm immunogenicity"] = new_immunogenicity_value
        df_collected_epitopes.at[i, "norm antigenicity"] = new_antigenicity_value
        df_collected_epitopes.at[i, "norm allergenicity"] = new_allergenicity_value

    immunogenicity_name = "norm immunogenicity"
    antigenicity_name = "norm antigenicity"
    allergenicity_name = "norm allergenicity"
    immunogenicity_max = df_collected_epitopes[immunogenicity_name].max()
    immunogenicity_min = df_collected_epitopes[immunogenicity_name].min()
    antigenicity_max = df_collected_epitopes[antigenicity_name].max()
    antigenicity_min = df_collected_epitopes[antigenicity_name].min()
    allergenicity_max = df_collected_epitopes[allergenicity_name].max()
    allergenicity_min = df_collected_epitopes[allergenicity_name].min()
    for i in range(df_collected_epitopes.shape[0]):
        immunogenicity_value = df_collected_epitopes[immunogenicity_name][i]
        antigenicity_value = df_collected_epitopes[antigenicity_name][i]
        allergenicity_value = df_collected_epitopes[allergenicity_name][i]
        if mhc == "I":
            new_immunogenicity_value = (immunogenicity_value - immunogenicity_min) / (
                    immunogenicity_max - immunogenicity_min)
        else:
            new_immunogenicity_value = (immunogenicity_value - immunogenicity_max) / (
                    immunogenicity_min - immunogenicity_max)
        new_antigenicity_value = (antigenicity_value - antigenicity_min) / (antigenicity_max - antigenicity_min)
        new_allergenicity_value = (allergenicity_value - allergenicity_max) / (allergenicity_min - allergenicity_max)
        df_collected_epitopes.at[i, immunogenicity_name] = new_immunogenicity_value
        df_collected_epitopes.at[i, antigenicity_name] = new_antigenicity_value
        df_collected_epitopes.at[i, allergenicity_name] = new_allergenicity_value
    return df_collected_epitopes


# Trains the logistic regression model and uses it for the scoring function. Scores the epitopes and ranks them.
def apply_scoring_function(df_normalized_epitopes, mhc):
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

    for i in range(df_normalized_epitopes.shape[0]):
        df_normalized_epitopes["log binding affinity"] = np.log(df_normalized_epitopes["binding affinity (nM)"])
        df_x = df_normalized_epitopes[
            ["norm immunogenicity", "norm antigenicity", "norm allergenicity", "log binding affinity"]].values
        predicted_potential = logr.predict_proba(df_x)[:, 1]
        df_normalized_epitopes["potential"] = predicted_potential
    df_normalized_epitopes_ranked = df_normalized_epitopes.sort_values(by=["potential"], ascending=False)
    return df_normalized_epitopes_ranked


# Takes the top 20 epitopes after ranking as been completed and applies the manual filtration
# using the exclusion criteria.
def get_filtered_epitopes(df_ranked_epitopes, mhc, h, ins, t, ifn):
    # if scoring:
    #     top_20_df = df_ranked_epitopes.head(20)
    # else:
    top_20_df = df_ranked_epitopes
    if mhc == "I":
        if h and ins and t:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40) & (top_20_df["toxicity"] == "Non-Toxin")]
        elif h and ins and not t:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40)]
        elif h and not ins and t:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["half-life"] > 1) & (top_20_df["toxicity"] == "Non-Toxin")]
        elif h and not ins and not t:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["half-life"] > 1)]
        elif not h and ins and t:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["instability"] < 40) & (top_20_df["toxicity"] == "Non-Toxin")]
        elif not h and ins and not t:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["instability"] < 40)]
        elif not h and not ins and t:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["toxicity"] == "Non-Toxin")]
        else:
            df_filtered_epitopes = top_20_df.copy()
    else:
        if h and ins and t and ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40) & (
                    top_20_df["toxicity"] == "Non-Toxin") & (top_20_df["IFNg"] == "POSITIVE")]
        elif h and ins and t and not ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40) & (
                    top_20_df["toxicity"] == "Non-Toxin")]
        elif h and ins and not t and ifn:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40) & (top_20_df["IFNg"] == "POSITIVE")]
        elif h and ins and not t and not ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40)]
        elif h and not ins and t and ifn:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["half-life"] > 1) & (top_20_df["toxicity"] == "Non-Toxin") & (
                        top_20_df["IFNg"] == "POSITIVE")]
        elif h and not ins and t and not ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["toxicity"] == "Non-Toxin")]
        elif h and not ins and not t and ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["IFNg"] == "POSITIVE")]
        elif h and not ins and not t and not ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1)]
        elif not h and ins and t and ifn:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["instability"] < 40) & (top_20_df["toxicity"] == "Non-Toxin") & (
                        top_20_df["IFNg"] == "POSITIVE")]
        elif not h and ins and t and not ifn:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["instability"] < 40) & (top_20_df["toxicity"] == "Non-Toxin")]
        elif not h and ins and not t and ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["instability"] < 40) & (top_20_df["IFNg"] == "POSITIVE")]
        elif not h and ins and not t and not ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["instability"] < 40)]
        elif not h and not ins and t and ifn:
            df_filtered_epitopes = top_20_df.loc[
                (top_20_df["toxicity"] == "Non-Toxin") & (top_20_df["IFNg"] == "POSITIVE")]
        elif not h and not ins and t and not ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["toxicity"] == "Non-Toxin")]
        elif not h and not ins and not t and ifn:
            df_filtered_epitopes = top_20_df.loc[(top_20_df["IFNg"] == "POSITIVE")]
        else:
            df_filtered_epitopes = top_20_df.copy()
    return df_filtered_epitopes.reset_index(drop=True)


# Uses PCOptim and PCOptim-CD to obtain the optimized list of filtered epitopes for population coverage analysis.
def get_optimized_epitopes(df_filtered_epitopes, mhc):
    filtered_epitopes_str = ""
    op_df = pd.DataFrame(columns=df_filtered_epitopes.columns)
    if mhc == "I":
        java_file = "PopCoverageOptimization.java"
        for i in range(df_filtered_epitopes.shape[0]):
            allele = df_filtered_epitopes["allele"][i]
            peptide = df_filtered_epitopes["peptide"][i]
            filtered_epitopes_str = f"{filtered_epitopes_str}{allele}\t{peptide}\n"
        filtered_epitopes_str = filtered_epitopes_str.rstrip("\n")
        mhc_i = subprocess.run(["java", java_file, filtered_epitopes_str], capture_output=True, text=True)
        output = mhc_i.stdout.rstrip("\n")
        print(output)
        if output != "":
            peptide_epitopes_list = output.split("\n")
            for combination in peptide_epitopes_list:
                allele_list = combination[combination.find("HLA"):].split(",")
                pep = combination[:combination.find("HLA") - 1]
                for allele in allele_list:
                    index = df_filtered_epitopes.index[
                        (df_filtered_epitopes["allele"] == allele) & (df_filtered_epitopes["peptide"] == pep)][0]
                    op_df.loc[len(op_df)] = df_filtered_epitopes.iloc[index]
    else:
        java_file = "CD4PopCoverageOptimization.java"
        allowed_alleles = ["HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
        for i in range(df_filtered_epitopes.shape[0]):
            allele = df_filtered_epitopes["allele"][i]
            allele_check = allele[:allele.find("*")]
            if allele_check not in allowed_alleles or "/" in allele:
                pass
            else:
                peptide = df_filtered_epitopes["peptide"][i]
                filtered_epitopes_str = f"{filtered_epitopes_str}{allele}\t{peptide}\n"
        filtered_epitopes_str = filtered_epitopes_str.rstrip("\n")
        mhc_ii = subprocess.run(["java", java_file, filtered_epitopes_str], capture_output=True, text=True)
        output = mhc_ii.stdout.rstrip("\n")
        print(output)
        if output != "":
            peptide_epitopes_list = output.split("\n")
            for combination in peptide_epitopes_list:
                allele_list = combination[combination.find("HLA"):].split(",")
                pep = combination[:combination.find("HLA") - 1]
                for allele in allele_list:
                    index = df_filtered_epitopes.index[
                        (df_filtered_epitopes["allele"] == allele) & (df_filtered_epitopes["peptide"] == pep)][0]
                    op_df.loc[len(op_df)] = df_filtered_epitopes.iloc[index]
    return op_df


# Helper function to convert the peptide allele list for input into PCOptim and PCOptim-CD
def convert_to_text(df):
    peptide_allele_dict = {}
    for i in range(df.shape[0]):
        peptide = df["peptide"][i]
        allele = df["allele"][i]
        if peptide not in peptide_allele_dict.keys():
            peptide_allele_dict[peptide] = ""
            peptide_allele_dict[peptide] = f"{peptide_allele_dict[peptide]}{allele},"
        else:
            peptide_allele_dict[peptide] = f"{peptide_allele_dict[peptide]}{allele},"
    file_lines = []
    for peptide in peptide_allele_dict.keys():
        alleles = peptide_allele_dict[peptide]
        alleles = alleles[:-1]
        file_lines.append(f"{peptide}\t{alleles}")
    return file_lines


# Helper function to run the population coverage analysis tool and generate population coverage graphs
def population_coverage_helper(pop_filename, reg, mhc, pep_filename, plot_folder, path):
    completed_run = subprocess.run(
        ["python", pop_filename, "-p", reg, "-c", mhc, "-f", pep_filename, "--plot", plot_folder], capture_output=True,
        text=True)
    output = completed_run.stdout
    print(output)
    output = output[output.find("World"):]
    output_list = output.split("\n")
    output_list = output_list[:19]
    if path == "pop_cov_path":
        output_list.insert(0, "Region\tCoverage\tAverage Hit\tpc90")
    else:
        output_list.insert(0, "Region\tOptimized Coverage\tAverage Hit\tpc90")
    output = "\n".join(output_list)
    table = StringIO(output)
    df = pd.read_table(table, sep="\t")
    return df


# Function to create and organize files for population coverage analysis. Uses function above to run the
# population covarage analysis tool that is downloaded.
def get_population_coverage(df_filtered_epitopes, df_optimized_epitopes, mhc, parent_dir):
    if mhc == "I":
        peptide_allele_filename = f"filtered_epitopes_mhci.txt"
        optimized_peptide_allele_filename = f"optimized_epitopes_mhci.txt"
        filtered_epitopes_txt = convert_to_text(df_filtered_epitopes)
        with open(peptide_allele_filename, "w") as fo:
            for line in filtered_epitopes_txt:
                fo.write(f"{line}\n")
        optimized_epitopes_txt = convert_to_text(df_optimized_epitopes)
        with open(optimized_peptide_allele_filename, "w") as fo:
            for line in optimized_epitopes_txt:
                fo.write(f"{line}\n")
    else:
        peptide_allele_filename = f"filtered_epitopes_mhcii.txt"
        optimized_peptide_allele_filename = f"optimized_epitopes_mhcii.txt"
        filtered_epitopes_txt = convert_to_text(df_filtered_epitopes)
        with open(peptide_allele_filename, "w") as fo:
            for line in filtered_epitopes_txt:
                fo.write(f"{line}\n")
        optimized_epitopes_txt = convert_to_text(df_optimized_epitopes)
        with open(optimized_peptide_allele_filename, "w") as fo:
            for line in optimized_epitopes_txt:
                fo.write(f"{line}\n")
    new_dir = Path(f"{parent_dir}/Population_Coverage_Plots/")
    os.makedirs(new_dir, exist_ok=True)
    pop_coverage_filename = parent_dir / "population_coverage" / "calculate_population_coverage.py"
    plot_output_folder = new_dir / f"Regular_Plots"
    os.makedirs(plot_output_folder, exist_ok=True)
    optimized_plot_output_folder = new_dir / f"Optimized_Plots"
    os.makedirs(optimized_plot_output_folder, exist_ok=True)
    regions = ["World", "East Asia", "Northeast Asia", "South Asia", "Southeast Asia", "Southwest Asia", "Europe",
               "East Africa", "West Africa", "Central Africa", "North Africa", "South Africa", "West Indies",
               "North America", "Central America", "South America", "Oceania"]
    region_string = ",".join(regions)
    regular_results = population_coverage_helper(pop_coverage_filename, region_string, mhc, peptide_allele_filename,
                                                 plot_output_folder, "pop_cov_path")
    optimized_results = population_coverage_helper(pop_coverage_filename, region_string, mhc,
                                                   optimized_peptide_allele_filename, optimized_plot_output_folder,
                                                   "op_pop_cov_path")
    return regular_results, optimized_results


def variant_filter(filename):
    variants = pd.read_csv(filename)
    exons = variants.loc[variants["Gene section"].str.contains("Exon")]
    exons_pm = exons.loc[exons["AA change"].notna()]
    exons_spm = exons_pm.loc[~exons_pm["AA change"].str.contains(r'\*', na=False)]
    gene_mutations_series = exons_spm.groupby('Gene ID')["AA change"].apply(list)
    sorted_gene_mutations_series = gene_mutations_series.sort_values(key=lambda x: x.map(len), ascending=False)
    return sorted_gene_mutations_series


def gene_filter_helper(gene):
    driver = make_driver()
    driver.get('https://www.proteinatlas.org/')

    sleep(1)

    search_box = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[2]/form/div/input')
    search_box.send_keys(gene)

    submit_button = driver.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[2]/form/div/button')
    submit_button.click()

    sleep(0.5)

    select_gene = driver.find_element(By.XPATH, '/html/body/table/tbody/tr/td[2]/div/table[2]/tbody[1]/tr/td[2]/a')
    select_gene.click()

    sleep(2)

    specificity_box = driver.find_element(By.XPATH, '/html/body/table/tbody/tr/td[2]/div/table[5]/tbody/tr[3]/td').text
    driver.quit()
    if "carcinoma" in specificity_box.lower():
        return True
    else:
        return False


def gene_checker_helper(cancer_name):
    driver = make_driver()
    driver.get('https://cancer.sanger.ac.uk/cosmic')

    login_email = "tbp8up@virginia.edu"
    login_password = "R@j@_181019"

    sleep(1)

    login_email_box = driver.find_element(By.XPATH, '/html/body/div[2]/div/div/form/dl/dd[1]/input')
    login_email_box.send_keys(login_email)

    sleep(0.5)

    login_password_box = driver.find_element(By.XPATH, '/html/body/div[2]/div/div/form/dl/dd[2]/input')
    login_password_box.send_keys(login_password)

    sleep(0.5)

    login_button = driver.find_element(By.XPATH, '/html/body/div[2]/div/div/form/dl/dd[3]/button')
    login_button.click()

    sleep(2)

    search_box = driver.find_element(By.XPATH, '/html/body/div[2]/article/div[2]/section[1]/form/div/input')
    search_box.send_keys(cancer_name)

    sleep(0.5)

    submit_button = driver.find_element(By.XPATH, '/html/body/div[2]/article/div[2]/section[1]/form/div/button')
    submit_button.click()

    sleep(1)

    cancer_button = driver.find_element(By.XPATH, '/html/body/div[2]/div/div[2]/div[5]/div/table/tbody/tr[1]/td[1]/a')
    cancer_button.click()

    sleep(1)

    top20 = []
    for x in range(18, 95, 4):
        top20.append(
            driver.find_element(By.CSS_SELECTOR, f'#top20_census > div > svg > a:nth-child({x}) > text > tspan').text)
    top20_genes = [gene.split(' ')[0] for gene in top20]
    driver.quit()
    return top20_genes


def filter_csv(variant_file):
    variants = pd.read_csv(variant_file)
    exons = variants.loc[variants["Gene section"].str.contains("Exon")]
    exons_pm = exons.loc[exons["AA change"].notna()]
    exons_spm = exons_pm.loc[~exons_pm["AA change"].str.contains(r'\*', na=False)]
    gene_mutations_series = exons_spm.groupby('Gene ID')["AA change"].apply(list)
    sorted_gene_mutations_series = gene_mutations_series.sort_values(key=lambda x: x.map(len),
                                                                     ascending=False)  ##should be ascending = False when actually running
    print(sorted_gene_mutations_series)
    return sorted_gene_mutations_series


def process_one_gene(gene_target, all_mutations, parent_dir, flags, mhc_list, log_out, log_lin):
    """Fetch sequence, generate peptides & features, return raw DataFrames."""
    # 1) Get gene sequence
    log_lin.append(f"Getting gene sequence for {gene_target}...")
    log_out.text("\n".join(log_lin))
    get_gene_sequence(gene_target)
    log_lin.append(f"Sequence obtained for {gene_target}")
    log_out.text("\n".join(log_lin))

    dfs_raw = []
    immuno_script = parent_dir / "immunogenicity" / "predict_immunogenicity.py"

    for mhc_class in mhc_list:
        log_lin.append(f"Processing MHC Class {mhc_class} for {gene_target}...")
        log_out.text("\n".join(log_lin))

        # generate mutant FASTA & peptides
        log_lin.append(f"Making mutant FASTA proteins")
        log_out.text("\n".join(log_lin))
        make_mutant_genes(all_mutations, f"{gene_target}.fasta", parent_dir)
        log_lin.append(f"Getting all epitopes and binding affinities")
        log_out.text("\n".join(log_lin))
        mutant_all = get_epitopes_ba(all_mutations, mhc_class, parent_dir, log_out, log_lin)
        log_lin.append(f"Filtering for mutant epitopes")
        log_out.text("\n".join(log_lin))
        mutant_mut = get_mutant_epitopes(all_mutations, mhc_class, mutant_all, parent_dir)
        log_lin.append(f"Generating possible peptide sequences for all mutants")
        log_out.text("\n".join(log_lin))
        get_peptides(all_mutations, mutant_mut, mhc_class, parent_dir)

        log_lin.append(f"Obtaining epitope characteristics")
        log_out.text("\n".join(log_lin))
        # write raw epitopes + features
        # raw_file = parent_dir / f"{gene_target}_all_variables_mhc{mhc_class.lower()}.xlsx"
        # with pd.ExcelWriter(raw_file, engine='openpyxl') as writer:
        for pm in all_mutations:
            current_df = mutant_mut[pm].copy()
            current_df["gene"] = gene_target
            current_df["mutation"] = pm
            current_df["mhc class"] = mhc_class
            cols = ["gene", "mutation", "mhc class"] + [c for c in current_df.columns if
                                                        c not in ("gene", "mutation", "mhc class")]
            current_df = current_df[cols]
            print(current_df)
            # (i) Immunogenicity
            if flags['i']:
                log_lin.append(f"Obtaining {pm} immunogenicity...")
                log_out.text("\n".join(log_lin))
                if mhc_class == 'I':
                    current_df = get_local_immunogenicity_mhci(immuno_script,
                                                               parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt",
                                                               current_df)
                else:
                    fasta = (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.fasta").read_text()
                    peptides = [l.strip() for l in
                                (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt").read_text().splitlines()]
                    current_df = get_immunogenicity_mhcii(peptides, fasta, current_df)
                log_lin.append(f"Done immunogenicity for {pm}")
                log_out.text("\n".join(log_lin))

            # (ii) Antigenicity
            if flags['an']:
                log_lin.append(f"Obtaining {pm} antigenicity...")
                log_out.text("\n".join(log_lin))
                peptides = [l.strip() for l in
                            (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt").read_text().splitlines()]
                fasta_path = parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.fasta"
                current_df = get_antigenicity(peptides, fasta_path, current_df)
                log_lin.append(f"Done antigenicity for {pm}")
                log_out.text("\n".join(log_lin))

            # (iii) Allergenicity
            if flags['al']:
                log_lin.append(f"Obtaining {pm} allergenicity...")
                log_out.text("\n".join(log_lin))
                peptides = [l.strip() for l in
                            (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt").read_text().splitlines()]
                fasta = (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.fasta").read_text()
                if mhc_class == 'I':
                    current_df = get_allergenicity_algpred(peptides, fasta, current_df)
                else:
                    current_df = get_allergenicity_netallergen(peptides, fasta, current_df,
                                                               parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.fasta")
                log_lin.append(f"Done allergenicity for {pm}")
                log_out.text("\n".join(log_lin))

            # (iv) Toxicity
            if flags['t']:
                log_lin.append(f"Obtaining {pm} toxicity...")
                log_out.text("\n".join(log_lin))
                peptides = [l.strip() for l in
                            (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt").read_text().splitlines()]
                fasta = (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.fasta").read_text()
                current_df = get_toxicity(peptides, fasta, current_df)
                log_lin.append(f"Done toxicity for {pm}")
                log_out.text("\n".join(log_lin))

            # (v) ProtParam
            if any(flags[k] for k in ('h', 'ins', 'ali', 'iso', 'g')):
                log_lin.append(f"Obtaining protparam for {pm}...")
                log_out.text("\n".join(log_lin))
                peptides = [l.strip() for l in
                            (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt").read_text().splitlines()]
                current_df = get_protparam(
                    peptides,
                    flags['h'], flags['ins'], flags['ali'], flags['iso'], flags['g'],
                    current_df
                )
                log_lin.append(f"Done protparam for {pm}")
                log_out.text("\n".join(log_lin))

            # (vi) IFN-γ
            if flags['ifn']:
                log_lin.append(f"Obtaining {pm} IFN-γ predictions...")
                log_out.text("\n".join(log_lin))
                peptides = [l.strip() for l in
                            (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.txt").read_text().splitlines()]
                fasta = (parent_dir / "Sequences" / f"{pm}peptides_{mhc_class}.fasta").read_text()
                current_df = get_ifn(peptides, fasta, current_df)
                log_lin.append(f"Done IFN-γ for {pm}")
                log_out.text("\n".join(log_lin))

            # current_df.to_excel(writer, sheet_name=pm, index=False)
            dfs_raw.append(current_df)

    return dfs_raw


def postprocess_gene(mhci_df, mhcii_df, flags, mhc_list, log_out, log_lin):
    """Take raw DFs → ranking → filtering → population coverage Excel files."""
    # ensure directories exist
    results_dir = Path(os.getcwd()) / "results"
    results_dir.mkdir(exist_ok=True)

    df_map = {"I": mhci_df, "II": mhcii_df}

    for mhc_class in mhc_list:
        df = df_map[mhc_class]
        print(f'Goodbye\n\n{df}')
        if df is None:
            # log_lin.append(f"No data for MHC-{mhc_class}; skipping.")
            # log_out.text("\n".join(log_lin))
            continue

        # 1) Scoring
        if flags["scoring"]:
            log_lin.append(f"Normalizing & scoring MHC-{mhc_class}...")
            log_out.text("\n".join(log_lin))
            df = apply_scoring_function(normalize_data(df, mhc_class), mhc_class)

        # 2) Filtering
        if flags["filtering"]:
            log_lin.append(f"Filtering MHC-{mhc_class} epitopes...")
            log_out.text("\n".join(log_lin))
            df = get_filtered_epitopes(
                df,
                mhc_class,
                flags.get('h', False),
                flags.get('ins', False),
                flags.get('t', False),
                flags.get('ifn', False)
            )

        # 3) Write one final .xlsx with population coverage sheets if requested
        final_path = results_dir / f"final_aggregated_mhc_{mhc_class}.xlsx"
        with pd.ExcelWriter(final_path, engine='openpyxl') as writer:
            # main epitopes sheet
            df.to_excel(writer, sheet_name="Epitopes", index=False)

            if flags["pop_cov"]:
                log_lin.append(f"Calculating population coverage for MHC-{mhc_class}...")
                log_out.text("\n".join(log_lin))

                opt_df = get_optimized_epitopes(df, mhc_class)
                cov_reg, cov_opt = get_population_coverage(df, opt_df, mhc_class, results_dir)

                cov_reg.to_excel(writer, sheet_name="PopCov_Regular", index=False)
                cov_opt.to_excel(writer, sheet_name="PopCov_Optimized", index=False)

        log_lin.append(f"Wrote final file for MHC-{mhc_class}: {final_path.name}")
        log_out.text("\n".join(log_lin))

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w') as zf:
        for file in results_dir.rglob('*'):
            if file.is_file():
                zf.write(file, arcname=file.relative_to(results_dir))
    zip_buffer.seek(0)
    return zip_buffer


def run_personalized(sorted_gene_mutations_series, case, num_genes, cancer, i, an, al, ali, g, iso, h, ins, t, ifn,
                     filtering, scoring, pop_cov, mhc_classes, log_out, log_lin):
    """Main entry: select genes, call process_one_gene, then postprocess, then zip results."""
    if log_lin is None:
        log_lin = []

    # Gene selection
    if case:
        add20 = gene_checker_helper(cancer)
        genes = sorted_gene_mutations_series.loc[
            sorted_gene_mutations_series.index.isin(add20)
        ].to_dict()
    else:
        genes = dict(sorted_gene_mutations_series.head(num_genes))

    flags = {
        'i': i, 'an': an, 'al': al, 'ali': ali, 'g': g,
        'iso': iso, 'h': h, 'ins': ins, 't': t, 'ifn': ifn,
        'filtering': filtering, 'scoring': scoring, 'pop_cov': pop_cov
    }

    parent_dir = Path(os.getcwd())
    results_dir = parent_dir / 'results'
    results_dir.mkdir(exist_ok=True)

    log_lin.append(f"Starting to gather epitopes")
    log_out.text("\n".join(log_lin))

    genes = {"MUC16": ["V14466L"], "TTN": ["T33079A"]}
    # 1) Raw processing
    raw_results_mhci = []
    raw_results_mhcii = []
    if not an:
        for gene, muts in genes.items():
            dfs_raw = process_one_gene(gene, muts, parent_dir, flags, mhc_classes, log_out, log_lin)
            for df in dfs_raw:
                cls = df["mhc class"].iat[0]
                if cls == "I":
                    raw_results_mhci.append(df)
                else:
                    raw_results_mhcii.append(df)

    if raw_results_mhci:
        agg_mhc_i = pd.concat(raw_results_mhci, ignore_index=True)
        mhc_i_path = results_dir / "aggregated_mhc_I.xlsx"
        agg_mhc_i.to_excel(mhc_i_path, index=False)
        agg_mhc_i_deduped = agg_mhc_i.drop_duplicates(subset=["peptide"], keep="first")
        with open(results_dir / "all_peptides_anti_mhci.fasta", "w") as fh:
            for idx, pep in enumerate(agg_mhc_i_deduped["peptide"], start=1):
                fh.write(f">peptide_{idx}\n")
                fh.write(f"{pep}\n")

    if raw_results_mhcii:
        agg_mhc_ii = pd.concat(raw_results_mhcii, ignore_index=True)
        mhc_ii_path = results_dir / "aggregated_mhc_II.xlsx"
        agg_mhc_ii.to_excel(mhc_ii_path, index=False)
        agg_mhc_ii_deduped = agg_mhc_ii.drop_duplicates(subset=["peptide"], keep="first")
        with open(results_dir / "all_peptides_anti_mhcii.fasta", "w") as fh:
            for idx, pep in enumerate(agg_mhc_ii_deduped["peptide"], start=1):
                fh.write(f">peptide_{idx}\n")
                fh.write(f"{pep}\n")

    # # 2) Post‐processing
    # else:
    #     for gene, (df_raw, gene_dir) in aggregated_files.items():
    #         postprocess_gene(gene, df_raw, flags, mhc_classes, parent_dir, log_out, log_lin)

    # 3) Zip up results folder
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w') as zf:
        for file in results_dir.rglob('*'):
            if file.is_file():
                zf.write(file, arcname=file.relative_to(results_dir))
    zip_buffer.seek(0)
    return zip_buffer


def vaxijen_to_df(vaxijen_txt):
    # # 1) Read the text from a file
    # with open(vaxijen_txt, "r") as file:
    #     text = file.read()
    text = vaxijen_txt.read().decode("utf-8")

    # Initialize lists to store parsed data
    peptides = []
    antigenicity_scores = []

    # Regular expression patterns
    peptide_pattern = r'>peptide_\d+\s+([A-Z]+)'
    antigenicity_pattern = r'Overall Prediction for the Protective Antigen = ([-+]?\d*\.?\d+)'

    # Find all peptide sequences
    peptide_matches = re.findall(peptide_pattern, text)

    # Find all antigenicity scores
    antigenicity_matches = re.findall(antigenicity_pattern, text)

    # Ensure we have the same number of peptides and scores
    if len(peptide_matches) != len(antigenicity_matches):
        print(f"Warning: Found {len(peptide_matches)} peptides but {len(antigenicity_matches)} scores")

    # Create the DataFrame
    df = pd.DataFrame({
        'peptide': peptide_matches,
        'antigenicity': [float(score) for score in antigenicity_matches]
    })

    return df


st.markdown(
    "<h1 style='text-align: center;'>Personalized Cancer Vaccine Design</h1>",
    unsafe_allow_html=True
)

# define your epitope options in one place
epitope_options = [
    "Aliphatic Index",
    "GRAVY Score", "Isoelectric Point", "Half-Life",
    "Instability Index", "Toxicity", "IFN-γ Release"
]

tab1, tab2, tab3 = st.tabs([
    "① Generate Aggregated MHC Files",
    "② Process Aggregated Files",
    "③ Documentation"
])

with tab1:
    st.subheader("Step 1: From Variant File → Aggregated MHC-I/II .xlsx")
    gene_file = st.file_uploader(
        "Upload your Partek Flow variant CSV",
        type=["csv"],
        help="This CSV drives candidate epitope generation."
    )
    # only process once they upload
    if gene_file:
        try:
            sorted_genes = filter_csv(gene_file)
        except Exception as e:
            st.error(f"Error processing Partek Flow CSV: {e}")
            sorted_genes = None
    else:
        sorted_genes = None

    case_selection = st.radio(
        "Vaccine Design Module",
        [
            "Module 1: Identify Generalized Variants",  # Human Protein Atlas
            "Module 2: Identify Cancer-Specific Variants"  # COSMIC Database
        ]
    )
    case = (case_selection == "Module 2: Identify Cancer-Specific Variants")
    cancer_name = (
        st.text_input("Enter cancer name (e.g. Lung, Melanoma…)")
        if case else None
    )

    if not case and sorted_genes is not None:
        max_genes = len(sorted_genes)
        st.markdown(f"Max genes available: **{max_genes}**")
        num_genes = st.number_input(
            "How many genes to include",
            min_value=1, max_value=max_genes,
            value=1, step=1
        )
    else:
        num_genes = None

    mhc_classes = st.multiselect(
        "Choose MHC class(es)",
        ["Class I", "Class II"]
    )
    epitope_props_raw = st.multiselect(
        "Select additional epitope characteristics",
        epitope_options
    )

    # filter_raw = st.checkbox("Filter candidates by selected properties")
    # popcov_raw = st.checkbox("Calculate population coverage")

    if st.button("Generate .xlsx Files"):
        if not gene_file:
            st.warning("Please upload a variant CSV first.")
        elif case and not cancer_name:
            st.warning("Please enter a cancer name.")
        elif not mhc_classes:
            st.warning("Please select at least one MHC class.")
        else:
            # map UI labels to your run_personalized flags
            flags = {
                prop.lower().replace(" ", "_").replace("-", "_").replace("γ", "gamma"): (prop in epitope_props_raw)
                for prop in epitope_options
            }
            print(flags)
            # normalize class selection
            if "Class I" in mhc_classes and "Class II" in mhc_classes:
                mhc_list = ["I", "II"]
            elif "Class I" in mhc_classes:
                mhc_list = ["I"]
            else:
                mhc_list = ["II"]
            print(mhc_list)
            with st.expander("Output Log", expanded=True):
                log_out = st.empty()
                log_lines = []

                out_files = run_personalized(
                    sorted_genes,
                    case=case,
                    cancer=cancer_name,
                    num_genes=num_genes,
                    i=True,
                    an=False,  # no antigenicity yet
                    al=True,
                    ali=flags["aliphatic_index"],
                    g=flags["gravy_score"],
                    iso=flags["isoelectric_point"],
                    h=flags["half_life"],
                    ins=flags["instability_index"],
                    t=flags["toxicity"],
                    ifn=flags["ifn_gamma_release"],
                    filtering=False,
                    scoring=False,
                    pop_cov=False,
                    mhc_classes=mhc_list,
                    log_out=log_out,
                    log_lin=log_lines
                )

            st.download_button(
                "Download Aggregated Files (ZIP)",
                data=out_files,
                file_name="aggregated_mhc_files.zip",
                mime="application/zip"
            )

with tab2:
    st.subheader("Step 2: Upload Aggregated .xlsx + Antigenicity → Ranking, Filtering & Population Coverage")
    antigenicity_files = st.file_uploader(
        "Upload new antigenicity values in TXT",
        type=["txt"],
        accept_multiple_files=True,
        help="Optional: if you have updated antigenicity data"
    )
    antigenicity_df_mhci = None
    antigenicity_df_mhcii = None
    if antigenicity_files:
        for txt in antigenicity_files:
            name = txt.name.lower()
            if "mhc_i" in name:
                antigenicity_df_mhci = vaxijen_to_df(txt)
            elif "mhc_ii" in name:
                antigenicity_df_mhcii = vaxijen_to_df(txt)
            else:
                st.warning(f"Could not detect MHC class in filename '{txt.name}'")
    antigenicity_dfs = [antigenicity_df_mhci, antigenicity_df_mhcii]

    mhc_xlsx_files = st.file_uploader(
        "Upload aggregated epitope files in XLSX",
        type=["xlsx"],
        accept_multiple_files=True,
        help="Required to update with antigenicity, filtering & ranking"
    )
    agg_df_mhci = pd.DataFrame()
    agg_df_mhcii = pd.DataFrame()
    if mhc_xlsx_files:
        for xlsx in mhc_xlsx_files:
            temp_df = pd.read_excel(xlsx, engine="openpyxl")
            try:
                cls = temp_df["mhc class"].iat[0]
            except Exception as e:
                st.error(f"Could not determine MHC class in '{xlsx.name}': {e}")
                continue

            if cls == "I":
                agg_df_mhci = pd.concat([agg_df_mhci, temp_df], ignore_index=True)
            elif cls == "II":
                agg_df_mhcii = pd.concat([agg_df_mhcii, temp_df], ignore_index=True)
            else:
                st.warning(f"Unknown MHC class '{cls}' in '{xlsx.name}'")
    agg_dfs = [agg_df_mhci, agg_df_mhcii]

    agg_anti_df_mhci = None
    agg_anti_df_mhcii = None
    if antigenicity_files and mhc_xlsx_files:
        if len(antigenicity_dfs) == len(agg_dfs):
            if (antigenicity_df_mhci is not None) and (len(agg_df_mhci) != 0):
                agg_anti_df_mhci = pd.merge(agg_df_mhci, antigenicity_df_mhci, on="peptide", how="left")
                print(agg_anti_df_mhci)
            elif (antigenicity_df_mhcii is not None) and (len(agg_df_mhcii) != 0):
                agg_anti_df_mhcii = pd.merge(agg_df_mhcii, antigenicity_df_mhcii, on="peptide", how="left")
                print(agg_anti_df_mhcii)
            else:
                st.warning(f"MHC Classes of antigenicity files and aggregate .xlsx files do not match. "
                           f"Please submit matching files.")
        else:
            st.warning(f"Please upload the same number of antigenicity files and aggregated .xlsx files")

    epitope_props_agg = st.multiselect(
        "Which additional characteristics exist in these .xlsx files?",
        epitope_options
    )
    # use_ml_scoring = False
    # if all(p in epitope_props_agg for p in ["Immunogenicity", "Allergenicity"]):
    use_ml_scoring = st.checkbox("Use machine learning-driven epitope scoring")
    filter_agg = st.checkbox("Filter epitopes")
    popcov_agg = st.checkbox("Calculate population coverage")

    if st.button("Process Aggregated Files"):
        if not mhc_xlsx_files:
            st.warning("Please upload the aggregated .xlsx files first.")
        else:
            flags = {
                prop.lower().replace(" ", "_").replace("γ", "gamma"): (prop in epitope_props_agg)
                for prop in epitope_options
            }
            flags.update({
                "scoring": use_ml_scoring,
                "filtering": filter_agg,
                "pop_cov": popcov_agg
            })

            mhc_list = []
            if agg_anti_df_mhci is not None:
                mhc_list.append("I")
            if agg_anti_df_mhcii is not None:
                mhc_list.append("II")

            with st.expander("Output Log", expanded=True):
                log_out = st.empty()
                log_lines = []

                out_files = postprocess_gene(
                    agg_anti_df_mhci,
                    agg_anti_df_mhcii,
                    flags,
                    mhc_list,
                    log_out, log_lines)

            st.download_button(
                "Download Final Output (ZIP)",
                data=out_files,
                file_name="final_vaccine_outputs.zip",
                mime="application/zip"
            )

with tab3:
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

