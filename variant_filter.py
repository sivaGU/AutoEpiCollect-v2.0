import pandas as pd
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver import ChromeOptions
import chromedriver_autoinstaller
from time import sleep


def gene_filter_helper(gene, cancer, filter_option):
    if filter_option == 1:
        options = webdriver.ChromeOptions()
        options.add_argument("--headless=new")
        driver = webdriver.Chrome()
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
        if specificity_box.contains("carcinoma", case=False):
            return True
        else:
            return False
    else:
        options = webdriver.ChromeOptions()
        options.add_argument("--headless=new")
        driver = webdriver.Chrome()
        driver.get('https://cancer.sanger.ac.uk/cosmic')

        sleep(1)

        search_box = driver.find_element(By.XPATH, '/html/body/div[1]/article/div[2]/section[1]/form/div/input')
        search_box.send_keys(cancer)

        submit_button = driver.find_element(By.XPATH, '/html/body/div[1]/article/div[2]/section[1]/form/div/button')
        submit_button.click()

        sleep(1.5)

        select_cancer = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[2]/div[5]/div/table/tbody/tr[1]/td[1]/a')
        select_cancer.click()

        sleep(2)

        genome_browser = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[3]/div/div[2]/section[2]/h2/a')
        driver.execute_script("arguments[0].scrollIntoView();", genome_browser)
        sleep(1)

        genome_area = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[3]/div/div[1]/nav/ul/li[1]/a')
        genome_area.click()
        sleep(1)

        top_20 = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[3]/div/div[2]/section[1]/div/div[2]/div[1]/div/svg/a[1]')
        sleep(1)


chromedriver_autoinstaller.install()
variants = pd.read_csv("variant_list.csv")
exons = variants.loc[variants["Gene section"].str.contains("Exon")]
exons_pm = exons.loc[exons["AA change"].notna()]
exons_spm = exons_pm.loc[~exons_pm["AA change"].str.contains(r'\*', na=False)]
gene_mutations_series = exons_spm.groupby('Gene ID')["AA change"].apply(list)
sorted_gene_mutations_series = gene_mutations_series.sort_values(key=lambda x: x.map(len), ascending=False)
print(sorted_gene_mutations_series)
print(sorted_gene_mutations_series.map(len))
gene_filter_helper(sorted_gene_mutations_series.index[0], "breast cancer", 2)



