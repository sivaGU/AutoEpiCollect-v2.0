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


