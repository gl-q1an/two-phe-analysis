# DESCRIPTION:This script is used to automatically transform the imported GWAS summary into 
#             those columns we want, and try to calculate those columns we want

import os
import numpy as np
import pandas as pd
import json
import logging

Header = "=============================================\n"
Header +="||                                         ||\n"
Header +="||      Automated GWAS Summary Converter   ||\n"
Header +="||                                         ||\n"
Header +="=============================================\n"

describe_dict = {
    "args_all":Header,
    "summary":"",
    "out":"",
    "outform":"",
    "N_exa":"",
    "Nca_exa":"",
    "Nco_exa":"",
    "N_order":"",
    "atuo":""
}

def identify_json(jsonfile):
    with open(jsonfile, 'r') as json_file:
        load = json.load(json_file)
        ref_list = {key: value for key, value in load.items() if "DESCRIPTION" not in key}
        res_list = {key: None for key in ref_list}
    return(ref_list,res_list)

if __name__=="__main__":
    exit()