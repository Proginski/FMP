#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 16:47:03 2021

@author: paul.roginski
"""
import pandas as pd
import csv

csv = pd.read_csv("/home1/paul.roginski/WD/Archaea_Data/csv_ist_withGC.csv")

for group in set(csv['csvs_GC']):

    print(group)
    
    # csv_to_concat = [pd.read_csv("/home1/paul.roginski/WD/" + str(csv['csvs'][i])) for i in range(0,len(csv['csvs'])) if csv['csvs_GC'][i] == group]
    csv_to_concat = ["/home1/paul.roginski/WD/" + str(csv['csvs'][i]) for i in range(0,len(csv['csvs'])) if csv['csvs_GC'][i] == group]        
    
    
    df = pd.concat((pd.read_csv(f, header = 0) for f in csv_to_concat))
    df_deduplicated = df.drop_duplicates()
    df_deduplicated.to_csv(str(group)+".csv")