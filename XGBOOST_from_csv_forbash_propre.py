#!/usr/bin/env python3
# -*- coding: latin-1 -*-
"""
Created on Wed Jul  7 11:09:08 2021

@author: paul.roginski

This script takes as input a csv file containing a 'category' column with the 
classes to predict and all the necessary features.

csv_path : specify the path to a csv file to make prediction with
kfold : (boolean) if k-fold cv should be performed or not
features_no_train : avoid using certain features for trainning. Especially
usefull for simple prediction where you can visualize these features in the 
output for manual analysis.
"""


import pandas as pd



def main(csv_path,
         kfold = True,
         features_no_train = ['category']):
    
    
    dataset = pd.read_csv(csv_path)
    
    # The 'category' column must not be used for trainning
    if 'category' not in features_no_train :
        features_no_train.append('category')
        

    ################## Prediction Algorithm ##################
    
    from xgboost import XGBClassifier
    
    # Split data into X and y
    X = dataset.drop(features_no_train, axis=1)
    Y = dataset['category']
    
    if kfold :
        # stratified k-fold cross validation evaluation of xgboost model
        from sklearn.model_selection import StratifiedKFold
        from sklearn.model_selection import cross_val_score    
    
        # TRAIN THE XGBOOST MODEL
        # fit model to training data
        model = XGBClassifier(use_label_encoder=False, eval_metric='logloss')
        kfold = StratifiedKFold(n_splits=10)
        
        # return(kfold)
        results = cross_val_score(model, X, Y, cv=kfold)
        print("Accuracy: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))

        print('DONE.')
        return(results.mean())
       
        
    else :
    
        from sklearn.model_selection import train_test_split
        # split data into train and test sets
        test_size = 0.2
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size)

        # TRAIN THE XGBOOST MODEL
        # fit model to training data
        model = XGBClassifier(use_label_encoder=False, eval_metric='logloss')
        model.fit(X_train, y_train)
       
        # To see the parameters used in a trained model:
        # print(model)
        
        # MAKE PREDICTIONS
        # Proba between 0 and 1
        y_pred = model.predict_proba(X_test)
        # Integer 0 OR 1
        predictions = [round(value) for value in y_pred[:,1]]
        
        # evaluate predictions
        from sklearn.metrics import accuracy_score
        accuracy = accuracy_score(y_test, predictions)
        print("Accuracy: %.2f%%" % (accuracy * 100.0))
        
    
    ##########################################################################
        
    
        # Generate the ouput based on the test sequences
        results = dataset.loc[X_test.index]
        results.insert(0, 'Prediction', predictions)
        
        # Generate the 'constistency' column
        def conditions(s):
            if (s['Prediction'] == 0) :
                if (s['category'] == 0) :
                    return "T_cds"
                else : return "F_cds"
            else:
                if (s['category'] == 1) :
                    return "T_igorf"
                else : return "F_igorf"
    
        results.insert(2, 'consistency', results.apply(conditions, axis=1))
    
    
        print('DONE.')
        return(results)




csv_path = "D:/WD/Prediciton/Data/ScerCDSnRAND.csv"

output = main(csv_path = csv_path,
              kfold = False,
              features_no_train=['Name',"seq_length_nucl","LysArg","GC_content","AT_content"])

# For bash use
# import sys
# output = main(csv_path = sys.argv[1],
#               features_no_train=['Name',"seq_length_nucl","LysArg","GC_content","AT_content"])



