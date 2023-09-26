import os, sys

import numpy as np
import pandas as pd

import sklearn
import sklearn.ensemble
import sklearn.metrics 
import sklearn.externals
import sklearn.model_selection

from argparse import ArgumentParser

from pathlib import Path
from typing import Optional, List

import xgboost
import gc

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

inputs = [
    'fractional_subjet_pt', 'min_pairwise_subjets_mass', 'mass', 'nsubjets', 'tau3_over_tau2', 'pt'
]

xgboost_params = {
    'n_estimators': 500, 'learning_rate': .1, 'max_depth': 4, 'min_child_weight': 0.01
}

pT_distribution_plot = False

def loading_and_processing_h5(file):
    
    df = pd.read_hdf(file)
    df = df.fillna(False)
    # skimming the df sgn/bkg -> pure hadronic jets/pure qcd   
    if 'ttX_mass' in str(file):
        # for signal: has_top_fromResonance, has_hadronicTop_topIsInside -> pure hadronic and jet from resonance
        df = df[df['has_top_fromResonance'] & df['has_hadronicTop_topIsInside'] ]
        df = df[inputs]
        df['label'] = 1
    else:
        if 'has_other' not in df.columns: df['has_other'] = False
        # for background: 
        #   -'has_noTopDaughters_topIsNotInside_and_has_gluon_or_quark_not_fromTop' -> pure QCD in case of samples with top
        #   -'has_gluon_or_quark_not_fromTop', 'has_other'  -> pure QCD in case of pure QCD samples
        df = df[df['has_gluon_or_quark_not_fromTop'] | df['has_other']]
        df = df[inputs]
        df['label'] = 0
    return df    

def pT_plot(df, path):

    # Separate the data based on the label
    pt_sgn = df[df["label"] == 1]["pt"]
    pt_bkg = df[df["label"] == 0]["pt"]

    # Plot
    plt.figure(figsize=(10, 6))

    plt.hist(pt_sgn, bins=30, alpha=0.5, label=f"pure hadronic top: {len(pt_sgn)}")
    plt.hist(pt_bkg, bins=30, alpha=0.5, label=f"pure QCD : {len(pt_bkg)}")

    plt.xlabel('HOTVR pT')
    plt.ylabel('HOTVR')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()

    # Save the figure
    fpath = path / 'pT_distribution_1_vs_0.png'
    plt.savefig(fpath, dpi=300)
    print(f"Saving file: {fpath}")

def creating_training_data(df, path):
    features_array, label_array = [], []

    label_array = df['label'].values.reshape(-1, 1)
    features_array = df[[input for input in inputs if input != 'pt']].values # to be modified such to avoid 'pt' in the dataframe

    print(f"Total samples: {features_array.shape}")
    # only 20% for training
    features_array, _, label_array, _, = sklearn.model_selection.train_test_split(
        features_array,label_array, test_size=0.8,shuffle=True, random_state=76177
    )

    path = path / 'trainBDT_outputs'
    path.mkdir(parents=True, exist_ok=True)
    path.chmod(0o777)
    fpath = path / 'trainData.npz'
    np.savez_compressed(fpath,features=features_array,labels=label_array)

def draw_loss(path, results):
    plt.figure(figsize=[7.3, 6.5])
    plt.plot(np.arange(xgboost_params['n_estimators']),results['validation_0']['logloss'],color='darkgray',linewidth=2, linestyle='-',label='Train')
    plt.plot(np.arange(xgboost_params['n_estimators']),results['validation_1']['logloss'],color='royalblue',linewidth=3, linestyle='--',label='Test')
    plt.ylabel("log(Loss)")
    plt.xlabel("#Estimators")
    plt.tight_layout(rect=[0., 0, 1, 0.91])
    plt.legend(ncol=2, bbox_to_anchor=(0, 1.0), loc="lower left")
    plt.savefig(path)
    plt.close()
    print(f"Saving file: {path}")
    path.chmod(0o777)

def draw_roc(path, rates):
    plt.figure(figsize=[7.3, 6.5])
    plt.plot(rates['true_positive_train'], rates['false_positive_train'], color='darkgray',linewidth=2, linestyle='-',label="Train")
    plt.plot(rates['true_positive_test'], rates['false_positive_test'], color='royalblue',linewidth=3, linestyle='--',label="Test")
    bkgIdx = np.argmin(np.abs(rates['false_positive_test']-0.1))
    plt.plot(
        [rates['true_positive_test'][bkgIdx],rates['true_positive_test'][bkgIdx]],[1e-3,rates['false_positive_test'][bkgIdx]],
        color='black',linewidth=1, linestyle='-'
        )
    plt.xlim([0.0,1.0])
    plt.ylim([1e-3,1.0])
    plt.yscale('log')
    plt.xlabel("Signal efficiency")
    plt.ylabel("Background efficiency")
    plt.gca().xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    plt.gca().yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    plt.grid(True, which='both', axis='both', color='black', linestyle='--', linewidth=1)
    plt.tight_layout(rect=[0., 0, 1, 0.91])
    plt.legend(ncol=2, bbox_to_anchor=(0, 1.0), loc="lower left")
    plt.savefig(path)
    plt.close()
    print(f"Saving file: {path}")
    path.chmod(0o777)

def draw_score(path, predictions):
    plt.figure(figsize=[7.3, 6.5])
    plt.hist(predictions['train_bkg'], bins=40, range=[0,1], histtype='stepfilled', color='royalblue', label="Train bkg", alpha=0.6)
    plt.hist(predictions['train_sgn'], bins=40, range=[0,1], histtype='stepfilled', color='orange', label="Train sig", alpha=0.6)
    plt.hist(predictions['test_bkg'], bins=40, range=[0,1], histtype='step', color='cornflowerblue', linewidth=2, linestyle='--', label="Test bkg")
    plt.hist(predictions['test_sgn'], bins=40, range=[0,1], histtype='step', color='gold', linewidth=2, linestyle='--', label="Test sig")
    plt.xlim([0.0,1.0])
    plt.xlabel("Score")
    plt.ylabel("Normalized events")
    plt.tight_layout(rect=[0., 0, 1, 0.91])
    plt.legend(ncol=2, bbox_to_anchor=(0, 1.0), loc="lower left")
    plt.savefig(path)
    plt.close()
    print(f"Saving file: {path}")
    path.chmod(0o777)

def main(
        input_dir,
        output_dir,
        year):

    # Ensure input/output dirs are Path instances.
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    dfs = pd.DataFrame()

    ######### LOADING/PROCESSING THE DATAFRAME
    # Given the large number of files to be processed, the total df (all jets with inputs/labels) is available as 'all_jets.h5'
    # if it is not, the total df is generated and saved as 'all_jets.h5'
    if os.path.exists(input_dir / f"total_jets_{year}_h5" / 'all_jets.h5'):
        # reading all_jets.h5 if exists
        dfs = pd.read_hdf(input_dir / f"total_jets_{year}_h5" / 'all_jets.h5')
    else:
        # Glob the files, and convert to a sorted list.
        # sgn and bkg files in separated list so it is easy to filter them when needed
        fnames_sgn = sorted(input_dir.glob(f"sgn_{year}_hotvr_h5/*/*.h5"))
        fnames_bkg = sorted(input_dir.glob(f"bkg_{year}_hotvr_h5/qcd*/*.h5"))  # preliminary study with only QCD files 
        fnames =  fnames_sgn + fnames_bkg
        print(f"Found {len(fnames)} files in {input_dir}")

        for idx_file, fname in enumerate(fnames):
            # print(
            #     "Done %4.1f%% - processing: %s",
            #     idx_file / len(fnames) * 100, fname)

            df = loading_and_processing_h5(fname)
            if df is None:
                # If the file is empty, skip it.
                continue
            dfs = dfs.append(df)

        # Shuffle the dataframe
        dfs = dfs.sample(frac=1).reset_index(drop=True)
        
        # Saving the dataframe of all jets after the processing
        # so it can be read faster as input h5 
        input_dir = input_dir / f"total_jets_{year}_h5"
        input_dir.mkdir(parents=True, exist_ok=True)
        input_dir.chmod(0o777)

        fpath = input_dir /  'all_jets.h5'
        dfs.to_hdf(fpath, key='df', mode='w')
        print(f"Saving file: {fpath}")
        fpath.chmod(0o777)

    if dfs.empty: 
        print("Unable to find 'all_jets.h5' or unable to process the individual files")
        sys.exit()
    
    print('### Successfully loaded the h5 with all the jets')
    print(dfs)
    #########

    # plotting purpose: pT distribution of the two classes
    if pT_distribution_plot:
        output_dir = output_dir / 'pT_distribution_sgnVSbkg'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_dir.chmod(0o777)

        pT_plot(dfs, output_dir)
    # dfs = dfs.drop(columns='pt') # to be modified such to avoid 'pt' in the dataframe
    inputs.remove('pt') # to be modified such to avoid 'pt' in the input list

    ######### TRAINING-TEST SPLITTING
    # create the training set if does not exist
    print('\n### Train-test splitting')
    if os.path.exists(output_dir / 'trainBDT_outputs' / 'trainData.npz')==0:
        creating_training_data(dfs, output_dir)
    else: 
        output_dir = output_dir / 'trainBDT_outputs'
        fpath = output_dir / 'trainData.npz'

        loaded = np.load(fpath)
        features_array = loaded['features']
        label_array = loaded['labels']

    print(f"Split samples {features_array.shape}")

    trainFeatures,testFeatures,trainLabel,testLabel, = sklearn.model_selection.train_test_split(
        features_array, label_array, test_size=0.2, shuffle=True, random_state=24155
    )

    # use test=train for testing
    testFeatures = trainFeatures
    testLabel = trainLabel

    # ensure a one-dimensional array (*, )
    trainLabel = trainLabel.ravel()
    testLabel = testLabel.ravel()

    print(f"Train/test: {trainFeatures.shape}/{testFeatures.shape}")
    # for i, inputName in enumerate(inputs):
    #     print(inputName,trainFeatures[0:10,i])

    del features_array, label_array
    gc.collect()
    #########

    ######### XGBOOST 
    print(f"\nXGBOOST PARAMETERS: {xgboost_params}")

    bdt = xgboost.XGBClassifier(
     #loss='deviance', #=for logistic regression
     n_estimators=xgboost_params['n_estimators'], 
     learning_rate=xgboost_params['learning_rate'],
     max_depth=xgboost_params['max_depth'], 
     min_child_weight=xgboost_params['min_child_weight'],
     use_label_encoder = False,
     verbose=2,
     objective='binary:logistic',
     n_jobs=1,
     subsample=0.75, # default value... it might need to be changed
     importance_type='gain',
     booster='gbtree',
     tree_method='approx', 
    )

    bdt.fit(
        trainFeatures, trainLabel, 
        eval_metric=["error", "logloss"], 
        eval_set=[(trainFeatures, trainLabel), (testFeatures, testLabel)], 
        verbose=False
    )

    output_dir = output_dir / 'trainBDT_outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir.chmod(0o777)
    fpath = output_dir / 'hadronicTopVSQCD_bdt.bin'
    bdt._Booster.save_model(str(fpath))

    results = bdt.evals_result()

    trainPred = bdt.predict_proba(trainFeatures)[:,1]
    testPred = bdt.predict_proba(testFeatures)[:,1]
    fprTrain, tprTrain, thresTrain = sklearn.metrics.roc_curve(trainLabel, trainPred)
    fprTest, tprTest, thresTest = sklearn.metrics.roc_curve(testLabel, testPred)
    rates = {
        'false_positive_train': fprTrain, 'true_positive_train': tprTrain,
        'false_positive_test': fprTest, 'true_positive_test': tprTest,
    }

    predTrainBkg = trainPred[trainLabel<0.5]
    predTrainSig = trainPred[trainLabel>0.5]
    predTestBkg = testPred[testLabel<0.5]
    predTestSig = testPred[testLabel>0.5]
    predictions = {
        'train_bkg': predTrainBkg, 'train_sgn': predTrainSig,
        'test_bkg': predTestBkg, 'test_sgn': predTestSig,
    }


    aucTrain, aucTest = sklearn.metrics.auc(fprTrain, tprTrain), sklearn.metrics.auc(fprTest, tprTest)
    accTrain = sklearn.metrics.accuracy_score(trainLabel, 1.*(trainPred>0.5), normalize=True)
    accTest = sklearn.metrics.accuracy_score(testLabel, 1.*(testPred>0.5), normalize=True)
    print(f"\nAUC: {aucTrain} {aucTest}")
    print(f"Accuracy: {accTrain} {accTest}\n")

    for bkgFrac in [0.3, 0.2, 0.1, 0.05, 0.01]:
        bkgIdx = np.argmin(np.abs(fprTest-bkgFrac))
        print(f"Bkg {100.*bkgFrac}% -> Sigeff = {100.*tprTest[bkgIdx]}, Thres = {thresTest[bkgIdx]}")
        
    fpath = output_dir / 'hadronicTopVSQCD_stat.dat'
    fstat = open(fpath,'w')
    fstat.write("%.3e; %.3e; %.5f; %.5f; %.3f; %.3f \n"%(
        results['validation_0']['logloss'][-1],
        results['validation_1']['logloss'][-1],
        aucTrain, aucTest,
        100.*accTrain, 100.*accTest
    ))
    fstat.close()

    fpath = fpath = output_dir / 'hadronicTopVSQCD_lr.png'
    draw_loss(fpath, results)

    fpath = fpath = output_dir / 'hadronicTopVSQCD_roc.png'
    draw_roc(fpath, rates)

    fpath = fpath = output_dir / 'hadronicTopVSQCD_dist.png'
    draw_score(fpath, predictions)

def parse_args(argv: Optional[List[str]] = None):
    parser = ArgumentParser()

    parser.add_argument('--input_dir', type=Path, required=True,
        help="Input directory, where to find the h5 files")
    parser.add_argument('--output_dir', type=Path,
        help="Top-level output directory. "
             "Will be created if not existing. "
             "If not provided, takes the input dir.")
    parser.add_argument('--year', type=int, required=True,
        help='Year of the samples.')

    args = parser.parse_args(argv)

    # If output directory is not provided, assume we want the output to be
    # alongside the input directory.
    if args.output_dir is None:
        args.output_dir = args.input_dir

    # Return the options as a dictionary.
    return vars(args)


if __name__ == "__main__":
    main(**parse_args())
