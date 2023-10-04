import os, sys, shutil
import time

from multiprocessing import Pool

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

from scipy.interpolate import interp1d

from resample import multi_dim_resample

inputs = [
    'fractional_subjet_pt', 'min_pairwise_subjets_mass', 'mass', 'nsubjets', 'tau3_over_tau2', 'pt'
]

pT_distribution_plot = False
resampling = True

binning=np.arange(0.,10000,100) #np.logspace(0,4,100)

def read_and_filter(file):
    
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

def loading_and_processing_h5(files, output_dir, chunk_size=50):
    # process the files in chunks and save each merged DataFrame to a temporary .h5 file
    current_chunk = 1
    total_chunks = len(files) // chunk_size + (1 if len(files) % chunk_size else 0)
    
    for idx_file in range(0, len(files), chunk_size):
        print(f"Processing files-chunk n.{current_chunk} out of {total_chunks}")
        file_list = files[idx_file:idx_file+chunk_size]

        # creating a pool of processes to speed up the loading/processing of the .h5 files
        with Pool() as p:
            dfs = p.map(read_and_filter, file_list)
        
        merged_df = pd.concat(dfs, ignore_index=True)
        # Save the merged DataFrame to a temporary .h5 file
        temp_dir = output_dir / 'temp' 
        temp_dir.mkdir(parents=True, exist_ok=True)
        # temp_dir.chmod(0o777)
        fpath = temp_dir / f'temp_{idx_file}.h5'
        merged_df.to_hdf(fpath, key='df')

        current_chunk += 1

    # Concatenate the temporary .h5 files into a single file
    temp_files = [temp_dir / f'temp_{idx_file}.h5' for idx_file in range(0, len(files), chunk_size)]
    dfs = [pd.read_hdf(file) for file in temp_files]
    df = pd.concat(dfs, ignore_index=True)
    df = df.fillna(False)

    # Save the final DataFrame to a .h5 file
    fpath = output_dir / 'all_jets.h5'
    df.to_hdf(fpath, key='df')

    # Delete the temp directory
    shutil.rmtree(temp_dir)

    return

def pT_plot(ptBkg, ptSig, fpath, binning=binning):

    plt.figure()
    plt.hist(ptBkg,bins=binning,alpha=0.5,label=f"pure QCD: {len(ptBkg)}", density=True)
    plt.hist(ptSig,bins=binning,alpha=0.5,label=f"pure hadronic top : {len(ptSig)}", density=True, histtype='step', color='red', linewidth=2)
    plt.legend()
    # plt.xscale('log')
    plt.xlim(0,2000)

    plt.xlabel('HOTVR pT')
    # plt.ylabel('HOTVR')
    plt.ylabel('Normalized Counts')
    plt.legend(loc='upper right')
    # plt.grid(True)
    plt.tight_layout()

    # Save the figure
    plt.savefig(fpath, dpi=300)
    print(f"Saving file: {fpath}")

def cut_based_tagger_ROC(features, label):
    frac_pt_thr, min_mass_pair_thr, mass_thr, nsubj_thr, tau3VStau2_thr = 0.8, 50, [140.,220.], 3, 0.56

    true_positive_array = features[
        (features[:,0]<frac_pt_thr) & (features[:,1]>min_mass_pair_thr)  
        & (features[:,2]>min(mass_thr)) & (features[:,2]<max(mass_thr)) 
        & (features[:,3]>=nsubj_thr) & (features[:,4]<tau3VStau2_thr) & (label==1)
    ]
    all_true = label[label==1]

    false_positive_array = features[
        (features[:,0]<frac_pt_thr) & (features[:,1]>min_mass_pair_thr) 
        & (features[:,2]>min(mass_thr)) & (features[:,2]<max(mass_thr)) 
        & (features[:,3]>=nsubj_thr) & (features[:,4]<tau3VStau2_thr) & (label==0)
    ]
    all_false = label[label==0]

    print('\n### Cut-based efficiencies')
    print(f"Signal efficiency: {true_positive_array.shape[0]/all_true.shape[0]}; Background efficiency: {false_positive_array.shape[0]/all_false.shape[0]}")

    return true_positive_array.shape[0]/all_true.shape[0], false_positive_array.shape[0]/all_false.shape[0]

def creating_training_data(arrays, fpath):

    np.savez_compressed(fpath,features=arrays[0],labels=arrays[1])

    return

def draw_loss(path, results, **kwargs):
    plt.figure(figsize=[7.3, 6.5])
    plt.plot(np.arange(kwargs['n_estimators']),results['validation_0']['logloss'],color='darkgray',linewidth=2, linestyle='-',label='Train')
    plt.plot(np.arange(kwargs['n_estimators']),results['validation_1']['logloss'],color='royalblue',linewidth=3, linestyle='--',label='Test')
    plt.ylabel("log(Loss)")
    plt.xlabel("#Estimators")
    # plt.tight_layout(rect=[0., 0, 1, 0.91])
    plt.legend(ncol=2, bbox_to_anchor=(0, 1.0), loc="lower left")
    plt.savefig(path)
    plt.close()
    print(f"Saving file: {path}")
    path.chmod(0o777)

def draw_roc(path, rates, **kwargs):
    plt.figure(figsize=[7.3, 6.5])
    plt.plot(rates['true_positive_train'], rates['false_positive_train'], color='darkgray',linewidth=2, linestyle='-',label="Train")
    plt.plot(rates['true_positive_test'], rates['false_positive_test'], color='royalblue',linewidth=3, linestyle='--',label="Test")
    if 'cut_based_efficiencies' in kwargs.keys(): 
        plt.scatter(kwargs['cut_based_efficiencies'][0], kwargs['cut_based_efficiencies'][1], color='red', marker='*', s=150, label='Cut based tagger')   
    
    f = interp1d(rates['true_positive_test'], rates['false_positive_test'], kind='linear', fill_value='extrapolate')
    y_train_at_cut_based_eff = f(kwargs['cut_based_efficiencies'][0]).item()
    # print(y_train_at_cut_based_eff)

    # Add the y-values to the y-ticks - DOES NOT WORK ?!
    # yticks = list(plt.gca().get_yticks())
    # yticks.extend(list(rates['false_positive_train'].flatten()))
    # yticks.extend([y_train_at_cut_based_eff, kwargs['cut_based_efficiencies'][1]])
    # yticks = sorted(list(set(yticks)))
    # special_ticks = list(rates['false_positive_train'].flatten()) + [y_train_at_cut_based_eff, kwargs['cut_based_efficiencies'][1]]
    # plt.yticks(yticks, [f'{tick:.3f}' if tick in special_ticks else '' for tick in yticks])

    # Draw horizontal lines
    plt.axhline(y=y_train_at_cut_based_eff, color='red', linestyle='--')
    plt.axhline(y=kwargs['cut_based_efficiencies'][1], color='red', linestyle='--')

    plt.xlim([0.0,1.0])
    plt.ylim([1e-3,1.0])
    plt.yscale('log')
    plt.xlabel("Signal efficiency")
    plt.ylabel("Background efficiency")
    # plt.gca().xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    # plt.gca().yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    # plt.grid(True, which='both', axis='both', color='black', linestyle='--', linewidth=1)
    plt.tight_layout(rect=[0., 0, 1, 0.91])
    plt.legend(ncol=2, bbox_to_anchor=(0, 1.0), loc="lower left")

    plt.savefig(path)
    plt.close()
    print(f"Saving file: {path}")
    path.chmod(0o777)

def draw_score(path, predictions):
    plt.figure(figsize=[7.3, 6.5])
    plt.hist(predictions['train_bkg'], bins=40, range=[0,1], histtype='stepfilled', color='cornflowerblue', label="Train bkg", alpha=0.6)
    plt.hist(predictions['train_sgn'], bins=40, range=[0,1], histtype='stepfilled', color='gold', label="Train sig", alpha=0.6)
    plt.hist(predictions['test_bkg'], bins=40, range=[0,1], histtype='step', color='royalblue', linewidth=2, linestyle='--', label="Test bkg")
    plt.hist(predictions['test_sgn'], bins=40, range=[0,1], histtype='step', color='orange', linewidth=2, linestyle='--', label="Test sig")
    plt.xlim([0.0,1.0])
    plt.xlabel("Score")
    plt.ylabel("Normalized events")
    # plt.tight_layout(rect=[0., 0, 1, 0.91])
    plt.legend(ncol=2, bbox_to_anchor=(0, 1.0), loc="lower left")
    plt.savefig(path)
    plt.close()
    print(f"Saving file: {path}")
    path.chmod(0o777)

def main(
        input_dir,
        output_dir,
        year,
        n_trees,
        learning_rate,
        max_depth,
        min_child_weight):

    # Ensure input/output dirs are Path instances.
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    dfs = pd.DataFrame()

    ######### LOADING/PROCESSING THE DATAFRAME
    # Given the large number of files to be processed, the total df (all jets with inputs/labels) is available as 'all_jets.h5'
    # if it is not, the total df is generated and saved as 'all_jets.h5'
    if os.path.exists(output_dir / 'all_jets.h5')==0:
        # Glob the files, and convert to a sorted list.
        # sgn and bkg files in separated list so it is easy to filter them when needed
        fnames_sgn = sorted(input_dir.glob(f"sgn_{year}_hotvr_h5/*/*.h5"))
        fnames_bkg = sorted(input_dir.glob(f"bkg_{year}_hotvr_h5/qcd*/*.h5"))  # preliminary study with only QCD files 
        fnames =  fnames_sgn + fnames_bkg
        print(f"Found {len(fnames)} files in {input_dir}")
        # files are read and processed in chunk and save in temp files that are then deleted
        loading_and_processing_h5(fnames, output_dir)
        
    # reading "all_jets.h5" if exists
    dfs = pd.read_hdf(output_dir / 'all_jets.h5')
    if dfs.empty: 
        print("Unable to find 'all_jets.h5' or unable to process the individual files")
        sys.exit()
    
    print('### Successfully loaded the h5 with all the jets')
    # Shuffle the dataframe
    dfs = dfs.sample(frac=1).reset_index(drop=True)
    print(dfs)
    #########

    output_dir = output_dir / 'trainBDT_outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir.chmod(0o777)

    ######### TRAINING-TEST SPLITTING
    # create the training set if does not exist
    print('\n### Train-test splitting')
    features_array, label_array = [], []

    if os.path.exists(output_dir / 'trainBDT_outputs' / 'trainData.npz')==0 or os.path.exists(output_dir / 'trainBDT_outputs' / 'trainData_resampled.npz')==0:    
        print('### Creating training dataset...')
        if resampling:
            print('### RE-SAMPLING')
            # jets resampling by pT
            ptSig = dfs[dfs['label']==1]['pt'].values
            ptBkg = dfs[dfs['label']==0]['pt'].values
            print(f"Pre-resampling [SGN]: {ptSig.shape[0]}")
            print(f"Pre-resampling [BKG]: {ptBkg.shape[0]}")

            pT_plot(ptBkg, ptSig, output_dir / f"pT_distribution_SGNvsBKG_pre_sampling.png")

            # separating bkg features from sgn features
            bkg_data = dfs[dfs['label'] == 0][[input for input in inputs]].values
            sgn_data = dfs[dfs['label'] == 1][[input for input in inputs]].values

            # the target distribution is the sgn pT
            # additional info on the multi_dim_resample in "resample.py"
            targetHist, _ = np.histogram(sgn_data[:, -1], bins=binning)
            bkg_resampled, sgn_resampled = multi_dim_resample([bkg_data, sgn_data], targetHist, binning)

            ptSig = sgn_resampled[:, 5]
            ptBkg = bkg_resampled[:, 5]

            print(f"Post-resampling [SGN]: {sgn_resampled.shape[0]}")
            print(f"Post-resampling [BKG]: {bkg_resampled.shape[0]}")

            pT_plot(ptBkg, ptSig, output_dir / f"pT_distribution_SGNvsBKG_post_sampling.png")

            features_array = np.vstack((bkg_resampled, sgn_resampled))
            label_array = np.concatenate((np.zeros(len(bkg_resampled)), np.ones(len(sgn_resampled))))

            creating_training_data((features_array, label_array), output_dir/f"trainData_resampled.npz")
        else:
            label_array = dfs['label'].values.reshape(-1, 1)
            features_array = dfs[[input for input in inputs if input != 'pt']].values # to avoid 'pt' in the dataframe

            creating_training_data((features_array, label_array), output_dir/f"trainData.npz")

    else: 
        fpath = ''
        if resampling: fpath = output_dir / 'trainData_resampled.npz'
        else: fpath = output_dir / 'trainData.npz'

        loaded = np.load(fpath)
        features_array = loaded['features']
        label_array = loaded['labels']

    trainFeatures,testFeatures,trainLabel,testLabel, = sklearn.model_selection.train_test_split(
        features_array, label_array, test_size=0.2, shuffle=True, random_state=24155
    )

    print(f"\nTotal samples {features_array.shape}")
    print(f"Split samples -> train: {trainFeatures.shape}, test: {testFeatures.shape}")

    # use test=train for testing
    # testFeatures = trainFeatures
    # testLabel = trainLabel

    # ensure a one-dimensional array (*, )
    trainLabel = trainLabel.ravel()
    testLabel = testLabel.ravel()

    # print(f"Train/test: {trainFeatures.shape}/{testFeatures.shape}")
    # for i, inputName in enumerate(inputs):
    #     print(inputName,trainFeatures[0:10,i])

    del features_array, label_array
    gc.collect()
    #########

    ######### Cut-based top tagger: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetTopTagging#Working_point_and_scale_factors
    cut_based_efficiencies = cut_based_tagger_ROC(testFeatures, testLabel)

    start_time = time.time()
    ######### XGBOOST 
    xgboost_params = {
    'n_estimators': n_trees, 'learning_rate': learning_rate, 'max_depth': max_depth, 'min_child_weight': min_child_weight
    }
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
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    output_dir = output_dir / "plots_and_stats"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir.chmod(0o777)
    fpath = output_dir / f"hadronicTopVSQCD_bdt_nTrees{xgboost_params['n_estimators']}_maxDepth{xgboost_params['max_depth']}_learningRate{xgboost_params['learning_rate']}_minChildWeight{xgboost_params['min_child_weight']}.bin"
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

    fpath = output_dir / f"hadronicTopVSQCD_stat_nTrees{xgboost_params['n_estimators']}_maxDepth{xgboost_params['max_depth']}_learningRate{xgboost_params['learning_rate']}_minChildWeight{xgboost_params['min_child_weight']}.csv"
    
    stats = {
        'train_last_logloss': [results['validation_0']['logloss'][-1]], 'test_last_logloss': [results['validation_1']['logloss'][-1]],
        'train_AUC': [aucTrain], 'test_AUC': [aucTest], 'train_accuracy': [100.*accTrain], 'test_accuracy': [100.*accTest],
        'execution_time': [elapsed_time]
    }
    df_stats = pd.DataFrame(stats)
    df_stats.to_csv(fpath, index=False)
    print(f"Saving stats file: {fpath}")

    fpath = output_dir / f"hadronicTopVSQCD_lr_nTrees{xgboost_params['n_estimators']}_maxDepth{xgboost_params['max_depth']}_learningRate{xgboost_params['learning_rate']}_minChildWeight{xgboost_params['min_child_weight']}.png"
    draw_loss(fpath, results, n_estimators=n_trees)

    fpath = output_dir / f"hadronicTopVSQCD_roc__nTrees{xgboost_params['n_estimators']}_maxDepth{xgboost_params['max_depth']}_learningRate{xgboost_params['learning_rate']}_minChildWeight{xgboost_params['min_child_weight']}.png"
    draw_roc(fpath, rates, cut_based_efficiencies=cut_based_efficiencies)

    fpath = output_dir / f"hadronicTopVSQCD_dist_nTrees{xgboost_params['n_estimators']}_maxDepth{xgboost_params['max_depth']}_learningRate{xgboost_params['learning_rate']}_minChildWeight{xgboost_params['min_child_weight']}.png"
    draw_score(fpath, predictions)



#################################################

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
    # xgboost parameters
    parser.add_argument('--n_trees', type=int, default=500,
        help='Number of estimators.')
    parser.add_argument('--learning_rate', type=float, default=.05,
        help='Learning Rate.')
    parser.add_argument('--max_depth', type=int, default=3,
        help='Learning Rate.')
    parser.add_argument('--min_child_weight', type=float, default=.01,
        help='Learning Rate.')

    args = parser.parse_args(argv)

    # If output directory is not provided, assume we want the output to be
    # alongside the input directory.
    if args.output_dir is None:
        args.output_dir = args.input_dir

    # Return the options as a dictionary.
    return vars(args)


if __name__ == "__main__":
    main(**parse_args())
