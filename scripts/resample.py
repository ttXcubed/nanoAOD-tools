import numpy as np
import matplotlib.pyplot as plt


def multi_dim_resample(samples, targetHist, binning=np.logspace(0,4,100), oversample=1., threshold=2):
    all_resampled = []
    for data in samples:
        pt_values = data[:, -1]
        other_data = data[:, :-1]
        
        hists = []
        weightList = []
        
        hist, _ = np.histogram(pt_values, bins=binning)
        hists.append(hist)

        weight = np.zeros(len(binning) - 1)
        for i in range(len(hist)):
            ntarget = targetHist[i]
            if ntarget < threshold:
                continue
            nsample = hist[i]
            if nsample < threshold:
                continue

            # Calculate resample weight
            weight[i] = 1. * ntarget / nsample
            
        weightList.append(weight)

        # Rescale weights
        maxweight = max(weightList[0])
        weightList[0] *= oversample / maxweight

        resampled_data = []
        resampled_pt_values = []
        
        for idx, pt in enumerate(pt_values):
            ibin = np.digitize(pt, binning)
            if ibin < 0 or ibin >= len(binning):
                continue

            # Oversample in case weights > 1
            for _ in range(int(weightList[0][ibin - 1])):
                resampled_data.append(data[idx, :-1])
                resampled_pt_values.append(pt)

            # Randomly resample depending on weight
            if np.random.uniform(0, 1) > weightList[0][ibin - 1] - int(weightList[0][ibin - 1]):
                continue

            resampled_data.append(data[idx, :-1])
            resampled_pt_values.append(pt)

        # Convert lists back to arrays
        resampled_data = np.array(resampled_data)
        resampled_pt_values = np.array(resampled_pt_values).reshape(-1, 1)

        # Combine the data and pt values
        combined_resampled = np.hstack((resampled_data, resampled_pt_values))
        
        all_resampled.append(combined_resampled)

    return all_resampled
