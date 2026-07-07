#June 2026

from tqdm import tqdm
import keras
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from tensorflow.keras.callbacks import EarlyStopping
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
import seaborn as sns
import math
%matplotlib inline
%config InlineBackend.figure_format = 'retina'

num = 30 # add on if thats the requirement. 5 is the minimum.
bin_num = 3.5
path = 'C:/Users/Asus/Final_flvl_csv/Extracted_data/'
filename = 'QQ_standard_flvl_nonan_data_{}_bin.csv'.format(bin_num)
pd.set_option('display.max_columns', None)
df_standard_flvl_nonan = pd.read_csv(path+filename)

all_predictions = []; mean_predictions_list = []; all_histories = []; all_peaks = [];
peaks_extrapolated = []; peaks_within = [];
a_list, b_list, c_list = [], [], [];
a_mean_list, b_mean_list, c_mean_list = [], [], [];
z_test_scaled_list = []; z_peak_mean_list = []; z_peak_real_mean_list = [];
all_peaks_argmax = []; scaler_list = [];
all_predictions_test=[]; mean_predictions_list_test =[];
peaks_within_test = []; all_peaks_argmax_test=[]; all_peaks_test=[];
all_peaks_mean = []; all_peaks_argmax_mean = []; z_test_list = [];
def reg_r2(x, a, b, c):
    return a*x**2 + b*x + c


#### ---- DNN model ---- ####
def build_model():
    model_loop = keras.models.Sequential([
        #changed from relu to elu
        keras.layers.Dense(64, activation='elu', use_bias=False),
        keras.layers.BatchNormalization(),
        keras.layers.Dropout(0.2),
        keras.layers.Dense(32, activation='elu', use_bias=False),
        keras.layers.BatchNormalization(),
        keras.layers.Dense(16, activation='elu', use_bias=False),
        keras.layers.BatchNormalization(),
        keras.layers.Dense(1)])
    return model_loop


#### ---- data ---- ####
combined_sorted = df_standard_flvl_nonan.copy()
z = combined_sorted['z'].values.reshape(-1,1)
eps = np.log10(combined_sorted['epsopt']).values.reshape(-1,1)
z_data_min = float(z.min())
z_data_max = float(z.max())
#It can be changed from z.min to zero;
z_fixed = np.linspace(z_data_min, z_data_max, 500).reshape(-1,1)

#### ---- loop ---- ####
for run in tqdm(range(num)):
    while True:
        #for the "stratify" issue.
        try:
            z_bins = pd.cut(z.squeeze(), bins=10, labels=False)
            z_train0, z_test, eps_train0, eps_test = train_test_split(
                z, eps, test_size=0.25, stratify=z_bins)
        except ValueError:
            z_train0, z_test, eps_train0, eps_test = train_test_split(
                z, eps, test_size=0.25)
        try:
            z_bins_train = pd.cut(z_train0.squeeze(), bins=10, labels=False)
            z_train1, z_validation, eps_train1, eps_validation = train_test_split(
                z_train0, eps_train0, test_size=0.25, stratify=z_bins_train)
        except ValueError:
            z_train1, z_validation, eps_train1, eps_validation = train_test_split(
                z_train0, eps_train0, test_size=0.25)

        scaler = StandardScaler()
        z_train_scaled = scaler.fit_transform(z_train1)
        z_validation_scaled = scaler.transform(z_validation)
        z_test_scaled = scaler.transform(z_test)
        z_fixed_scaled = scaler.transform(z_fixed)
        early_stop = EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True)

        model_loop = build_model()
        #didnt add any constraints for the loss so the model can be data-trained as natural as possible.
        model_loop.compile(optimizer=keras.optimizers.Adam(1e-3), loss=tf.keras.losses.Huber(delta=0.12)) #was mse before - can be changed.
        history_loop = model_loop.fit(z_train_scaled, eps_train1,
                       validation_data=(z_validation_scaled, eps_validation),
                       epochs=200, batch_size=128, callbacks=[early_stop], verbose=0)

        preds_test = model_loop.predict(z_test_scaled, verbose=0) #put this here for the plot
        preds_fixed = model_loop.predict(z_fixed_scaled, verbose=0)

        try:
            parameters, _ = curve_fit(reg_r2, z_fixed_scaled.reshape(-1), preds_fixed.squeeze())
            a, b, c = parameters
            if a < 0:
                z_peak_scaled_cf = -b / (2*a)
                z_min_scaled = z_fixed_scaled.min().squeeze();z_max_scaled = z_fixed_scaled.max().squeeze()
                z_peak_real_cf = scaler.inverse_transform([[z_peak_scaled_cf]])[0][0]
                #for the z > 10 being after z > 0 in the previous sessions
                if z_peak_real_cf > 10:
                    print(f"Run {run}: non-physical peak [z = {round(z_peak_real_cf,3)}], retrying")
                    continue
                elif z_peak_real_cf <= 0:
                    print(f"Run {run}: negative redshift [z = {round(z_peak_real_cf,3)}], retrying")
                    continue
                #if 0 < z < 10
                else:
                    if z_min_scaled <= z_peak_scaled_cf <= z_max_scaled:
                        peaks_within.append(z_peak_real_cf)
                        print(f"Run {run} [curve_fit]: peak at z = {z_peak_real_cf:.3f}")
                                #argmax peak
                        y_smoothed = gaussian_filter1d(preds_fixed.squeeze(), sigma=10)
                        ipk = int(np.argmax(y_smoothed))
                        z_peak_argmax = float(z_fixed.squeeze()[ipk])
                        all_peaks_argmax.append(z_peak_argmax)
                        print(f"Run {run} [argmax]:    peak at z = {z_peak_argmax:.3f}")
                                #curvefit lists
                        all_peaks.append(z_peak_real_cf)
                        all_predictions.append(preds_fixed.squeeze())
                        all_predictions_test.append(preds_test.squeeze())
                        z_test_list.append(z_test)
                        mean_predictions_list.append(np.mean(all_predictions, axis=0))
                        mean_predictions_list_test.append(np.mean(all_predictions_test, axis=0))
                        all_histories.append(history_loop)
                        a_list.append(a); b_list.append(b); c_list.append(c)
                        # z_test_scaled_list.append(z_test_scaled)
                        scaler_list.append(scaler)
                    else:
                        peaks_extrapolated.append(z_peak_real_cf)
                        print(f"Run {run} [curve_fit]: peak at z = {z_peak_real_cf:.3f} (extrapolated)")
                                #argmax peak
                        y_smoothed = gaussian_filter1d(preds_fixed.squeeze(), sigma=10)
                        ipk = int(np.argmax(y_smoothed))
                        z_peak_argmax = float(z_fixed.squeeze()[ipk])
                        all_peaks_argmax.append(z_peak_argmax)
                        print(f"Run {run} [argmax]:    peak at z = {z_peak_argmax:.3f} (extrapolated)")
                                #curvefit lists
                        all_peaks.append(z_peak_real_cf)
                        all_predictions.append(preds_fixed.squeeze())
                        all_predictions_test.append(preds_test.squeeze())
                        z_test_list.append(z_test)
                        mean_predictions_list.append(np.mean(all_predictions, axis=0))
                        mean_predictions_list_test.append(np.mean(all_predictions_test, axis=0))
                        all_histories.append(history_loop)
                        a_list.append(a); b_list.append(b); c_list.append(c)
                        # z_test_scaled_list.append(z_test_scaled)
                        scaler_list.append(scaler)
                break
            else:
                print(f"Run {run}: parabola opens upward [a={round(a,3)}], retrying")
                continue
        except Exception as e:
            print(f"Run {run}: curve_fit failed: {e}, retrying")
            continue

mean_predictions = np.mean(all_predictions, axis=0) #total mean
mean_predictions_test = np.mean(all_predictions_test,axis=0)
std_predictions = np.std(all_predictions, axis=0)

#### ---- evolution of the mean peak per run ---- ####
for i in range(len(mean_predictions_list)):
    z_fixed_scaled_i = scaler_list[i].transform(z_fixed)
    # mean_fixed_i = np.mean(all_predictions[i], axis=0)
    try:
        parameters_mean, _ = curve_fit(reg_r2, z_fixed_scaled_i.reshape(-1),
                                        mean_predictions_list[i])
        a_mean, b_mean, c_mean = parameters_mean
        if a_mean < 0:
            z_peak_mean = -b_mean / (2*a_mean)
            z_peak_real_mean = scaler_list[i].inverse_transform([[z_peak_mean]])[0][0]
            z_min_i = z_fixed_scaled_i.min()
            z_max_i = z_fixed_scaled_i.max()

            if z_min_i <= z_peak_mean <= z_max_i:
                print(f"Mean Run +{i}: [curvefit] peak at z = {z_peak_real_mean:.3f}")
            else:
                print(f"Mean Run +{i}: [curvefit] peak at z = {z_peak_real_mean:.3f} (extrapolated)")

            if z_peak_real_mean > 10:
                print(f"Mean Run +{i}: non-physical, skipping")
            elif z_peak_real_mean <= 0:
                print(f"Mean Run +{i}: negative redshift, skipping")
            else:
                y_smoothed_mean = gaussian_filter1d(mean_predictions_list[i].squeeze(), sigma=10)
                ipk_mean = int(np.argmax(y_smoothed_mean))
                z_peak_argmax_mean = float(z_fixed.squeeze()[ipk_mean])
                all_peaks_argmax_mean.append(z_peak_argmax_mean)
                print(f"\n Mean Run +{i}: [argmax] peak at z = {z_peak_argmax_mean:.3f}")

                all_peaks_mean.append(z_peak_real_mean)
                a_mean_list.append(a_mean); b_mean_list.append(b_mean); c_mean_list.append(c_mean)
        else:
            print(f"Mean Run +{i}: parabola opens upward")
    except Exception as e:
        print(f"Mean Run +{i}: curve_fit failed: {e}")

#### ---- argmax peaks for mean ---- ####
P = np.stack(all_predictions, axis=0)
y_med_grid = np.median(P, axis=0)
std_grid = np.std(P, axis=0)
y_smoothed_final = gaussian_filter1d(y_med_grid, sigma=10)
ipk_final = int(np.argmax(y_smoothed_final))
z_peak_argmax_final = float(z_fixed.squeeze()[ipk_final])
print(f"\nFinal mean [argmax]: peak at z = {z_peak_argmax_final:.3f}")


### ---- for QQ only ---- ###
sort_idx = np.argsort(z.squeeze())
z_sorted = z.squeeze()[sort_idx]
eps_sorted = eps.squeeze()[sort_idx]
qq_y_smoothed_final = gaussian_filter1d(eps_sorted, sigma=10)
qq_ipk_final = int(np.argmax(qq_y_smoothed_final))
qq_z_peak_argmax_final = float(z_sorted[qq_ipk_final])
print(f"\nQQ [argmax]: peak at z = {qq_z_peak_argmax_final:.3f}")


#### ---- peaks info ---- ####
all_peaks = np.array(all_peaks)
all_peaks_argmax = np.array(all_peaks_argmax)
print(f"\nTotal runs with a peak: {len(all_peaks)}")
print(f"\n###---- curve_fit peaks ----###")
print(f"Mean peak: {np.mean(all_peaks):.3f} ± {np.std(all_peaks):.3f}")
print(f"Final Peak: {np.mean(all_peaks):.3f} ± {np.std(all_peaks):.3f} "
      f"~ ({np.mean(all_peaks)-np.std(all_peaks):.3f}, "
      f"{np.mean(all_peaks)+np.std(all_peaks):.3f})")
print(f"Peaks z=1-3.0: {np.sum((all_peaks>=1)&(all_peaks<=3.0555))} out of {len(all_peaks)}")
print(f"Peaks z=0.5-4.0: {np.sum((all_peaks>=0.5555)&(all_peaks<=4.0555))} out of {len(all_peaks)}")
print(f"\n###---- argmax peaks ----###")
print(f"Total: {len(all_peaks_argmax)}")
print(f"Mean peak: {np.mean(all_peaks_argmax):.3f} ± {np.std(all_peaks_argmax):.3f}")
print(f"Peaks z=1-3.5: {np.sum((all_peaks_argmax>=1)&(all_peaks_argmax<=3.5))} out of {len(all_peaks_argmax)}")

#### ---- plotting histogram for peaks ---- ####
fig, axes = plt.subplots(1, 2, figsize=(15, 3))
axes[0].hist(all_peaks, bins=15, edgecolor='black')
axes[0].axvline(np.mean(all_peaks)+np.std(all_peaks), color='orange',
                linestyle='--', label=f'Mean + STD = {np.mean(all_peaks)+np.std(all_peaks):.3f}')
axes[0].axvline(np.mean(all_peaks), color='red', linestyle='--',
                label=f'Mean = {np.mean(all_peaks):.3f}')
axes[0].axvline(np.mean(all_peaks)-np.std(all_peaks), color='orange',
                linestyle='--', label=f'Mean - STD = {np.mean(all_peaks)-np.std(all_peaks):.3f}')
axes[0].set_xlabel('Peak redshift')
axes[0].set_ylabel('Count')
axes[0].set_title('Distribution of peak locations with curve_fit')
axes[0].legend()

axes[1].hist(all_peaks_argmax, bins=15, edgecolor='black')
axes[1].axvline(np.mean(all_peaks_argmax)+np.std(all_peaks_argmax), color='orange',
                linestyle='--', label=f'+ STD = {np.mean(all_peaks_argmax)+np.std(all_peaks_argmax):.3f}')
axes[1].axvline(np.mean(all_peaks_argmax), color='red', linestyle='--',
                label=f'Mean = {np.mean(all_peaks_argmax):.3f}')
axes[1].axvline(np.mean(all_peaks_argmax)-np.std(all_peaks_argmax), color='orange',
                linestyle='--', label=f'- STD = {np.mean(all_peaks_argmax)-np.std(all_peaks_argmax):.3f}')
axes[1].set_xlabel('Peak redshift')
axes[1].set_ylabel('Count')
axes[1].set_title('Distribution of peak locations with argmax')
axes[1].legend()
plt.tight_layout()
sns.set_style("whitegrid")
plt.show()

#### ---- loss subplots ---- ####
n_runs = len(all_histories); n_cols = 5; n_rows = math.ceil(n_runs / n_cols)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, n_rows*3))
axes = axes.flatten()
for i, history in enumerate(all_histories):
    axes[i].plot(all_histories[i].history["loss"], label="train loss")
    axes[i].plot(all_histories[i].history["val_loss"], label="validation loss")
    axes[i].set_xlabel('Epoch'); axes[i].set_ylabel('Loss')
    axes[i].set_title(f'Run {i+1}')
    axes[0].legend()
for j in range(i+1, len(axes)):
    axes[j].set_visible(False)
plt.tight_layout()
sns.set_style("whitegrid")
plt.show()

#### ---- mean loss plot ---- ####
max_len = max(len(h.history['val_loss']) for h in all_histories)
def pad_history(hist, max_len):
    return list(hist) + [hist[-1]] * (max_len - len(hist))

mean_val_loss = np.mean([pad_history(h.history['val_loss'], max_len) for h in all_histories], axis=0)
mean_train_loss = np.mean([pad_history(h.history['loss'], max_len) for h in all_histories], axis=0)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(mean_train_loss, label="Mean train loss")
ax.plot(mean_val_loss, label="Mean validation loss")
ax.set_xlabel('Epoch'); ax.set_ylabel('Loss')
ax.legend()
plt.tight_layout()
sns.set_style("whitegrid")
plt.show()

#### ---- curve and argmax plot ---- ####
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(z.squeeze(), eps.squeeze(), c='b', alpha=0.5, s=50,
           label='QQ observation points')
plt.fill_between(z_fixed.squeeze(), min(y_med_grid - 3*std_grid), max(y_med_grid + 3*std_grid), alpha=0.1, 
                 color='orange', label=r' min/max $\pm 3\sigma$')
plt.fill_between(z_fixed.squeeze(), min(y_med_grid - 2*std_grid), max(y_med_grid + 2*std_grid), alpha=0.1, 
                 color='blue', label=r' $ min/max \pm 2\sigma$')
plt.fill_between(z_fixed.squeeze(), min(y_med_grid - std_grid), max(y_med_grid + std_grid), alpha=0.1, 
                 color='green', label=r'$ min/max \pm \sigma$')
ax.fill_between(z_fixed.squeeze(), y_med_grid - 3*std_grid, y_med_grid + 3*std_grid, alpha=0.3, 
                color='orange', label=r' $\pm 3\sigma$')
ax.fill_between(z_fixed.squeeze(), y_med_grid - 2*std_grid, y_med_grid + 2*std_grid, alpha=0.3, 
                color='blue', label=r' $\pm 2\sigma$')
ax.fill_between(z_fixed.squeeze(), y_med_grid - std_grid, y_med_grid + std_grid, alpha=0.3, 
                color='green', label=r'$\pm \sigma$')

ax.scatter(z_fixed.squeeze(), mean_predictions_list[0], alpha=0.3, c='red', s=30, label='Mean predictions (fixed grid for [0])')
for i in range(len(all_predictions_test)):
    if i == 0:
        ax.scatter(z_test_list[i].squeeze(), all_predictions_test[i], alpha=0.1, c='cyan', s=40, label='Predictions [test] across all runs')
    else:
        ax.scatter(z_test_list[i].squeeze(), all_predictions_test[i], alpha=0.1, c='cyan', s=40)
ax.plot(z_fixed.squeeze(), y_smoothed_final, color='purple', linewidth=2, label=f'Ensemble median [argmax]')
ax.axvline(qq_z_peak_argmax_final, color='pink', linestyle='--', label=f'QQ Peak = {qq_z_peak_argmax_final:.3f} [argmax]')
if len(a_mean_list) > 0:
    ax.plot(z_fixed.squeeze(), 
            reg_r2(scaler_list[0].transform(z_fixed).squeeze(), a_mean_list[0], b_mean_list[0], c_mean_list[0]), 
            c='black', alpha=0.5, label='curve_fit mean[0] predict')
ax.axvline(np.mean(all_peaks)+np.std(all_peaks), color='orange',
           linestyle='--', label=f'[curvefit] Mean + STD = {np.mean(all_peaks)+np.std(all_peaks):.3f}')
ax.axvline(np.mean(all_peaks), color='black', linestyle='--',
           label=f'[curvefit] Mean Peak = {np.mean(all_peaks):.3f}')
ax.axvline(np.mean(all_peaks)-np.std(all_peaks), color='orange',
           linestyle='--', label=f'[curvefit] Mean - STD = {np.mean(all_peaks)-np.std(all_peaks):.3f}')
ax.axvline(z_peak_argmax_final, color='magenta', linestyle='--',
           label=f'[argmax] Mean Peak = {z_peak_argmax_final:.3f}')
ax.axvspan(0, z_data_min, alpha=0.06, color='gray', label='Extrapolation region')
ax.axvline(z_data_min, color='gray', linestyle=':', linewidth=1.5)
ax.set_xlabel('$z$', size=14)
ax.set_ylabel(r'Log $\epsilon$', size=14)
ax.legend(fontsize=10, ncol=2, loc='best')
plt.tight_layout()
sns.set_style("whitegrid")
plt.show()