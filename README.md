# Micro_Interference4SLGW_identification
This repository is used to share the code and data for microlensing diffraction and SL identification

### Step 0: Plot_result.py in Plot_code folder can plot all of the figures in this paper.
One can run this code to reproduce the result quickly.



### Step 1: Run code in the `SLGW_simulation` folder to generate unlensed gravitational waves and only macro-lensed gravitational wave data.
1. **Run `GWDistribution.py`**: Generate the parameters for gravitational wave events.
2. **Run `Generate_GW_file_FD.py`**: Create the gravitational waveforms, save them to files, and then perform template matching to calculate the signal-to-noise ratio (SNR).
3. **Run `ReadAndEstimate.py`**: Perform parameter estimation on the saved gravitational wave files.
4. **Run `LensParameterGen.py`**: Generate the lens parameters. The `sersic_profile.py` file is called here to set the microlensing density based on the Sersic profile.
5. **Run `Generate_GW_file_Lensed_4_macro_4_SNR.py`**: Compute the SNR considering only macrolensing amplification, and save the waveform data into files.

---

### Step 2: Compute microlensing diffraction integrals in the `Microlensing_Simulation` folder.
1. **Run `Totalmicro.cpp`**: Calculate the microlensing diffraction integrals.
(The microlensing data can be obtained from https://pan.bnu.edu.cn/l/X1QPKG.)
---

### Step 3: Return to the `SLGW_simulation` folder and generate the strongly lensed gravitational waveforms affected by microlensing.
1. **Run `Generate_GW_file_Lensed_4_micro.py`**: Generate the waveforms considering microlensing effects.
2. **Run `ReadAndEstimate_with_micro*.py`**: Perform parameter estimation for the microlensed images.

---

### Step 4: Run cWB-related code in the `cWB_core` folder and identify pairs.
1. **Run `Plot_SNR_VS_match`**: Generate the theoretical match for all lensing cases.
2. **Run `Gen_input_cat_file_unlens`**: Generate FRAMES data, input files, and DQ files for the unlensed cases, used as input for cWB.
3. **Run `Gen_input_cat_file_micro_Total`**: Generate FRAMES data, input files, and DQ files for all microlensed cases.
4. **Run `cal_match_unlens`**: Compute the match for unlensed cases.
5. **Run `cal_match_micro_Total` and `cal_match_micro_Total_final`**: Compute the match for microlensed cases.
6. **Run `OverlapAna.py`**: Analyze the Bayes factor for pairwise comparisons of unlensed cases.
7. **Run `Overlap_lens.py`**: Analyze the Bayes factor for pairwise comparisons of lensed cases.
8. **Run `Plot_Overlap_and_find_high_bayes_event.py`**: Plot the Bayes factors for unlensed and lensed cases, and identify events with high Bayes factors.

---

### Step 5: Identify the host galaxy in the `Host_identification` folder.
1. **Run `Lens_Galaxy_Sample.py`**: Generate galaxy-galaxy strong lensing systems based on the JWST star catalog.
2. **Run `Lens_light_population.pu`**: Calculate the effective radius (R_e) and magnitude of the lens galaxy using the fundamental plane.
3. **Run `Source_light_population_improve.py`**: Identify the magnitudes, effective radius (R_e), and Sersic index (R_sersic_n) of the host galaxies for the three selected SLGW events from the JWST star catalog, and save cases with star formation rates greater than 1.
4. **Run `Dec_Rec_Full_CSST/JWST_host_pop.py`**: Simulate the host galaxy-galaxy strong lensing system and reconstruct the lensing parameters.
5. **Run `TD_consistent_JWST_host_pop.py`**: Calculate the time delay and sample the four-image region based on the star formation rate.
6. **Run `Area_in_quadruple_region_host (unhost).py`**: Calculate the area of the quadruple region for all host (and unhost) galaxies.
7. **Run `Bayes_factor4TD.py`**: Calculate the Bayes factors for the time delay in unlensed and host lens galaxies.
