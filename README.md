# Analysis code

Regional cortical brain age is independent of dopaminergic denervation in Parkinson's disease
Hu X, Zhu S, Qian X, Li W.  Brain Communications.

Code is released under the MIT licence. It reproduces every analysis and figure reported in
the paper. Manuscript-preparation utilities are not included; they are not needed to
reproduce the results.


> Note: some figure scripts read previously generated figure files from a directory
> referred to as `figures_previous`. Set it to wherever your earlier figure outputs live,
> or regenerate those figures from the analysis scripts first.

## Contents

### Cortical processing and feature extraction
- `build_features.py`
- `build_features_v2.py`
- `check_hc_done.sh`
- `check_json_site.py`
- `check_ordering.py`
- `check_thickness.sh`
- `feature_reconstruct.py`
- `fix_dns.sh`
- `gen_done_list.sh`
- `gen_final_lists.sh`
- `gen_hc_idlist.sh`
- `gen_idlist.sh`
- `gen_pd_p2.sh`
- `inspect_civet_out.sh`
- `install_apptainer.sh`
- `inventory_newhc.sh`
- `launch_civet.sh`
- `launch_hc.sh`
- `launch_newhc.sh`
- `launch_pd_all.sh`
- `probe_civet_container.sh`
- `probe_civet_opts.sh`
- `progress_newhc.sh`
- `progress_pd.sh`
- `proxy_test.sh`
- `quickstat.sh`
- `run_civet.sh`
- `run_civet_batch.sh`
- `run_civet_test.sh`
- `scan_ledd.py`
- `setup_ext4_scratch.sh`
- `status_civet.sh`
- `verify_batch.sh`

### Regional brain-age estimation, age-bias correction and harmonization
- `build_clinical_table.py`
- `build_master_clinical.py`
- `build_master_clinical_extended.py`
- `compare_age.py`
- `dat_brainage_analysis.py`
- `eval_hc.py`
- `explore_clinical.py`
- `extract_scanner.py`
- `final_analysis.py`
- `harmonize_analyze.py`
- `lme_scanner.py`

### Primary analyses
- `analyze_pd.py`
- `dat_brainage_analysis.py`
- `dat_variance_partition.py`
- `final_analysis.py`
- `incremental_validity.py`
- `pd_subtype_brainage.py`
- `pd_subtype_v2.py`
- `reviewer_sensitivity.py`

### Extended analyses
- `build_master_clinical_extended.py`
- `revision_A2_scale_robustness.py`
- `revision_A_ready.py`
- `revision_B2_gwr.py`
- `revision_B2b_incremental.py`
- `revision_B2c_gwr_site.py`
- `revision_B5B8_outcomes.py`
- `revision_B7_copathology.py`
- `revision_longitudinal.py`

### Figures
- `brain_surface.py`
- `build_brainnet_volume.py`
- `build_revision_figures.py`
- `build_revision_supplement.py`
- `compose_brainfig.py`
- `compose_fig3.py`
- `make_S2_S3.py`
- `make_S5_S6.py`
- `make_fig5_dat.py`
- `make_figures.py`
- `make_figures_final.py`
- `make_figures_nature.py`
- `make_flow_diagram.py`
- `make_graphical_abstract_R1.py`
- `make_remaining_figs.py`

## Reproducing brain-PAD from model output

The network-to-atlas mapping is defined in `revision_B2_gwr.py` (dictionary `NET`, AAL indices).
The age-bias correction is a linear fit of predicted on chronological age estimated in the
261 PPMI healthy controls; brain-PAD is the residual. See `revision_A_ready.py`.

## Figure pipeline order

Later scripts refine the output of earlier ones, so run them in this order:

1. `build_revision_figures.py`      main Figures 2-5
2. `fix_revision_fig1.py`           Figure 1
3. `build_revision_supplement.py`   Supplementary Figure 6 and all supplementary tables
4. `fix_revision_supp_figs.py`      Supplementary Figures 1 and 7
5. `fix_revision_supp_figs2.py`     Supplementary Figures 5 and 7
6. `make_graphical_abstract_R1.py`  graphical abstract

## Data

Individual-level clinical and imaging data are available from the Parkinson's Progression
Markers Initiative (www.ppmi-info.org) under its data use agreement and cannot be
redistributed here. No participant-level data are contained in this repository. The scripts
expect PPMI files in the layout described in each header.

## Model weights

The regional brain-age model weights are those distributed by the original authors
(regional_Brain_age) and were used without modification or retraining.
