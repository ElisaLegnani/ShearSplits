# ShearSplits

Constructs physically motivated blue/red galaxy subsamples for cosmic shear analyses, aiming to isolate IA contamination.

Splits DES Y3 redshift bins into blue/red subsamples, using Balrog injections and Bagpipes SED fits to derive physical properties (stellar mass, sSFR). Uses a DecisionTree / RandomForest classifiers with leaf reassignment to optimize the split for maximal SNR in detecting distinct IA signals across subsamples.

Slides with progress [here](https://docs.google.com/presentation/d/1c1TtCute3TPkhdmHo9TaczSFcOCW1qLBMDqlfJuds6w/edit?usp=sharing).


Full Balrog-Bagpipes redshift and physical properties distributions

<img class="post-img" height=200 align="left" src="/figures/nzs_hist_full.png"/>
<img class="post-img" height=200 align="left" src="/figures/stell_ssfr_hist_full.png"/>
