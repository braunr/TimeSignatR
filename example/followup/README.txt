Data and scripts in this directory are designed to reproduce the plots referenced in the PNAS letter Braun &al (2019) "Accurate prediction of circadian time across platforms".

The purpose of this analysis was fourfold:

1) To demonstrate, in subjects with available melatonin data, that local midnight approximates DLMO25% in these highly selected subjects, making local time a proxy for circadian phase (hrs since DLMO25) [plots/letterFig1.pdf];

2) To demonstrate, in subjects with available melatonin data, that the output produced by TimeSignature as originally trained and published (Braun 2018) are an accurate prediction of circadian phase (hrs since DLMO25), as would be expected due to (1) above [plots/letterFig2.pdf, "oTS"];

3) To demonstrate, in subjects with available melatonin data, that when trained using circadian phase (rather than local time) as the outcome, TimeSignature produces highly accurate predictions of circadian phase across data from multiple studies and platforms [plots/letterFig2.pdf, top & "mTS"]; 

4) To demonstrate that TimeSignature produces more accurate and more stable (across studies) predicitions of circadian phase than both the PLSR model and the two-timepoint "differential PLSR" model [plots/letterFig2.pdf, bottom].

This directory contains:

* TS-DLMO-PLSR-comparison.R = the script to produce the plots
* matlab_PLSR_files/ = directory containing input/output data files for use with Laing &al's matlab PLSR code
* plots/ = directory containing the output plots for the letter, along with captions


