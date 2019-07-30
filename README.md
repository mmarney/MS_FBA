# Add_adducts
The Add_adducts function is responible for finding the compounds potential mass with adducts based on the ionization mode of the instrument.  The funstion takes in a vector of masses and an ion mode and outputs a dataframe of the compounds mass with the potential adducts as the columns.
# isotope_count 
The isotope_count function is responisble for predicting the isotope ratio of a compound based on the molecular formula.  The count of each element is multiplied by the natural isotope distrubution of that element's M+1 isotope.  These values are added and result in the expected isotope ratio.
# plot_mzxml_spectum
This is a quick way to plot the mass to charge value vs intentsity of a specific scan in a file.  The mzxml_file and the scan number are requied for this funstion.
# ranked_macthes_df
This is the function to decode the potential matches column of potential_matches data frame produced by MS_FBA.  It will take one potential match cell and output a ranked list of those matches in a much easier format to be read by people.  
# MS_FBA
The MS_FBA function evaluates XCMS Online results against a custom compound list. xcms_tsvfile,files,file_prefixes,compounds_tsv,ion_mode,pvalue = 0.005,ppm =5 There are eight variables that must be set for the MS_FBA function to be completed
1. xcms_tsvfile: This is located in the results folder downloaded after the XCMS Online run is completed the exact title will be XCMS.annotated.diffreport.tsv
2. files: This needs to be a character vector of the files names that were produced from your Mass Spec. These are the same files that you uploaded to XCMS Online. Easiest to have the files in the working directory so the path isn't needed every time.
3. file_prefixes: This is also a character vector, but it is just the prefixes for how you split up your sample classes.
4. compounds_tsv: This is the tsv of the expected compounds that you are wishing to compare the XCMS Online results against. This tsv must be formatted 'compound\tmonoisotopic mass\tmolecular formula\n' followed by the compounds. This file can be created from any Flux Balance analysis program, manually from any list of compounds, or directly in the PyFBA script that is provided.
5. ion_mode: This is either 'N' or 'P' for if the experiment was run in negative or positive mode. This variable changes the potential adducts filtered through
6. pvalue: This variable determines that features that should be considered as a viable feature from the XCMS Online results. If the Pvalue is lower than the Pvalue then the feature will be included in the search otherwise it will be disregarded.
7. ppm: Parts Per Million is used to create a mass matching window for features to peaks in the raw files and features to potential compounds.
8. rtmulti: The retention time multiplier is responsible for the retention time window created. Since XCMS Online preforms Peak Alignment, it is important to expand the retention time window to ensure that the top of the peak is always found for each feature in each file. The larger the window the longer the program takes to run but the more certain you will be that you found the correct peak top. The retention time window is created by taking the RT difference from XCMS Online results and multiplying it by the rtmulti then adding the resulting time window on either side of the pre-existing window. example ( if the RT window was from 1.9 to 2 minutes and the rtmulti was 2 the new RT window would be 1.7 to 2.2

