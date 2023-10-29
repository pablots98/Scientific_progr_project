# Main scripts

In this folder we find the files with the main codes for this project, as well as PDFs with the code and better comments for each algorithm


Heatmap: This code generates the heatmap based on all the samples and 123 genes that are important for senescence. It had to be done in R due to problems running it in MATLAB. 


Global_threshold: Works with the global thresholding technique, the simplest, establishing a single threshold per sample for all genes. You can also find a PDF file with the livescript and better code comments.


LocalT2_threshold: This code works with the localT2 technique, establishing a threshold per gene. You can also find a PDF file with the livescript and better code comments.


Localgini: This code works with the localT2 technique, establishing a threshold per gene based on the gini coefficient. You can also find a PDF file with the livescript and better code comments.


StanDep: This code returns an array with what are considered core reactions. A percentage of accuracy could not be achieved since it was not possible to find an effective way to go from the list of housekeeping genes to creating one of housekeeping reactions to be able to compare. Still, it is expected to achieve this in the future.
