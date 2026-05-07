About:
This Shiny app calculates qPCR Fold Change and Copy Number.

How to Use:
Prepare your data in CSV format (see example below).
Upload the CSV file using the file input on the left.
Select your Housekeeping gene and Reference (Control) group.
View results in the main panel, or download as Excel/CSV files.


Access the app: 
https://janan91.shinyapps.io/Shiny_app_CT_To_CopyNumber_FoldChange/

Input File Format:
The first two columns must be: Replicate and Treatment.
Each additional column must be a gene (header = gene name; values = Ct values).
Example:
Replicate,Treatment,GeneA,GeneB,GAPDH
R1,Control,24.1,25.3,21.2
R2,Control,23.8,25.1,21.0
R1,Treated,22.7,24.2,20.5
R2,Treated,22.5,24.1,20.3
      
Steps Performed:
HK Stability: Calculates SD of Ct for each gene.
ΔCt: For each target gene, ΔCt = Ct(Housekeeping) - Ct(Target).
Copy Number: 2^ΔCt * 1000 for each gene/sample.
Fold Change: For each gene/sample, Copy Number / Average Copy Number of Reference group.
