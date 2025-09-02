## 20250824

- Exome history data obtained
- From exome history, used a [keywords matching algorithm](https://github.com/veritasvita3/CoreGenomicsProject/tree/main/Scripts/keyword_match.py) to filter renal-cardio complication within patients with the same disease
- Used Biobert to further filter using machine learning as given in the [script](Scriptshttps://github.com/veritasvita3/CoreGenomicsProject/tree/main/Scripts/mlBioBert_embeddings.py)

## 20250902

- Filtered TCGA database with 3d neurosphere specimen, open access, CNV data
- Installed gdc-client on wsl
- Downloaded CNV files via gdc-client
- Started inspecting aberrant loci for genes 
- Installed biomaRt on R-studio
- Using the ensembl names of genes from particular loci, genes are [annotated](https://github.com/veritasvita3/CoreGenomicsProject/tree/main/Files/annotated_gene_list_1.xlsx)