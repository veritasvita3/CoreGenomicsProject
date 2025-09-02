import pandas as pd
import re

df = pd.read_excel(r"C:\Users\verit\OneDrive\CoreGenomicsPrivate\exome_clinical_history.xlsx")


for col in ["Clinical History", "FAMILY HISTORY"]:
    df[col] = (
        df[col].astype(str)
        .str.lower()
        .str.replace(r"[^a-z0-9\s]", " ", regex=True)   
        .str.replace(r"\s+", " ", regex=True)           
        .str.strip()
    )


renal_terms = [
    "renal", "kidney", "nephro", "nephritis", "nephropathy",
    "ckd", "esrd", "aki", "dialysis"
]

cardiac_terms = [
    "cardiac", "heart", "coronary", "angina", "myocardial",
    "mi", "infarction", "chf", "cardiomyopathy", "arrhythmia",
    "atrial fibrillation"
]


renal_pattern = re.compile(r"\b(" + "|".join(renal_terms) + r")\b", re.IGNORECASE)
cardiac_pattern = re.compile(r"\b(" + "|".join(cardiac_terms) + r")\b", re.IGNORECASE)

def find_keywords(text, pattern):
    if not text or text.lower() == "nan":
        return []
    return list(set(re.findall(pattern, text)))


df["renal"] = df.apply(
    lambda row: list(set(
        find_keywords(row["Clinical History"], renal_pattern) +
        find_keywords(row["FAMILY HISTORY"], renal_pattern)
    )),
    axis=1
)

df["cardio"] = df.apply(
    lambda row: list(set(
        find_keywords(row["Clinical History"], cardiac_pattern) +
        find_keywords(row["FAMILY HISTORY"], cardiac_pattern)
    )),
    axis=1
)


df_filtered = df[(df["renal"].str.len() > 0) | (df["cardio"].str.len() > 0)]


df_out = df_filtered[["CASE ID", "Clinical History", "FAMILY HISTORY", "renal", "cardio"]]


df_out.to_excel(r"C:\Users\verit\OneDrive\CoreGenomicsPrivate\renal_cardio_exome_clinical_history.xlsx", index=False)

print("Done")
