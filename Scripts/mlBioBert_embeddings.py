import pandas as pd
from sentence_transformers import SentenceTransformer, util

df = pd.read_excel("exome_clinical_history.xlsx")

df["full_history"] = df["Clinical History"].astype(str) + " " + df["FAMILY HISTORY"].astype(str)

model = SentenceTransformer("pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb")

history_embeddings = model.encode(df["full_history"].tolist(), convert_to_tensor=True)

queries = {
    "renal": "kidney disease, renal failure, nephropathy, dialysis",
    "cardio": "heart disease, coronary artery disease, cardiomyopathy, heart failure"
}

query_embeddings = {label: model.encode(text, convert_to_tensor=True) for label, text in queries.items()}

df["renal_score"] = util.cos_sim(history_embeddings, query_embeddings["renal"]).cpu().numpy().max(axis=1)
df["cardio_score"] = util.cos_sim(history_embeddings, query_embeddings["cardio"]).cpu().numpy().max(axis=1)

df["renal_flag"] = df["renal_score"] > 0.45
df["cardio_flag"] = df["cardio_score"] > 0.45
df["renal_cardio_flag"] = df["renal_flag"] | df["cardio_flag"]

df_out = df[df["renal_cardio_flag"]][["CASE ID", "Clinical History", "FAMILY HISTORY", "renal_score", "cardio_score", "renal_flag", "cardio_flag"]]
df_out.to_excel("exome_clinical_history_ml_output_filtered.xlsx", index=False)

print("ML-based detection done. File saved as ml_output_filtered.xlsx")
