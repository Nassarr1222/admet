import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt

st.title("🧪 ADMET Analysis Tool")

uploaded_file = st.file_uploader("📄 Upload your SMILES .txt file", type="txt")

def get_admet(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return [None]*7
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hba = Descriptors.NumHAcceptors(mol)
    hbd = Descriptors.NumHDonors(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    lipinski = "✅" if (mw <= 500 and logp <= 5 and hba <= 10 and hbd <= 5) else "❌"
    return [mw, logp, tpsa, hba, hbd, rotb, lipinski]

if uploaded_file:
    smiles_list = [line.strip() for line in uploaded_file if line.strip()]
    smiles_str_list = [s.decode("utf-8") for s in smiles_list]
    df = pd.DataFrame({"SMILES": smiles_str_list})
    admet_data = [get_admet(s) for s in df["SMILES"]]
    admet_df = pd.DataFrame(admet_data, columns=["MW", "LogP", "TPSA", "HBA", "HBD", "RotB", "Lipinski"])
    final_df = pd.concat([df, admet_df], axis=1)

    st.write("### 📊 ADMET Results")
    st.dataframe(final_df)

    st.download_button("📥 Download Excel", data=final_df.to_csv(index=False), file_name="admet_results.csv")

    st.write("### 🔬 LogP Plot")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(final_df["SMILES"], final_df["LogP"], color="skyblue")
    ax.set_xlabel("LogP")
    ax.set_title("Lipophilicity (LogP) of Compounds")
    st.pyplot(fig)