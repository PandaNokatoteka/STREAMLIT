import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import os
import requests
import scRNA_seq_pipeline

st.set_page_config(page_title="scRNA-seq Pipeline", layout="wide")
st.title("🔬 Διαδραστική Ανάλυση scRNA-seq Δεδομένων")

st.markdown("Ανεβάστε δεδομένα ή χρησιμοποιήστε το demo dataset (pancreas) από Google Drive.")

# --- Demo dataset config ---
DEMO_URL = "https://drive.google.com/file/d/1XybfO8QZ0G3gigwwzHk-n8gopuT4iEHB/view?usp=drive_link"
DEMO_FILENAME = "pancreas_data.h5ad"

if data_source == "Demo Dataset":
    if not os.path.exists(DEMO_FILENAME):
        st.info("📥 Γίνεται λήψη του demo αρχείου από Google Drive...")
        with open(DEMO_FILENAME, "wb") as f:
            response = requests.get(DEMO_URL)
            f.write(response.content)
        st.success("✅ Το demo dataset κατέβηκε με επιτυχία.")
    adata_path = DEMO_FILENAME

else:
    uploaded_file = st.file_uploader("Ανεβάστε αρχείο .h5ad", type=["h5ad"])
    if uploaded_file is not None:
        adata_path = uploaded_file

if 'adata_path' in locals():
    st.sidebar.header("⚙️ Ρυθμίσεις Pipeline")
    normalize = st.sidebar.checkbox("Normalize Data", value=True)
    log1p = st.sidebar.checkbox("Log1p Transform", value=True)
    min_genes = st.sidebar.slider("Ελάχιστα γονίδια ανά κύτταρο", 0, 200, 50)
    resolution = st.sidebar.slider("Clustering Resolution", 0.1, 2.0, 1.0, step=0.1)

    if st.button("🔁 Εκτέλεση Pipeline"):
        with st.spinner("⏳ Εκτελείται το pipeline..."):
            adata = run_pipeline(
                adata_path=adata_path,
                normalize=normalize,
                log1p=log1p,
                min_genes=min_genes,
                resolution=resolution
            )

        st.success("✅ Ολοκληρώθηκε η ανάλυση!")

        st.subheader("📊 PCA Projection")
        sc.pl.pca(adata, show=False)
        st.pyplot(plt.gcf())

        st.subheader("🧭 UMAP Projection")
        sc.pl.umap(adata, color=['leiden'], show=False)
        st.pyplot(plt.gcf())

        st.subheader("🧬 Διαφορική Έκφραση")
        if "rank_genes_groups" in adata.uns:
            sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False)
            st.pyplot(plt.gcf())
        else:
            st.warning("⚠️ Δεν εντοπίστηκε ανάλυση διαφορικής έκφρασης.")
else:
    st.info("📥 Παρακαλώ ανεβάστε ή επιλέξτε δεδομένα για ανάλυση.")
