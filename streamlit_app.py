
import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
from scRNA_seq_pipeline import run_pipeline


st.set_page_config(page_title="scRNA-seq Pipeline", layout="wide")

st.title("🔬 Διαδραστική Ανάλυση scRNA-seq Δεδομένων")
st.markdown("Ανεβάστε δεδομένα ή χρησιμοποιήστε το demo dataset για ανάλυση.")

# Επιλογή αρχείου ή χρήση demo
data_source = st.radio("Επιλογή δεδομένων:", ("Demo Dataset", "Ανέβασμα αρχείου (.h5ad)"))

if data_source == "Demo Dataset":
    adata_path = "pancreas_data.h5ad"
else:
    uploaded_file = st.file_uploader("Ανεβάστε .h5ad αρχείο", type=["h5ad"])
    if uploaded_file is not None:
        adata_path = uploaded_file

if 'adata_path' in locals():
    st.success("Δεδομένα φορτώθηκαν επιτυχώς!")

    # Επιλογές για preprocessing και ανάλυση
    st.sidebar.header("⚙️ Ρυθμίσεις Pipeline")
    normalize = st.sidebar.checkbox("Normalize Data", value=True)
    log1p = st.sidebar.checkbox("Log1p Transform", value=True)
    min_genes = st.sidebar.slider("Ελάχιστα γονίδια ανά κύτταρο", 0, 200, 50)
    resolution = st.sidebar.slider("Clustering Resolution", 0.1, 2.0, 1.0, step=0.1)

    if st.button("🔁 Εκτέλεση Pipeline"):
        with st.spinner("Εκτελείται το pipeline..."):
            adata = run_pipeline(
                adata_path=adata_path,
                normalize=normalize,
                log1p=log1p,
                min_genes=min_genes,
                resolution=resolution
            )

        st.success("✅ Ολοκληρώθηκε η ανάλυση!")

        # Παρουσίαση αποτελεσμάτων
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
            st.warning("Δεν εντοπίστηκε ανάλυση διαφορικής έκφρασης.")
else:
    st.info("Παρακαλώ επιλέξτε αρχείο δεδομένων για ανάλυση.")
