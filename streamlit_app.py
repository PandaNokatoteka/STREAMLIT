import os
import scanpy as sc
import streamlit as st
import gdown
import matplotlib.pyplot as plt

# Βασικές ρυθμίσεις εμφάνισης
st.set_page_config(layout="wide")
st.title("🔬 Single-cell RNA-seq Viewer")

DEMO_FILE_ID = "1XybfO8QZ0G3gigwwzHk-n8gopuT4iEHB"
DEMO_FILENAME = "pancreas_data.h5ad"


@st.cache_data(show_spinner=False)
def get_demo_data():
    """Λήψη και ανάγνωση demo αρχείου με ασφάλεια"""
    if not os.path.exists(DEMO_FILENAME):
        st.info("📥 Κατεβάζουμε demo αρχείο από Google Drive...")
        url = f"https://drive.google.com/uc?id={DEMO_FILE_ID}"
        try:
            gdown.download(url, DEMO_FILENAME, quiet=False)
        except Exception as e:
            st.error(f"❌ Αποτυχία λήψης αρχείου: {e}")
            st.stop()

    try:
        adata = sc.read_h5ad(DEMO_FILENAME)
        return adata
    except Exception as e:
        st.error(f"❌ Σφάλμα κατά την ανάγνωση του αρχείου: {e}")
        st.stop()


# Από εδώ ξεκινά η εφαρμογή
st.sidebar.header("🧬 Ρυθμίσεις δεδομένων")

# Δυνατότητα χρήσης demo ή ανέβασμα custom αρχείου
use_demo = st.sidebar.checkbox("Χρήση demo dataset", value=True)

if use_demo:
    adata = get_demo_data()
else:
    uploaded_file = st.sidebar.file_uploader("Ανέβασε .h5ad αρχείο", type=["h5ad"])
    if uploaded_file is not None:
        try:
            adata = sc.read_h5ad(uploaded_file)
            # Προετοιμασία UMAP αν δεν υπάρχει
if "X_umap" not in adata.obsm:
    st.info("🔄 Υπολογίζουμε UMAP coordinates...")
    try:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    except Exception as e:
        st.error(f"❌ Σφάλμα κατά τον υπολογισμό UMAP: {e}")
        st.stop()

except Exception as e:
    st.error(f"❌ Σφάλμα κατά την ανάγνωση του αρχείου: {e}")
    st.stop()
    else:
        st.warning("📂 Περιμένουμε αρχείο .h5ad...")
        st.stop()


# Εμφάνιση βασικών πληροφοριών
st.subheader("📊 Βασικές πληροφορίες dataset")
st.markdown(f"**Σχήμα δεδομένων:** {adata.X.shape[0]} κελιά × {adata.X.shape[1]} γονίδια")
st.markdown(f"**Διαθέσιμα AnnData.obs:** {', '.join(adata.obs.columns)}")
st.markdown(f"**Διαθέσιμα AnnData.var:** {', '.join(adata.var.columns)}")

# Επιλογή cluster (αν υπάρχει)
if "louvain" in adata.obs.columns:
    st.subheader("🧩 Cluster visualization")
    fig, ax = plt.subplots()
    sc.pl.umap(adata, color="louvain", ax=ax, show=False)
    st.pyplot(fig)

# Επιλογή γονιδίων για plotting
st.subheader("🎯 Gene expression")
gene_list = st.multiselect("Διάλεξε γονίδια:", list(adata.var_names), default=adata.var_names[:1])

for gene in gene_list:
    if gene in adata.var_names:
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=gene, ax=ax, show=False)
        st.pyplot(fig)
    else:
        st.warning(f"⚠️ Το γονίδιο {gene} δεν υπάρχει στο dataset.")

