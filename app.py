import streamlit as st
from utils import construir_cucurbituril_completo, showm
import tempfile

st.set_page_config(layout="wide")
st.title("Generador de Cucurbituril")

st.sidebar.markdown("### Parámetros")
n_units = st.sidebar.slider("Número de unidades de glicolurilo", min_value=3, max_value=12, value=6)
style = st.sidebar.selectbox("Estilo de visualización", ["stick", "sphere", "line", "cartoon"])

uploaded_file = st.file_uploader("Carga el archivo MOL de una unidad de glicolurilo (.mol)", type=["mol"])

if uploaded_file:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".mol") as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = tmp.name

    try:
        mol_cb = construir_cucurbituril_completo(tmp_path, n_unidades=n_units)
        view = showm(mol_cb, style=style)
        st.components.v1.html(view._make_html(), height=500, width=600)
    except Exception as e:
        st.error(f"Ocurrió un error: {e}")
else:
    st.info("Por favor, carga un archivo MOL para comenzar.")
