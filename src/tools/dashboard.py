import streamlit as st
import pandas as pd
from datetime import datetime

@st.cache_data
def get_data():
    return pd.DataFrame({
        "Name": ["Alice", "Bob", "Charlie"],
        "Score": [85, 92, 78],
        "Timestamp": [datetime.now().strftime("%H:%M:%S")] * 3,
    })

st.title("Live Dashboard")
df = get_data()
st.table(df.style.format({"Score": "{:.1f}"}))