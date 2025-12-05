import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import sys
from pathlib import Path

# Add project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / '01_Scripts'))

from Python.config import resolve_path
from Python.viz_bump_charts import load_data, filter_by_scope
from Python.pattern_definitions import get_pattern_colors

# Page Config
st.set_page_config(layout="wide", page_title="DRP1 Pathway Explorer")

# Load Data (Cached)
@st.cache_data
def get_data():
    return load_data()

df = get_data()

# Sidebar Controls
st.sidebar.title("Configuration")

mutation = st.sidebar.selectbox("Mutation", ["G32A", "R403C"], index=0)
scope = st.sidebar.selectbox("Scope", ["focused", "significant", "all"], index=0)
y_type = st.sidebar.radio("Y-Axis", ["NES", "Rank"], index=0)

# Filter Data
df_scope = filter_by_scope(df, scope)
pattern_col = f'Pattern_{mutation}'

# Highlighting Mechanism
st.sidebar.markdown("### Highlight Pathways")
search_term = st.sidebar.text_input("Search (Regex supported)", "")

# Filter available options based on search
if search_term:
    matches = df_scope[df_scope['Description'].str.contains(search_term, case=False, na=False)]
    options = matches['pathway_id'].tolist()
else:
    options = df_scope['pathway_id'].tolist()

# Multiselect for specific highlights
selected_ids = st.sidebar.multiselect(
    "Select specific pathways",
    options=options,
    format_func=lambda x: df_scope[df_scope['pathway_id'] == x]['Description'].values[0][:50]
)

# Plotting Logic
st.title(f"Pathway Trajectory: {mutation}")
st.markdown(f"**Scope:** {scope} | **Pathways:** {len(df_scope)}")

# Plotly Chart
fig = go.Figure()

stages = ['Early', 'TrajDev', 'Late']
nes_cols = [f'NES_{s}_{mutation}' for s in stages]
x_vals = [0, 1, 2]
x_labels = ['Early', 'TrajDev', 'Late']

# 1. Background Traces (Faint)
# Group by pattern for efficient rendering
patterns = df_scope[pattern_col].unique()
colors = get_pattern_colors()

for pattern in patterns:
    df_pat = df_scope[df_scope[pattern_col] == pattern]
    # Exclude selected if any
    if selected_ids:
        df_pat = df_pat[~df_pat['pathway_id'].isin(selected_ids)]
        opacity = 0.1 # Very faint if highlighting
    else:
        opacity = 0.4 # Normal if no highlight

    # For speed, use fewer traces or sampling if 'all'
    if scope == 'all' and len(df_pat) > 500 and not selected_ids:
        df_pat = df_pat.sample(500) # Sample for performance if no filter

    for _, row in df_pat.iterrows():
        y = row[nes_cols].values
        if y_type == 'Rank':
            # Need to compute ranks dynamically or use pre-computed
            # For this demo, let's stick to NES or implement simple ranking
            pass 
        
        fig.add_trace(go.Scatter(
            x=x_vals, y=y,
            mode='lines',
            line=dict(color=colors.get(pattern, 'gray'), width=1),
            opacity=opacity,
            hoverinfo='skip', # optimize
            showlegend=False
        ))

# 2. Highlighted Traces (Thick, Opaque)
if selected_ids:
    df_hi = df_scope[df_scope['pathway_id'].isin(selected_ids)]
    for _, row in df_hi.iterrows():
        pattern = row[pattern_col]
        y = row[nes_cols].values
        
        desc = row['Description']
        hover_text = f"<b>{desc}</b><br>{pattern}<br>NES: {y}"
        
        fig.add_trace(go.Scatter(
            x=x_vals, y=y,
            mode='lines+markers',
            line=dict(color=colors.get(pattern, 'gray'), width=4),
            opacity=1.0,
            name=desc[:30],
            text=hover_text,
            hoverinfo='text'
        ))

fig.update_layout(
    xaxis=dict(tickmode='array', tickvals=x_vals, ticktext=x_labels),
    yaxis=dict(title="NES"),
    height=600,
    template='plotly_white'
)

st.plotly_chart(fig, use_container_width=True)

# Data Table
if selected_ids:
    st.markdown("### Selected Data")
    cols = ['pathway_id', 'Description', pattern_col] + nes_cols
    st.dataframe(df_hi[cols])
else:
    st.info("Select pathways in the sidebar to highlight them and see details.")
