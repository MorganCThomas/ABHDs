"""Small plotting helpers shared by the notebooks."""
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from . import config as C


def scatter_abhd(coords, labels, is_abhd, title, ax=None, annotate=True):
    """Static scatter: ABHD-named entries in red (optionally labelled), everything else in grey."""
    ax = ax or plt.subplots(figsize=(7, 5))[1]
    bg = [i for i, a in enumerate(is_abhd) if not a]
    ax.scatter(coords[bg, 0], coords[bg, 1], c="lightgrey", s=18, label="other (non-ABHD name)")
    fg = [i for i, a in enumerate(is_abhd) if a]
    ax.scatter(coords[fg, 0], coords[fg, 1], c="crimson", s=22, label="ABHD-named")
    if annotate:
        for i in fg:
            ax.annotate(labels[i], coords[i], fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set(title=title, xlabel="component 1", ylabel="component 2")
    ax.legend(frameon=False, fontsize=8)
    return ax


def interactive_scatter(coords, hover, groups, title, path):
    """Plotly HTML with hover text; `groups` maps group name -> (indices, colour). Saved under figures/."""
    fig = go.Figure()
    for name, (idx, colour) in groups.items():
        fig.add_trace(go.Scatter(x=coords[idx, 0], y=coords[idx, 1], mode="markers", name=name,
                                 marker=dict(color=colour, size=8), text=[hover[i] for i in idx],
                                 hovertemplate="%{text}<extra></extra>"))
    fig.update_layout(title=title, template="simple_white", width=900, height=650)
    C.FIGURES.mkdir(exist_ok=True)
    fig.write_html(C.FIGURES / path)
    return fig
