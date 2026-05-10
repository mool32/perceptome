"""Module co-variation network analysis."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from ..utils import cosine


def compute_network(
    scores_per_cell,
    cell_type_column=None,
    obs=None,
    method='spearman',
    threshold=0.1,
):
    """Compute module co-variation network from single-cell data.

    Parameters
    ----------
    scores_per_cell : DataFrame
        Individual cells x modules.
    cell_type_column : str, optional
        Compute per-type networks if provided (requires obs).
    obs : Series, optional
        Cell type assignments.
    method : str
        'spearman' or 'pearson'.
    threshold : float
        Correlation threshold for edges.

    Returns
    -------
    dict with keys:
        'correlation_matrix': DataFrame (modules x modules)
        'adjacency': DataFrame (thresholded binary)
        'hub_degrees': Series (degree per module)
        'communities': dict (community assignments via spectral clustering)
        'rich_club_phi': dict (phi at each degree threshold)
        'rich_club_intact': bool (are top 5 hubs fully connected?)
    """
    modules = list(scores_per_cell.columns)
    X = scores_per_cell.values

    # Correlation matrix
    if method == 'spearman':
        corr_mat, _ = spearmanr(X)
    else:
        corr_mat = np.corrcoef(X.T)

    if corr_mat.ndim == 0:
        corr_mat = np.array([[corr_mat]])

    corr_df = pd.DataFrame(corr_mat, index=modules, columns=modules)

    # Adjacency (thresholded)
    adj = (np.abs(corr_mat) > threshold).astype(int)
    np.fill_diagonal(adj, 0)
    adj_df = pd.DataFrame(adj, index=modules, columns=modules)

    # Hub degrees
    degrees = adj.sum(axis=1)
    hub_degrees = pd.Series(degrees, index=modules, name='degree').sort_values(ascending=False)

    # Rich club analysis
    rich_club_phi = {}
    unique_degrees = sorted(set(degrees))
    for k in unique_degrees:
        if k == 0:
            continue
        rich_nodes = np.where(degrees >= k)[0]
        n_rich = len(rich_nodes)
        if n_rich < 2:
            continue
        edges_among_rich = sum(
            adj[i, j] for i in rich_nodes for j in rich_nodes if i < j
        )
        max_edges = n_rich * (n_rich - 1) / 2
        rich_club_phi[int(k)] = float(edges_among_rich / max_edges) if max_edges > 0 else 0

    # Check if top 5 hubs are fully connected
    top5_idx = np.argsort(degrees)[-5:]
    edges_top5 = sum(adj[i, j] for i in top5_idx for j in top5_idx if i < j)
    rich_club_intact = bool(edges_top5 == 10)  # 5 choose 2 = 10

    # Simple community detection via thresholded correlation
    from sklearn.cluster import SpectralClustering
    try:
        n_clusters = min(5, len(modules))
        affinity = np.abs(corr_mat)
        np.fill_diagonal(affinity, 0)
        clustering = SpectralClustering(
            n_clusters=n_clusters,
            affinity='precomputed',
            assign_labels='kmeans',
            random_state=42,
        ).fit(affinity)
        communities = {modules[i]: int(clustering.labels_[i]) for i in range(len(modules))}
    except Exception:
        communities = {m: 0 for m in modules}

    result = {
        'correlation_matrix': corr_df,
        'adjacency': adj_df,
        'hub_degrees': hub_degrees,
        'communities': communities,
        'rich_club_phi': rich_club_phi,
        'rich_club_intact': rich_club_intact,
    }

    # Per-cell-type networks
    if cell_type_column and obs is not None:
        per_ct = {}
        for ct in obs.unique():
            mask = obs == ct
            if mask.sum() < 50:
                continue
            ct_result = compute_network(scores_per_cell.loc[mask], method=method, threshold=threshold)
            per_ct[ct] = ct_result
        result['per_cell_type'] = per_ct

    return result


def compare_networks(network_disease, network_control):
    """Compare disease vs control networks.

    Parameters
    ----------
    network_disease : dict
        From compute_network (disease).
    network_control : dict
        From compute_network (control).

    Returns
    -------
    dict with keys:
        'np_score': float (total network perturbation)
        'edge_changes': DataFrame (per edge: delta correlation)
        'hub_changes': DataFrame (per module: delta degree)
        'decoupled_pairs': list (module pairs that lost correlation)
        'recoupled_pairs': list (module pairs that gained correlation)
        'rich_club_change': float (delta phi for core)
        'narrative': str (human-readable summary)
    """
    corr_d = network_disease['correlation_matrix']
    corr_c = network_control['correlation_matrix']

    modules = list(corr_d.index)
    delta_corr = corr_d - corr_c

    # Network perturbation score
    np_score = float(np.sqrt(np.sum(delta_corr.values ** 2)) / 2)  # /2 for symmetry

    # Edge changes
    edge_rows = []
    decoupled = []
    recoupled = []
    for i in range(len(modules)):
        for j in range(i + 1, len(modules)):
            d = delta_corr.iloc[i, j]
            c_val = corr_c.iloc[i, j]
            d_val = corr_d.iloc[i, j]
            edge_rows.append({
                'module_1': modules[i],
                'module_2': modules[j],
                'corr_control': float(c_val),
                'corr_disease': float(d_val),
                'delta': float(d),
            })
            # Decoupled: was correlated, now not
            if abs(c_val) > 0.2 and abs(d_val) < 0.1:
                decoupled.append((modules[i], modules[j]))
            # Recoupled: was not correlated, now is
            if abs(c_val) < 0.1 and abs(d_val) > 0.2:
                recoupled.append((modules[i], modules[j]))

    edge_df = pd.DataFrame(edge_rows)

    # Hub changes
    deg_d = network_disease['hub_degrees']
    deg_c = network_control['hub_degrees']
    hub_changes = (deg_d - deg_c).sort_values(ascending=False)
    hub_changes.name = 'delta_degree'

    # Rich club change
    phi_d = network_disease.get('rich_club_phi', {})
    phi_c = network_control.get('rich_club_phi', {})
    common_k = set(phi_d.keys()) & set(phi_c.keys())
    if common_k:
        max_k = max(common_k)
        rc_change = phi_d[max_k] - phi_c[max_k]
    else:
        rc_change = 0.0

    # Narrative
    n_decoupled = len(decoupled)
    n_recoupled = len(recoupled)
    top_gained = hub_changes.head(3)
    top_lost = hub_changes.tail(3)

    narrative = (
        f"Network perturbation score: {np_score:.3f}. "
        f"{n_decoupled} module pairs decoupled, {n_recoupled} recoupled. "
        f"Gained connectivity: {', '.join(top_gained.index[:3])}. "
        f"Lost connectivity: {', '.join(top_lost.index[:3])}. "
        f"Rich club change: {rc_change:+.3f}."
    )

    return {
        'np_score': np_score,
        'edge_changes': edge_df,
        'hub_changes': hub_changes,
        'decoupled_pairs': decoupled,
        'recoupled_pairs': recoupled,
        'rich_club_change': float(rc_change),
        'narrative': narrative,
    }
