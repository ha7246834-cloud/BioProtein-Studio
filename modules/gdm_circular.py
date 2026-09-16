from __future__ import annotations

import copy
import io
import math
import re
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from Bio import Phylo

from .gdm_style import style_from_preset, assign_colors


def _parse_tree(tree_text: str, keep: Iterable[str] | None = None, midpoint_root: bool = False):
    if not str(tree_text or '').strip():
        raise ValueError('Newick tree is required.')
    tree = Phylo.read(io.StringIO(str(tree_text).strip()), 'newick')
    tree = copy.deepcopy(tree)

    if keep is not None:
        wanted = {str(x) for x in keep}
        for terminal in list(tree.get_terminals()):
            if terminal.name not in wanted:
                try:
                    tree.prune(terminal)
                except Exception:
                    pass

    if len(tree.get_terminals()) < 2:
        raise ValueError('Circular phylogeny needs at least two matching terminal sequences.')

    if midpoint_root:
        try:
            tree.root_at_midpoint()
        except Exception:
            # Midpoint rooting is a display option. If branch lengths are not
            # suitable, retain the original rooting rather than altering data.
            pass
    return tree


def support_text(clade) -> str:
    """Return the support token without pretending one method is another.

    FastTree commonly stores local support as 0..1 numeric confidence. IQ-TREE
    may store a combined SH-aLRT/UFBoot token such as ``95/100`` in the internal
    node name. Both are preserved verbatim/faithfully for plotting.
    """
    conf = getattr(clade, 'confidence', None)
    if conf is not None:
        try:
            val = float(conf)
            return f'{val:.2f}' if val <= 1.0 else f'{val:.0f}'
        except Exception:
            pass

    name = str(getattr(clade, 'name', '') or '').strip()
    if re.fullmatch(r'\d+(?:\.\d+)?(?:/\d+(?:\.\d+)?)?', name):
        return name
    return ''


def _default_group_count(n: int) -> int:
    n = max(2, int(n))
    # A readable automatic target: enough groups to reveal major topology
    # without manufacturing dozens of pseudo-clades in large families.
    return min(8, max(2, int(round(math.sqrt(n)))))


def _first_tip_index(clade, tip_index: dict[str, int]) -> int:
    vals = [tip_index.get(t.name, 10**9) for t in clade.get_terminals()]
    return min(vals) if vals else 10**9


def auto_topology_clades(tree, target_groups: int | None = None):
    """Partition a tree into non-overlapping topology-derived display groups.

    These groups are deliberately labelled ``Auto Clade`` rather than assigned
    biological subgroup names. The algorithm recursively splits the largest
    current subtree until the requested/automatic number of topology groups is
    reached. It never uses gene names, species names, motifs, or annotations.
    """
    terminals = list(tree.get_terminals())
    n = len(terminals)
    if n < 2:
        return [], {}

    target = int(target_groups) if target_groups else _default_group_count(n)
    target = max(2, min(target, n, 12))
    tip_index = {t.name: i for i, t in enumerate(terminals)}

    groups = [c for c in tree.root.clades if c.count_terminals() > 0]
    if not groups:
        groups = [tree.root]

    # Usually Newick trees are bifurcating. For polytomies, do not exceed the
    # target merely to force an arbitrary split.
    while len(groups) < target:
        candidates = [
            g for g in groups
            if len(getattr(g, 'clades', [])) >= 2 and g.count_terminals() >= 2
        ]
        if not candidates:
            break
        split = max(candidates, key=lambda g: (g.count_terminals(), -_first_tip_index(g, tip_index)))
        children = [c for c in split.clades if c.count_terminals() > 0]
        if len(groups) - 1 + len(children) > target:
            break
        idx = groups.index(split)
        groups[idx:idx + 1] = children

    groups.sort(key=lambda g: _first_tip_index(g, tip_index))
    labels = []
    leaf_to_group = {}
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for i, clade in enumerate(groups):
        suffix = alphabet[i] if i < len(alphabet) else str(i + 1)
        label = f'Auto Clade {suffix}'
        labels.append((label, clade))
        for tip in clade.get_terminals():
            leaf_to_group[tip.name] = label

    return labels, leaf_to_group


def _depths(tree):
    d = tree.depths()
    if not d or max(d.values(), default=0) <= 0:
        d = tree.depths(unit_branch_lengths=True)
    return d


def _angles(tree):
    tips = list(tree.get_terminals())
    n = len(tips)
    tip_theta = {
        tip: (2.0 * math.pi * i / n)
        for i, tip in enumerate(tips)
    }
    theta = {}

    def visit(clade):
        if clade.is_terminal():
            theta[clade] = tip_theta[clade]
            return theta[clade]
        vals = [visit(child) for child in clade.clades]
        # Tree traversal keeps descendants contiguous; arithmetic mean therefore
        # gives a stable internal-node angle except for the root wrap, whose
        # exact angle is visually irrelevant because it has no parent branch.
        theta[clade] = float(sum(vals) / len(vals)) if vals else 0.0
        return theta[clade]

    visit(tree.root)
    return theta, tips


def _group_for_clade(clade, leaf_to_group: dict[str, str]):
    groups = {leaf_to_group.get(t.name) for t in clade.get_terminals()}
    groups.discard(None)
    return next(iter(groups)) if len(groups) == 1 else None


def circular_phylogeny_figure(
    tree_text: str,
    order: Iterable[str] | None = None,
    style=None,
    target_groups: int | None = None,
    midpoint_root: bool = False,
    show_support: bool = True,
    method_label: str = '',
):
    """Render a publication-oriented circular/radial phylogeny.

    Returns ``(figure, clade_table)``. Clade colours are topology-derived display
    groups only; they are not asserted to be functional/evolutionary subgroups.
    """
    tree = _parse_tree(tree_text, keep=order, midpoint_root=midpoint_root)
    tips = list(tree.get_terminals())
    n = len(tips)

    base = style_from_preset('Journal Classic')
    if style:
        base.update(style)

    group_defs, leaf_to_group = auto_topology_clades(tree, target_groups)
    palette = list(base.get('motif_palette') or ['#4C78A8', '#F58518', '#54A24B', '#E45756'])
    group_colors = assign_colors([name for name, _ in group_defs], palette)

    depth = _depths(tree)
    theta, tips = _angles(tree)
    rmax = max(depth.values(), default=1.0) or 1.0
    label_pad = max(rmax * 0.10, 0.08)
    ring_r = rmax + label_pad * 0.25
    label_r = rmax + label_pad

    size = min(16.0, max(8.5, 7.2 + 0.055 * n))
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    neutral = base.get('tree_color', '#222222')
    lw = float(base.get('line_width', 1.0))

    # Draw parent arcs and radial child branches. A branch inherits a group
    # colour only when all descendant terminals belong to the same auto clade;
    # deeper shared backbone branches remain neutral.
    for clade in tree.find_clades(order='preorder'):
        if clade not in theta:
            continue
        r0 = float(depth.get(clade, 0.0))
        children = [c for c in clade.clades if c in theta]
        if children:
            child_thetas = [theta[c] for c in children]
            lo, hi = min(child_thetas), max(child_thetas)
            arc = np.linspace(lo, hi, max(8, int(abs(hi - lo) * 30) + 2))
            group = _group_for_clade(clade, leaf_to_group)
            color = group_colors.get(group, neutral)
            ax.plot(arc, np.full_like(arc, r0), color=color, lw=lw, solid_capstyle='round')
            for child in children:
                cg = _group_for_clade(child, leaf_to_group)
                ccolor = group_colors.get(cg, neutral)
                ax.plot(
                    [theta[child], theta[child]],
                    [r0, float(depth.get(child, r0))],
                    color=ccolor,
                    lw=lw,
                    solid_capstyle='round',
                )

        if show_support and not clade.is_terminal():
            lab = support_text(clade)
            if lab:
                ax.text(
                    theta[clade],
                    r0 + label_pad * 0.05,
                    lab,
                    fontsize=float(base.get('support_size', 7)),
                    color=base.get('support_color', '#666666'),
                    ha='center',
                    va='center',
                )

    # Outer clade ring and labels.
    for tip in tips:
        group = leaf_to_group.get(tip.name)
        color = group_colors.get(group, neutral)
        t = theta[tip]
        ax.scatter([t], [ring_r], s=18, color=color, zorder=4)

        deg = math.degrees(t)
        rotation = deg - 90
        ha = 'left'
        if 90 < deg < 270:
            rotation += 180
            ha = 'right'
        ax.text(
            t,
            label_r,
            str(tip.name or ''),
            rotation=rotation,
            rotation_mode='anchor',
            ha=ha,
            va='center',
            fontsize=float(base.get('label_size', 9)),
            color='#222222',
        )

    # Draw a coloured outer arc for each topology group.
    for name, clade in group_defs:
        members = [t for t in clade.get_terminals() if t in theta]
        if not members:
            continue
        vals = sorted(theta[t] for t in members)
        arc = np.linspace(vals[0], vals[-1], max(12, len(vals) * 6))
        ax.plot(arc, np.full_like(arc, ring_r), color=group_colors[name], lw=4.0, alpha=0.9)

    ax.set_ylim(0, label_r + label_pad * 1.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)

    method = str(method_label or '').strip()
    title = 'Circular Phylogeny — Automatic Topology Clades'
    ax.set_title(title, fontsize=float(base.get('title_size', 13)) + 1, pad=24, weight='bold')
    if method:
        fig.text(0.5, 0.025, f'Inference/support source: {method}', ha='center', fontsize=9)
    fig.text(
        0.5,
        0.008,
        'Auto Clades are topology-derived display groups, not automatically validated biological subgroups.',
        ha='center',
        fontsize=8,
    )

    handles = [
        Line2D([0], [0], color=group_colors[name], lw=4, label=name)
        for name, _ in group_defs
    ]
    if handles:
        ax.legend(
            handles=handles,
            loc='center',
            bbox_to_anchor=(0.5, 0.5),
            frameon=False,
            fontsize=8,
            ncol=1 if len(handles) <= 5 else 2,
        )

    rows = []
    for name, clade in group_defs:
        members = [str(t.name or '') for t in clade.get_terminals()]
        rows.append({
            'auto_clade': name,
            'members': len(members),
            'root_support': support_text(clade),
            'member_ids': ';'.join(members),
        })
    clade_table = pd.DataFrame(rows)

    fig.tight_layout(rect=[0.02, 0.05, 0.98, 0.97])
    return fig, clade_table
