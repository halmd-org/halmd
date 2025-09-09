#!/usr/bin/env python3
# compare_h5md.py

import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt


# ---------- Helpers ----------
def get_values(f, grp_path, obs):
    """Return full value array (no reshaping) or None if missing."""
    if grp_path not in f or obs not in f[grp_path]:
        return None
    node = f[grp_path][obs]
    if "value" not in node:
        return None
    return np.asarray(node["value"], dtype=float)


def get_time_series(f, grp_path, obs):
    """Return (time, values) for 1-D observables in observables/*."""
    if grp_path not in f or obs not in f[grp_path]:
        return None, None
    node = f[grp_path][obs]
    if "value" not in node:
        return None, None
    v = np.asarray(node["value"], dtype=float).reshape(-1)
    if "time" in node:
        t = np.asarray(node["time"], dtype=float).reshape(-1)
    elif "step" in node:
        t = np.asarray(node["step"], dtype=float).reshape(-1)
    else:
        t = np.arange(v.shape[0], dtype=float)
    assert len(t) == len(v), "time and values must have same length"
    return t, v


def print_numeric_check(name1, name2, grp_label, obs, t1, v1, t2, v2, tol):
    """Check 1-D series (observables)."""
    assert v1.shape == v2.shape, f"{obs}: arrays must have same length"
    n = v1.size
    diffs = np.abs(v1 - v2)
    same = np.allclose(v1, v2, atol=tol, rtol=0.0)
    i_max = diffs.argmax()
    t_at_max = float(t1[i_max]) if n > 0 else np.nan

    status = "OK" if same else "MISMATCH"
    reason = f"(max |Δ| {diffs[i_max]:.3e}{' ≤ ' if same else ' > '}{tol:g})"
    print(f"- {grp_label} | {obs}: {status} {reason}")
    if n > 0:
        mean1, mean2 = np.mean(v1), np.mean(v2)
        print(f"    length={n}, max |Δ| at index {i_max} (t≈{t_at_max:.6g})")
        print(f"    {name1} avg={mean1:.6g}, {name2} avg={mean2:.6g}, |Δavg|={abs(mean2-mean1):.6g}")
    return bool(same)


def print_energy_per_particle_check(name1, name2, grp_label, e1, e2, tol):
    """Compare potential energy per particle: arrays (nsteps, nparticles)."""
    assert e1.shape == e2.shape, "Potential energy arrays must match in shape"
    nsteps, npart = e1.shape
    diffs = np.abs(e1 - e2)
    max_diff = float(diffs.max())
    mean_diff = float(diffs.mean())
    same = max_diff <= tol

    status = "OK" if same else "MISMATCH"
    print(f"- {grp_label} | potential_energy: {status} (max |Δ| {max_diff:.3e} ≤ {tol:g}?)")
    print(f"    mean |Δ|={mean_diff:.3e}, steps={nsteps}, particles={npart}")
    return bool(same)


def print_trajectory_check(name1, name2, grp_label, pos1, pos2, tol):
    """Compare trajectories: arrays (nsteps, nparticles, ndim)."""
    assert pos1.shape == pos2.shape, "Trajectory arrays must match in shape"
    nsteps, npart, ndim = pos1.shape
    diffs = np.linalg.norm(pos1 - pos2, axis=-1)  # per-particle distance per step
    rmsd_per_step = np.sqrt(np.mean(diffs**2, axis=1))
    max_rmsd = float(rmsd_per_step.max())
    mean_rmsd = float(rmsd_per_step.mean())
    same = max_rmsd <= tol

    status = "OK" if same else "MISMATCH"
    print(f"- {grp_label} | trajectory: {status} (max RMSD {max_rmsd:.3e} ≤ {tol:g}?)")
    print(f"    mean RMSD={mean_rmsd:.3e}, steps={nsteps}, particles={npart}")
    return bool(same)


import matplotlib.cm as cm

def plot_comparison(name1, name2, grp_label, series, show=True, save_prefix=None):
    """
    Plot comparisons:
    - Observables (temperature, pressure, internal_energy, potential_energy 1D): overlay
    - Particles/potential_energy (2D): mean per particle + heatmap of |Δ|
    - Particles/position (3D): RMSD per step
    """
    nplots = len(series)
    fig, axes = plt.subplots(nplots, 1, figsize=(10, 4 * nplots))
    if nplots == 1:
        axes = [axes]

    for ax, (obs, data) in zip(axes, series.items()):
        if obs in ("temperature", "pressure", "internal_energy", "potential_energy"):
            # ---------- Observables ----------
            if isinstance(data, tuple) and len(data) == 4:
                t1, v1, t2, v2 = data
                ax.plot(t1, v1, label=name1)
                ax.plot(t2, v2, "--", label=name2)
                ax.set_title(f"{grp_label} | {obs}")
                ax.set_xlabel("time")
                ax.legend()
            elif isinstance(data, tuple) and data[0].ndim == 2:
                e1, e2 = data
                mean1, mean2 = e1.mean(axis=1), e2.mean(axis=1)
                diffs = np.abs(e1 - e2)

                # plot mean trace
                ax.plot(mean1, label=f"{name1} mean")
                ax.plot(mean2, "--", label=f"{name2} mean")
                ax.set_title(f"{grp_label} | Potential Energy (mean per particle)")
                ax.set_xlabel("step")
                ax.legend()

                # add a heatmap of differences in a new figure
                fig2, ax2 = plt.subplots(figsize=(8, 4))
                im = ax2.imshow(diffs.T, aspect="auto", origin="lower",
                                cmap=cm.viridis, interpolation="nearest")
                ax2.set_title(f"{grp_label} | |ΔE| per particle")
                ax2.set_xlabel("step")
                ax2.set_ylabel("particle index")
                fig2.colorbar(im, ax=ax2, label="|ΔE|")
                plt.tight_layout()
                if save_prefix:
                    out = f"{save_prefix}_{grp_label.replace(' ', '')}_PE_heatmap.png"
                    fig2.savefig(out, dpi=150)
                    print(f"Saved figure: {out}")
                if show:
                    plt.show()
                else:
                    plt.close(fig2)

        # ---------- Trajectories ----------
        elif obs == "position":
            pos1, pos2 = data
            diffs = np.linalg.norm(pos1 - pos2, axis=-1)
            rmsd_per_step = np.sqrt(np.mean(diffs**2, axis=1))
            ax.plot(rmsd_per_step)
            ax.set_title(f"{grp_label} | Trajectory RMSD per step")
            ax.set_xlabel("step")
            ax.set_ylabel("RMSD")

    plt.tight_layout()
    if save_prefix:
        out = f"{save_prefix}_{grp_label.replace(' ', '')}.png"
        fig.savefig(out, dpi=150)
        print(f"Saved figure: {out}")
    if show:
        plt.show()
    else:
        plt.close(fig)



def parse_groups(s):
    """Parse --groups string: path[:label][,path[:label],...]."""
    groups = {}
    for chunk in s.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" in chunk:
            path, label = chunk.split(":", 1)
        else:
            path = chunk
            label = path.split("/")[-1]
        groups[path.strip()] = label.strip()
    return groups


# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Compare H5MD observables and particle data between two files."
    )
    ap.add_argument("files", nargs=2, help="Two HDF5 files to compare")
    ap.add_argument(
        "--groups",
        default="observables/A:ObsA,observables/B:ObsB,particles/A:PartA,particles/B:PartB",
        help="Comma-separated group paths with optional labels."
    )
    ap.add_argument(
        "--observables",
        nargs='+',
        default=("temperature", "pressure", "potential_energy", "internal_energy", "position"),
        help="Observables/quantities to compare."
    )
    ap.add_argument("--tol", type=float, default=1e-2, help="Absolute tolerance.")
    ap.add_argument("--no-plot", action="store_true", help="Skip plotting.")
    ap.add_argument("--save-prefix", default=None, help="If set, save figures with this prefix.")
    args = ap.parse_args()

    file1, file2 = args.files
    name1, name2 = Path(file1).name, Path(file2).name

    for p in (file1, file2):
        if not Path(p).exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            return 5

    groups = parse_groups(args.groups)
    obs_try = args.observables
    tol = args.tol

    any_compared = False
    any_mismatch = False

    with h5py.File(file1, "r") as f1, h5py.File(file2, "r") as f2:
        for grp_path, grp_label in groups.items():
            common = []
            for obs in obs_try:
                if grp_path.startswith("observables/"):
                    t1, v1 = get_time_series(f1, grp_path, obs)
                    t2, v2 = get_time_series(f2, grp_path, obs)
                    if v1 is not None and v2 is not None:
                        common.append(obs)
                else:  # particles/
                    v1 = get_values(f1, grp_path, obs)
                    v2 = get_values(f2, grp_path, obs)
                    if v1 is not None and v2 is not None:
                        common.append(obs)

            if not common:
                print(f"\n[{grp_label}] No common observables to compare.")
                continue

            print(f"\n--- Numeric check @ {grp_label} ({name1} vs {name2}, tol={tol:g}) ---")
            series = {}
            for obs in common:
                if grp_path.startswith("observables/"):
                    t1, v1 = get_time_series(f1, grp_path, obs)
                    t2, v2 = get_time_series(f2, grp_path, obs)
                    series[obs] = (t1, v1, t2, v2)
                    ok = print_numeric_check(name1, name2, grp_label, obs, t1, v1, t2, v2, tol)
                else:  # particles
                    v1 = get_values(f1, grp_path, obs)
                    v2 = get_values(f2, grp_path, obs)
                    if obs == "potential_energy":
                        series[obs] = (v1, v2)
                        ok = print_energy_per_particle_check(name1, name2, grp_label, v1, v2, tol)
                    elif obs == "position":
                        series[obs] = (v1, v2)
                        ok = print_trajectory_check(name1, name2, grp_label, v1, v2, tol)
                    else:
                        continue

                any_compared = True
                if not ok:
                    any_mismatch = True

            if not args.no_plot:
                plot_comparison(
                    name1, name2, grp_label, series,
                    show=(args.save_prefix is None), save_prefix=args.save_prefix
                )

    if not any_compared:
        print("ERROR: No observables were compared.", file=sys.stderr)
        return 4
    if any_mismatch:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
