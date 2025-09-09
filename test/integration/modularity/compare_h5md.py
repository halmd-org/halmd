#!/usr/bin/env python3
# compare_h5md.py
import argparse
import os
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib.pyplot as plt


# functions to help
def get_time_values(f, grp_path, obs):
    """Return (time, values) or (None, None) if missing."""
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

    n = min(len(t), len(v))
    return t[:n], v[:n]


def print_numeric_check(name1, name2, grp_label, obs, t1, v1, t2, v2, tol):

    """print summary and return True if within tol, else False."""

    n = min(v1.size, v2.size)
    a, b = v1[:n], v2[:n]
    diffs = np.abs(a - b)
    same = np.allclose(a, b, atol=tol, rtol=0.0)        #returns True if two arrays are element-wise equal within a tolerance.
    max_diff = float(diffs.max()) if n > 0 else np.nan
    i_max = int(diffs.argmax()) if n > 0 else -1        #index where the max abs difference occurs
    t_display = t1 if len(t1) >= n else np.arange(n)
    t_at_max = float(t_display[i_max]) if n > 0 else np.nan

    status = "OK" if same else "MISMATCH"
    reason = f"(max |diff| {max_diff:.3e}{' ≤ ' if same else ' > '}{tol:g})"
    print(f"- {grp_label} | {obs}: {status} {reason}")
    if n == 0:
        print("no overlapping points (one series empty?).")
    else:
        mean1, mean2 = np.mean(a), np.mean(b)
        print(f"length compared: {n}, max |diff| at index {i_max} ($t \sim$ {t_at_max:.6g})")
        print(f"{name1} avg.={mean1:.6g}, {name2} avg.={mean2:.6g}, $\Delta =${abs(mean2-mean1):.6g}")
    return bool(same)


def plot_comparison(name1, name2, grp_label, obs_list, series, show=True, save_prefix=None):
    """
    series is dict[obs] = (t1, v1, t2, v2)
    Creates rows: one per observable; 2 columns: overlay and absolute difference.
    """
    nrows = len(obs_list)
    fig, axes = plt.subplots(nrows, 2, figsize=(11, 3.0 * nrows), sharex="col")
    if nrows == 1:
        axes = np.array([axes])

    fig.suptitle(f"{grp_label} — {name1} vs {name2}", fontsize=13)

    for r, obs in enumerate(obs_list):
        t1, v1, t2, v2 = series[obs]

        # Left: overlay
        axL = axes[r, 0]
        if len(t1) == len(v1) == len(t2) == len(v2) and np.allclose(t1, t2, atol=1e-12, rtol=0):
            x1, x2 = t1, t2
            axL.set_xlabel(r"time [$\sqrt{m/k}$]")
        else:
            x1 = np.arange(len(v1))
            x2 = np.arange(len(v2))
            axes[-1, 0].set_xlabel("index")

        axL.plot(x1, v1, linewidth=2, label=name1)
        axL.plot(x2, v2, linewidth=2, linestyle="--", label=name2)
        axL.set_ylabel(obs.replace("_", " ").title(), fontsize=11)
        axL.grid(True, linestyle="--", alpha=0.5)
        if r == 0:
            axL.legend(fontsize=9)

        # Right: |\delta| with interpolation of file2 to file1's x
        axR = axes[r, 1]
        x = t1 if len(t1) > 0 else np.arange(len(v1))
        if len(t2) > 1:
            v2_interp = np.interp(x, t2 if len(t2) > 0 else np.arange(len(v2)), v2,
                                  left=np.nan, right=np.nan)
        else:
            n = min(len(v1), len(v2))
            x = x[:n]
            v2_interp = v2[:n]
        v1_trim = v1[:len(x)]
        diff = np.abs(v2_interp - v1_trim)
        axR.plot(x, diff, linewidth=1.8)
        axR.axhline(0, lw=1, ls=":")
        axR.set_ylabel(f"|$\Delta$ {obs.replace('_', ' ').title()}|", fontsize=11)
        axR.grid(True, linestyle="--", alpha=0.5)
        axes[-1, 1].set_xlabel(r"time [$\sqrt{m/k}$]" if len(t1) == len(v1) else "index")

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if save_prefix:
        out = f"{save_prefix}_{grp_label.replace(' ', '')}.png"
        fig.savefig(out, dpi=150)
        print(f"saved figure: {out}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_groups(s):
    """
    parse --groups string. Default is "observables/A:Group A,observables/B:Group B".
    format: path[:label][,path[:label],...]
    """
    groups = {}
    for chunk in s.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" in chunk:
            path, label = chunk.split(":", 1)
        else:
            path = chunk
            label = path.split("/")[-1]  # A, B, etc.
        groups[path.strip()] = label.strip()
    return groups


def main():
    ap.add_argument("files", nargs=2, help="Two HDF5 files ...")
    ap.add_argument("--groups", default="observables/A:Group A,observables/B:Group B", ...)
    ap.add_argument("--observables", default="temperature,pressure,potential_energy,internal_energy", ...)
    ap.add_argument("--tol", type=float, default=1e-2, ...)
    ap.add_argument("--no-plot", action="store_true", ...)
    ap.add_argument("--save-prefix", default=None, ...)

    args = ap.parse_args()

    file1, file2 = args.files
    name1, name2 = Path(file1).name, Path(file2).name

    # sanity check (clearer errors in CI)
    for p in (file1, file2):
        if not Path(p).exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            return 5

    groups = parse_groups(args.groups)
    obs_try = [o.strip() for o in args.observables.split(",") if o.strip()]
    tol = args.tol

    any_compared = False
    any_mismatch = False

    with h5py.File(file1, "r") as f1, h5py.File(file2, "r") as f2:
        for grp_path, grp_label in groups.items():
            # collect common observables for this group
            common = []
            for obs in obs_try:
                t1, v1 = get_time_values(f1, grp_path, obs)
                t2, v2 = get_time_values(f2, grp_path, obs)
                if t1 is not None and v1 is not None and t2 is not None and v2 is not None:
                    common.append(obs)

            if not common:
                print(f"\n[{grp_label}] No common observables to compare.")
                continue

            print(f"\n--- Numeric check @ {grp_label} ({name1} vs {name2}, tol={tol:g}) ---")
            # compute and print numeric check
            series = {}
            for obs in common:
                t1, v1 = get_time_values(f1, grp_path, obs)
                t2, v2 = get_time_values(f2, grp_path, obs)
                series[obs] = (t1, v1, t2, v2)

                ok = print_numeric_check(name1, name2, grp_label, obs, t1, v1, t2, v2, tol)
                any_compared = True
                if not ok:
                    any_mismatch = True

            if not args.no_plot:
                plot_comparison(
                    name1, name2, grp_label, common, series,
                    show=(args.save_prefix is None), save_prefix=args.save_prefix
                )

    if not any_compared:
        print("ERROR: No observables were compared (no common data found).", file=sys.stderr)
        return 4 # nothing was compared
    if any_mismatch:
        # Make CTest fail cleanly
        return 2 # at least one observable faled the tolerance
    return 0 # all compared obs where within tolerance


if __name__ == "__main__":
    sys.exit(main())
