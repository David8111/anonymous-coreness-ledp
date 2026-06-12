#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


DISPLAY_NAME = {
    "kcored": "KcoreD",
    "kcoreH": "KcoreH",
    "kcoreHT": "KcoreHT",
    "kcoreHA": "KcoreHA",
}

ORDER = ["kcored", "kcoreH", "kcoreHT", "kcoreHA"]
COLORS = {
    "kcored": "#f2a3a0",
    "kcoreH": "#bddb8a",
    "kcoreHT": "#a9c9df",
    "kcoreHA": "#cbb7d8",
}


def load_rows(path: Path):
    rows = {}
    with path.open() as fin:
        reader = csv.DictReader(fin)
        for row in reader:
            rows[row["algo"]] = row
    return rows


def draw_metric(ax, rows, field, title):
    xs = range(len(ORDER))
    values = [float(rows[algo][field]) for algo in ORDER]
    colors = [COLORS[algo] for algo in ORDER]
    labels = [DISPLAY_NAME[algo] for algo in ORDER]

    ax.bar(xs, values, color=colors, edgecolor="#444444", linewidth=1.0)
    ax.set_title(title)
    ax.set_xticks(list(xs), labels, rotation=20, ha="right")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)


def main():
    parser = argparse.ArgumentParser(description="Plot a single Facebook metrics figure.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--png", required=True)
    parser.add_argument("--pdf", required=True)
    args = parser.parse_args()

    rows = load_rows(Path(args.input))

    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "axes.titlesize": 18,
        "axes.labelsize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
    })

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    draw_metric(axes[0], rows, "mean_factor", "Mean Factor")
    draw_metric(axes[1], rows, "mae", "MAE")
    draw_metric(axes[2], rows, "rmse", "RMSE")
    draw_metric(axes[3], rows, "time_s", "Time (s)")

    fig.suptitle("Facebook overall comparison", y=1.02, fontsize=20)
    fig.tight_layout()

    Path(args.png).parent.mkdir(parents=True, exist_ok=True)
    Path(args.pdf).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.png, dpi=300, bbox_inches="tight")
    fig.savefig(args.pdf, bbox_inches="tight")


if __name__ == "__main__":
    main()
