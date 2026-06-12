# Anonymous Minimal Replica for Coreness Estimation

This directory is a minimal anonymous package for reproducing the core implementation used in our coreness experiments.

Included algorithms:
- `KcoreD`
- `KcoreH`
- `KcoreHT`
- `KcoreHA`

Everything unrelated to these four algorithms has been removed. Only one example dataset and one example figure are kept.

## Contents

- `graph.*`, `utility.*`: graph loading and shared helpers
- `single_coreness_function.h`, `kcored/`: distributed peeling baseline (`KcoreD`)
- `algorithms.*`: minimal implementations of `KcoreH`, `KcoreHT`, `KcoreHA`, metrics, and the four-algorithm runner
- `graphs/facebook/`: the only bundled example graph and its 80-way partition files
- `scripts/run_facebook_demo.sh`: end-to-end demo
- `scripts/plot_facebook_metrics.py`: generates a single Facebook summary figure

## Build

From this directory:

```bash
make
```

This produces:

```bash
./counting
```

## Minimal Run

Run the bundled Facebook example:

```bash
bash scripts/run_facebook_demo.sh
```

This command:

1. builds the binary,
2. runs `KcoreD`, `KcoreH`, `KcoreHT`, and `KcoreHA` on Facebook with `epsilon = 1`,
3. writes a CSV summary, and
4. generates one figure.

Outputs:

- `results/facebook_run.log`
- `results/facebook_metrics.csv`
- `figures/facebook_metrics.png`
- `figures/facebook_metrics.pdf`

## Direct Command

If you want to run the binary manually:

```bash
./counting 1.0 graphs/facebook/facebook.txt graphs/facebook/facebook_partitioned_80 results/facebook_metrics.csv
```

Arguments:

1. privacy budget `epsilon`
2. graph path
3. partition directory for `KcoreD`
4. optional CSV output path

## Figure

Only one figure is kept in this anonymous replica:

- `figures/facebook_metrics.pdf`

It reports Facebook overall comparison across the four algorithms using:

- Mean Factor
- MAE
- RMSE
- Runtime
