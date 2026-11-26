# Decorrelation Before Variable Selection: Left and Right Transforms and Their Trade-Off

R code that reproduces the examples, tables and figure in

Mahdi Nouraie and Samuel Muller, *Decorrelation Before Variable Selection: Left and Right Transforms and Their Trade-Off*. Preprint: [arXiv:2611.20420](https://arxiv.org/abs/2611.20420).

The paper classifies decorrelating transforms for variable selection (Puffer, DECO, HOLP, ZCA, PCA, Cholesky, weighted Lasso) as left or right transforms, and shows what each one gives up.

## Files

| File | Reproduces |
|---|---|
| `reproduce_paper.R` | The (T1, N, T2) example, the properties table (Table 3) and the checks of Proposition 3 |
| `simulation.R` | Table 1 and Figure 1 (grouped designs A and B, n = 100, p = 50, 500 data sets) |
| `realdata.R` | Table 2 (Boston housing predictors with a simulated response, 500 data sets) |

## Requirements

R, with no extra packages. `realdata.R` uses the `MASS` package, which ships with R (it provides the Boston data).

## Usage

```
Rscript reproduce_paper.R
Rscript simulation.R     # also writes fig_simulation.pdf
Rscript realdata.R
```

Each script sets its random seed, so the output matches the paper.

## Citation

If you use this code, please cite the paper.

## License

MIT
