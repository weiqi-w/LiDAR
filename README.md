# LiDAR

LiDAR infers spatial cell–cell communication by fitting ligand-specific distance-decay profiles to downstream transcriptional responses in receptor-positive spatial units. It estimates ligand exposure, spatial range and sample-level exposure–response coupling, and computes communication scores for spatial units and sender–receiver pairs.

## Installation

LiDAR requires R 4.3 or later. Install the dependencies and package from GitHub.

```r
install.packages(c("remotes", "BiocManager"))
repos <- BiocManager::repositories()
remotes::install_github("RausellLab/CelliD", repos = repos, upgrade = "never")
remotes::install_github("hfang-bristol/dnet", repos = repos, upgrade = "never")

remotes::install_github("weiqi-w/LiDAR", repos = repos, upgrade = "never")
library(LiDAR)
```

For a local installation, replace the last installation command with:

```r
remotes::install_local("LiDAR_0.1.0.tar.gz", repos = repos, upgrade = "never")
```

## Try the synthetic example

Run the full workflow on 100 spatial units and inspect the results for one ligand.

```r
source(system.file("examples", "synthetic.R", package = "LiDAR", mustWork = TRUE))
GetDecay(z)
PlotDecay(z, "g1")
PlotCoupling(z, "g1")
```

## Analyze a tissue section

Prepare a Seurat object (`seurat_obj`, default assay `SCT`, layer `data`) and a coordinate table (`coords`) with numeric `x` and `y` columns and spatial-unit IDs as row names. Use matching gene identifiers in the expression data and priors.

Download the target-by-ligand matrix and ligand–receptor network for your species from the NicheNet v2 dataset deposited by Robin Browaeys (2022). [NicheNet v2 priors on Zenodo](https://zenodo.org/records/7074291)

| Species | Ligand–target matrix | Ligand–receptor network |
| --- | --- | --- |
| Human | `ligand_target_matrix_nsga2r_final.rds` | `lr_network_human_21122021.rds` |
| Mouse | `ligand_target_matrix_nsga2r_final_mouse.rds` | `lr_network_mouse_21122021.rds` |

```r
lidardata <- CreateLiDAR(
  seurat_obj,
  ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds"),
  lr_network = readRDS("lr_network_human_21122021.rds"),
  coord = coords, coord_cols = c("x", "y")
)
lidardata <- ImputeExpr(lidardata)
lidardata <- BuildNetworks(lidardata)
lidardata <- InferResponse(lidardata)
lidardata <- BinLigandExpr(lidardata, bin = 50, n_bins = 10)
lidardata <- FitPropagation(lidardata)
lidardata <- ComputeScores(lidardata, mode = "aggregate")
lidardata <- ComputeScores(lidardata, mode = "pairwise")
```

- **Networks.** The default neighbor counts are 500 and 700. These defaults can be used for typical Visium sections with 1,000–4,000 spots.
- **Distance.** The example uses the original coordinate units. The first bin, `[0, 50)`, has weight 1, and the decay weight reaches zero at 500 coordinate units. Calling `BinLigandExpr(lidardata)` uses the smallest positive inter-unit distance as the bin width and reaches zero at ten times that distance.
- **Permutation testing.** To obtain `p.empirical`, use `FitPropagation(lidardata, permutation_p = TRUE)`. This enables 999 permutations, each with a refit of the decay model.
- **Scoring.** Both modes default to Pearson `p.value < 0.001` and `coupling > 0`. Set `pval = NULL` or `coupling = NULL` to disable the corresponding filter. Scores use their original scale by default.

### Plot one ligand

For a human section in which IL1B passes both filters:

```r
ligand <- "IL1B"
if (ligand %in% colnames(lidardata@Lidarscore@aggregated_score)) {
  lidardata_plot <- ComputeScores(lidardata, ligand_list = ligand, Log = TRUE)
  ExportPlots(
    lidardata_plot, ligand_list = ligand,
    features = c("LigandExpr", "Exposure", "Response", "Score"),
    ratio = 0.8, format = "pdf", out_dir = "IL1B_maps"
  )
}
```

Spatial maps show Z-scores by default. The LiDAR score panel displays the Z-score of `log10(1 + LiDAR score)`. Set `zscore = FALSE` to plot stored values.

## Results

| Result | Access |
| --- | --- |
| Decay parameters and exposure–response coupling | `GetDecay(lidardata)` |
| Exposure and response per spatial unit | `GetExposure(lidardata)`, `GetResponse(lidardata)` |
| LiDAR score per spatial unit and ligand | `lidardata@Lidarscore@aggregated_score` |
| Sender-by-receiver matrix | `lidardata@Lidarscore@pairwise_score[[ligand]]` |
| Per-ligand failures | `attr(lidardata, "failure_log")` |

See the [workflow guide](inst/doc/workflow.md) for model and plotting details.

## Testing

LiDAR 0.1.0 passed R CMD check on Windows 11 (R 4.3.0) and Rocky Linux 9.7 (R 4.3.3).

## License

MIT. See [LICENSE](LICENSE).
