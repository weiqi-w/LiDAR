# LiDAR

**Ligand, Diffusion, Activity and Receptor**

A computational framework that uses inferred downstream transcriptional
responses to constrain models of ligand propagation through tissue.
LiDAR characterizes ligand-specific propagation patterns and quantifies
the spatial correspondence between modeled ligand exposure and downstream
transcriptional response in potential receiver cells. The framework further
enables the profiling of ligand-mediated cell-cell communication at the
original spatial resolution of the input data.

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("weiqi-w/LiDAR")
```

Alternatively, install from source by setting the working directory to
the package source directory and running:

```r
remotes::install_local(".", dependencies = NA, upgrade = "never")
```

R >= 4.1.0 is required. The full workflow additionally requires
`Seurat`, `CelliD` and `dnet`.

## Usage

**Basic usage**

*Load packages*

```r
library(LiDAR)
```

*Load a spatial dataset*

Provide a Seurat object containing normalized expression data and
spatial coordinates.

```r
seu <- readRDS("path/to/spatial_data.rds")
```

*Load reference databases*

Place the species-matched ligand-target matrix and ligand-receptor
network in a local directory.

Database download links: <!-- To be added -->

```r
db <- load_nichenet_data(
  species = "human",
  data_dir = "path/to/database/human"
)
```

The directory should contain `ligand_target_matrix.rds` and
`lr_network.rds`.

*Prepare the LiDAR object*

Specify the assay containing normalized expression data.

```r
set.seed(123)

x <- PrepLidarData(
  seu,
  assay = "SCT",
  ligand_target_matrix = db$ligand_target_matrix,
  lr_network = db$lr_network
)
```

*Infer downstream transcriptional responses*

```r
x <- ImputeExp(x, nn = 30)
x <- CreateNetworks(x, c(500, 700))
x <- InferLigandAct(x, parallel = FALSE)
```

*Model ligand propagation and quantify spatial coupling*

```r
x <- CalcLigandDist(
  x,
  dist_threshold = 10,
  include_self = TRUE
)

x <- OptLigandDiff(x)
```

*Calculate communication scores*

Calculate aggregate scores:

```r
x <- CompLidarScore(x, mode = "aggregate")
```

Calculate pairwise communication scores and associated sender and
receiver scores for a ligand of interest:

```r
x <- CompLidarScore(
  x,
  mode = "pairwise",
  ligand_list = c("YOUR_LIGAND")
)
```

Replace `YOUR_LIGAND` with a ligand available in the fitted results.
Pairwise score matrices can require substantial memory.

*Save results*

```r
saveRDS(x, "lidardata.rds")
```

## Documentation

<!-- Documentation and tutorial links to be added -->

Function documentation is available within R:

```r
?PrepLidarData
?InferLigandAct
?OptLigandDiff
?CompLidarScore
```

## Development status

LiDAR is currently a development release. 

See [VALIDATION.md](VALIDATION.md) for details.

## References

<!-- Manuscript citation and DOI to be added -->

## Contact

Weiqi Wang  
Wangwq22@mails.tsinghua.edu.cn

## License

MIT License.