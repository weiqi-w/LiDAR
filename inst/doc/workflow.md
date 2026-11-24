# Method and workflow guide

## Coordinates and bins

Explicit coord_cols takes priority. Otherwise the implementation selects a
recognized pair (X/Y, x/y, col/row or supported image-coordinate names).
Coordinates are aligned by spatial-unit names. No unit conversion is performed.
Check coordinate columns and units before interpreting physical distances.

With bin=50 and n_bins=10, bins are [0,50), [50,100), ..., [500,550).
Weights are evaluated at left edges 0,50,...,500, with B=500. The first
weight is one, the last zero; distances beyond support also contribute zero.
At exactly 50 the second bin applies. Self contribution is always retained.
With bin=NULL the minimum strictly positive inter-unit distance is used.

## Response inference

ImputeExpr replaces zeros only, using inverse-distance embedding neighbors.
BuildNetworks creates simple undirected graphs without duplicate edges or loops.
SelectSeeds uses the 99th percentile of strictly positive target-gene weights
for each ligand and retains ties. Cell nodes are not in the seed percentile.
InferResponse aligns dnet rows to graph vertex order and columns to seed order
when dnet omits names; it validates dimensions and any existing names.
Candidate selection uses binary receptor availability. An invalid candidate
is skipped; a ligand without a usable candidate is recorded as failed.

## Decay, exposure and ranges

For 0 <= b <= B, D(b) = (1 - (b/B)^n) / (1 + (b/(B*ki))^n).
No max-min normalization is applied. Exposure is the sum of source ligand
expression weighted by the bin decay weights. CalcExposure performs the
matrix multiplication. FitPropagation handles receiver selection and fitting;
SummarizeDecay summarizes a fit and already-calculated exposure.
CalcRanges computes r50 and r01 from the continuous curve rather than bin edges.
Finite valid nonconverged optimization results are accepted with a warning.
fit_timeout defaults to 300 seconds per fit and is cooperative rather than a
hard process kill. Invalid inputs and fit failures are recorded in failure_log.

## Permutation inference

permutation_p defaults to FALSE. If TRUE, responses are permuted across the
fixed potential receivers; expression, coordinates and receiver membership
remain fixed. Each permutation repeats the same bounded optimization, recording
the optimized Pearson correlation. Empirical P = (1 + number of null correlations
at least observed correlation) / (M + 1). Default M=999; optionally use 9999.
A failed permutation makes the empirical test incomplete and P=NA, rather than
silently changing the denominator. Null values and diagnostic information are
stored in attr(GetDecay(x), "permutation_results").

This is an optimization-aware random-label null, assuming exchangeability of
receiver responses. It does not preserve response spatial autocorrelation.
Report that assumption and consider spatial null models for biological inference.
The ordinary Pearson P value after optimization is descriptive; empirical P
values are reported separately and are not adjusted across ligands.

## Scores and plotting

Scores preserve maximum scaling times 10, original products and sums. Aggregate
scores are not forced to equal pairwise totals. Sender and receiver scores are
outgoing and incoming pairwise sums. See ComputeScores for explicit ligand
selection; both ordinary Pearson P < 0.001 and coupling > 0 apply to all selections. Set pval or coupling to NULL to disable the corresponding filter. Log=TRUE stores log10(1 + score) after calculating each score or total.

Spatial plots use Z-scores by default. Use zscore=FALSE for stored values. Different
features use independent color scales; shared_scales=TRUE shares limits for
the same feature across exported ligands. reverse_y defaults to FALSE.
Decay plots show the continuous curve. ExportPlots returns a manifest and log,
and does not overwrite existing files by default.

## Reproducible example

Run source(system.file('examples', 'synthetic.R', package='LiDAR')). The example
uses synthetic expression, coordinates and positive prior weights. No user
networks or biological data are distributed. For real analysis, replace these
inputs with a normalized Seurat object and priors read using LoadPriors.
