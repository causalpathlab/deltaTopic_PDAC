# Figures

## Figure 1. Overview of Delta Topic model (Yichen) between a week and a month

- Figure 1a. Schematic PGM

- Figure 1b. Scatter plot examples (spliced vs. unspliced) <- based on the best "delta" genes 

	- Only showing top genes out of a giant webpage that contains all the gene-level figures
	
	- A much better option: shiny app

- Figure 1c. The same scatter plot examples colored by Peng et al. cell types

- Figure 1d. Delta topic model provides higher resolution views

## Figure 2. Delta topic fitted to PDAC data -- make sense out of it (+ Pattie)

- Figure 2a. structure plot

- Figure 2b. rho matrix (shared between the unspliced and spliced)

- Figure 2c. delta matrix (spliced-specific)

- Figure 2d. known cell types for each topic

- Input: topic-specific delta/rho matrices

- Figure 2e. GSEA for each topic

## Figure 3. Topic-specific differential velocity (pseudo-bulk + Yongjin)

- Input: topic fractions (cell x topic matrix) **where in the astrocyte?**

- Figure 3a-c. Top examples of genes with different velocity between cancer vs. non-cancer for some topics

- Figure 3d. Counts.. vs. FDR cutoff?

## Figure 4. Deconvolution by dynamic states (asking if we can give new insight into existing bulk data + together)

- 4a. Fig regression (for each individual $i$):

$$y_{gi} \sim s_{i} \sum_{t} \delta_{gt} \pi_{ti}$$

$\pi_{ti}> 0$ and $\sum_{t} \pi_{ti} = 1$

- 4b. Show that the estimated dynamic topic fractions are relevant to clinical outcomes (US data)

- 4c. (Canadian data)

- 4d. (Australian data)


# Supplementary Figures

## Supp. Figure 1. UMAP of PCA is sensitive to the choice of number of PCs on the spliced count

- SFig 1a. UMAP with PC 8

- SFig 1b. UMAP with PC 16

- SFig 1c. UMAP with PC 32

- SFig 1d. UMAP with PC 64

## Supp. Figure 2. UMAP of PCA is sensitive to the choice of number of PCs on the combined data

- Sfig 2a. UMAP with PC 8

- Sfig 2b. UMAP with PC 16

- Sfig 2c. UMAP with PC 32

- Sfig 2d. UMAP with PC 64

## Supp. Figure 3. UMAP of deltaTopic is robust to the choice of number of topics

- Sfig 3a. UMAP with 8

- Sfig 3b. UMAP with 16

- Sfig 3c. UMAP with 32

- Sfig 3d. UMAP with 64


# Supplementary materials

- project githubs: (1) `deltaTopic` (2) `deltaTopic_PDAC`

- a full list of figure 3. topic-specific DEGs

- a full list of figure 4. deconvolved cell type fraction for all the ICGC PDAC samples

- `zenodo` for everything in tar ball

