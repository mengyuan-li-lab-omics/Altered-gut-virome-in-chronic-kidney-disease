# Altered-gut-virome-in-chronic-kidney-disease

**Comprehensive Workflow for Single-Cell RNA Sequencing (scRNA-seq) Data Analysis**  

Single-cell RNA sequencing (scRNA-seq) has revolutionized our understanding of cellular heterogeneity, gene regulation, and developmental trajectories. Below, I outline my step-by-step analytical workflow, covering raw data processing, quality control, clustering, and advanced downstream analyses. This pipeline integrates established tools (e.g., Cell Ranger, Seurat, Scanpy) and custom scripts to ensure robust and reproducible results.  

### **1. Raw Data Preprocessing**  
The analysis begins with raw sequencing data (FASTQ files) generated from platforms like 10x Genomics. Key steps include:  
- **Alignment & Gene Counting**: Using `Cell Ranger` (v7.1.0), I align reads to a reference genome (e.g., GRCh38) and generate a feature-barcode matrix. Default parameters are adjusted for read depth (>50,000 reads/cell) and minimum gene detection (>500 genes/cell).  
- **Demultiplexing**: For multiplexed samples (e.g., CellPlex), I employ `cellranger multi` to assign cells to individual samples.  

### **2. Quality Control (QC) & Filtering**  
Low-quality cells and artifacts are removed to ensure data reliability:  
- **Metrics**: Cells with high mitochondrial gene ratios (>10%) or low unique molecular identifiers (UMIs) are filtered out using `Seurat` (R) or `Scanpy` (Python).  
- **Doublet Detection**: I use `Scrublet` or `DoubletFinder` to identify and remove doublets (threshold: doublet score >0.25).  
- **Normalization**: Data are normalized via `SCTransform` (Seurat) or `scanpy.pp.normalize_total` to correct for library size biases.  

### **3. Dimensionality Reduction & Clustering**  
To uncover cell populations:  
- **Feature Selection**: Highly variable genes (HVGs) are identified (`FindVariableFeatures` in Seurat or `scanpy.pp.highly_variable_genes`).  
- **PCA**: Linear dimensionality reduction is performed (top 50 PCs), followed by `RunUMAP` or `sc.tl.umap` for visualization.  
- **Clustering**: Graphs are constructed using `FindNeighbors` (Seurat) or `sc.pp.neighbors` (Scanpy), and clusters are resolved via Louvain/Leiden algorithms (resolution: 0.4–1.2).  

### **4. Cell Type Annotation**  
Clusters are annotated using:  
- **Reference-Based Methods**: `SingleR` (with Human Primary Cell Atlas) or `Azimuth` for automated labeling.  
- **Marker Genes**: Manual inspection of canonical markers (e.g., *CD3D* for T cells, *CD79A* for B cells) via `FindAllMarkers` (Seurat).  

### **5. Advanced Analyses**  
- **Trajectory Inference**: Pseudotime analysis is conducted with `Monocle3` or `PAGA` to model differentiation pathways.  
- **Cell-Cell Communication**: `CellPhoneDB` predicts ligand-receptor interactions between clusters.  
- **Differential Expression**: `DESeq2` or `MAST` identifies genes upregulated in specific conditions (FDR <0.05, log2FC >0.25).  

### **6. Data Integration & Batch Correction**  
For multi-sample studies:  
- **Harmony** or `scanpy.pp.bbknn` corrects batch effects while preserving biological variation.  
- **Integration Anchors** (Seurat) align datasets using reciprocal PCA.  

### **7. Visualization & Reproducibility**  
- **Plots**: UMAPs, violin plots, and heatmaps are generated with `ggplot2` (R) or `scanpy.pl` (Python).  
- **Code Management**: All scripts are version-controlled via GitHub, and dependencies are containerized using Docker.  

### **Conclusion**  
This workflow balances rigor and flexibility, leveraging open-source tools to dissect cellular diversity. Future directions include spatial transcriptomics integration and multi-omics approaches.  

---  
**Tools Used**: Cell Ranger, Seurat, Scanpy, Monocle3, SingleR, CellPhoneDB.  
