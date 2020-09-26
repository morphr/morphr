<img src="https://github.com/morphr/morphr/blob/master/assets/images/logo.png" width="350"/>

# Morphr (c) 2020

[**Morphr**](https://github.com/morphr/morphr) implements standard techniques for the representation of shapes, calculation of mean templates, statistical analysis (clustering, principal component analysis (PCA) and estimating deformations to mean shapes). It accepts case studies for medical, anatomical, and paleontological shape data. The package also implements convenient graphical user interface built using R shiny, allowing researchers with limited coding experience to conduct shape analysis case studies easily. 

The idea behind **Morphr** is not only to provide tools for analysis of continuous shapes and curves, but also a set of integrated functionalities that can be used by morphologist with limited coding experience to conduct shape analysis case studies easily. The package also utilizes ideas and functions from elastic functional data analysis pacakge: [**fdasrvf**](https://cran.r-project.org/web/packages/fdasrvf/index.html).


### Installation

* Install [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
* Install latest development version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("morphr/morphr")
```

### Examples

For demonstration, we present a case study for paleontological shapes. For a dataset of hadrosaurid bone specimens, we show results for building a mean template from a collection of shapes, summarize patterns of shape variation in PCA, and show hierarchical clusters obtained from pairwise distances between the collection. 

* Reconstruct Shape
```{r}
path = 'Path for ucf or svg file list'
verify_shapes(path)
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/Input_Shape.png width=50% />

* Construct Mean shape
```{r}
X = main_closed(path)
qarray = list()
for(i in 1:length(X)){
  qarray[[i]] = curve_to_q(X[[i]])
}
All_mean_shape = find_mean_shape(qarray)
saveRDS(All_mean_shape, file = "All_mean_shape.rds")
qmean = All_mean_shape[[1]]
qmean_new = project_curve(qmean)
pmean = q_to_curve(qmean_new)
plot_curve(pmean,'r',filename = 'Mean Shape')
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/Mean_Shape.png width=50% />


* Calculate Geodesic Distance Heatmap
```{r}
color_file_path = 'Path for colorcode and filename'
all_geodesic = geodesic_distance_all(qarray)
saveRDS(all_geodesic, file = "all_geodsic.rds")
geo_dist = all_geodesic[[6]]
output_result[lower.tri(output_result, diag=FALSE)] <- unlist(geo_dist)
   output_result <- t(output_result)
   for(i in 1:(length(X)-1)){
     for (j in (i+1):length(X)){
        output_result[j,i] = output_result[i,j]
      }
   }
rc_name <- c()
Taxoncolorcodes <- readr::read_csv(color_file_path,col_names = FALSE)
filenames <- Taxoncolorcodes$X1
for(i in 1:length(X)){
    rc_name <- c(rc_name,filenames[i])
}
rownames(output_result) <- rc_name
colnames(output_result) <- rc_name
plotly::plot_ly(x = rc_name, y = rc_name ,z = output_result, type = "heatmap")
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/geo_distance.png width=50% />

* PCA Plot on Eigen Axis 1 and 2
```{r}
alpha_t_array = All_mean_shape[[3]]
PCA_plot(alpha_t_array, qmean, qarray, color_file_path)
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/PCA_Plot.png width=50% />

* PCA Variation Along Eigen Axis 1
```{r}
plot_pca_variation(alpha_t_array, qmean, eigdir = 1,qarray, dt = 0)
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/PCA_Variation.png width=50% />

* MDS Plot
```{r}
mdsplot(alpha_t_array, geo_dist, X, color_file_path)
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/MDS_Plot.png width=50% />

* Dendrogram
```{r}
plot_dendrogram(geo_dist, color_file_path)
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/Dendrogram.png width=50% />

* Deformation Field
```{r}
deformation_field_all(alpha_t_array, pmean, qmean,X)
```
<img src=https://github.com/morphr/morphr/raw/master/assets/images/Deformation_Field.png width=50% />
