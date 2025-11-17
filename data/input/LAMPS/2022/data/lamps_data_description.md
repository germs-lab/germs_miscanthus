
# Dataset Description

This project's data consists of bacterial and arbuscular mycorrhizal fungal (AMF) microbiome composition as characterized by 16S and AMF sequencing, respectively. In addition, a suite of associated soil physical, chemical, and biological measurements are used. The purpose is to link microbial community composition with soil biochemical properties and structural stability.

## Data Overview

-   **Sequencing**: Two *phyloseq* objects containing tables of amplicon sequence variant (ASV) abundance and taxonomy from (1) 16S rRNA gene for bacteria/archaea, and (2) AMF-targeted primers.
-   **qPCR**: quantification of major microbial groups (bacteria, fungi, AMF) via qPCR assays.
-   **Biochemical analyses**: Soil organic carbon fractions, microbial biomass, extracellular polysaccharides, and carbohydrates .
-   **Physical properties**: Bulk density, porosity, texture, soil aggregation, and water holding characteristics.

## Variables

| Column                   | Description                                                   | Data type    | Units            |
|--------------|-------------------------------|--------------|--------------|
| `Sample_id`              | Unique identifier for each sample                             | categorical  |                  |
| `Site`                   | Field site identifier                                         | categorical  |                  |
| `Treatment`              | Vegetation growing in the plot where the sample was collected | categorical  |                  |
| `Replicate`              | Replicate number within site/treatment                        | integer      |                  |
| `pH`                     | Soil pH                                                       | quantitative |                  |
| `BD`                     | Bulk density                                                  | quantitative | g·cm-3           |
| `Total_Porosity`         | Soil porosity                                                 | quantitative | \%               |
| `EC`                     | Electrical conductivity                                       | quantitative | 𝜇S·cm-1          |
| `Clay`                   | Clay fraction                                                 | quantitative | \%               |
| `Sand`                   | Sand fraction                                                 | quantitative | \%               |
| `Silt`                   | Silt fraction                                                 | quantitative | \%               |
| `SOC`                    | Soil organic carbon                                           | quantitative | g·kg-1           |
| `SEOC`                   | Salt-extractable organic carbon                               | quantitative | mg·kg-1          |
| `MBC`                    | Microbial biomass carbon                                      | quantitative | mg·kg-1          |
| `WSA`                    | Water-stable aggregates                                       | quantitative | \%               |
| `MWD`                    | Mean weight diameter of soil aggregates                       | quantitative | mm               |
| `EPS`                    | Extracellular polysaccharides                                 | quantitative | µg·g-1 soil      |
| `Total_Carbs`            | Total carbohydrates                                           | quantitative | mg·g-1 soil      |
| `Mean_Pressure_5_15cm`   | Mean soil pressure at 5–15 cm depth                           | quantitative | kPa              |
| `Root_Density`           | Root density                                                  | quantitative | g·m-3            |
| `Mean_WC_33kPa`          | Mean water content at -33 kPa (field capacity)                | quantitative | g·g-1            |
| `16S_gene_copies`        | 16S rRNA gene copy number (qPCR)                              | quantitative | copies·g-1 soil  |
| `ITS_gene_copies`        | ITS region gene copy number (qPCR)                            | quantitative | copies·g-1 soil  |
| `AMF_gene_copies`        | AMF-targeted gene copy number (qPCR)                          | quantitative | copies·g-1 soil  |
| `Fucose`                 | Neutral sugar – fucose                                        | quantitative | 𝜇g·g-1 soil      |
| `Arabinose`              | Neutral sugar – arabinose                                     | quantitative | 𝜇g·g-1 soil      |
| `Rhamnose`               | Neutral sugar – rhamnose                                      | quantitative | 𝜇g·g-1 soil      |
| `Galactose`              | Neutral sugar – galactose                                     | quantitative | 𝜇g·g-1 soil      |
| `Weak_Acid_Glucose`      | Glucose fraction (weak acid hydrolysis)                       | quantitative | 𝜇g·g-1 soil      |
| `Strong_Acid_Glucose`    | Glucose fraction (strong acid hydrolysis)                     | quantitative | 𝜇g·g-1 soil      |
| `Total_Glucose`          | Total glucose content                                         | quantitative | 𝜇g·g-1 soil      |
| `Xylose`                 | Neutral sugar – xylose                                        | quantitative | 𝜇g·g-1 soil      |
| `Mannose`                | Neutral sugar – mannose                                       | quantitative | 𝜇g·g-1 soil      |
| `Plant_Derived_Carbos`   | Carbohydrates derived from plant sources                      | quantitative | 𝜇g·g-1 soil      |
| `Microbe_Derived_Carbos` | Carbohydrates derived from microbial sources                  | quantitative | 𝜇g·g-1 soil      |
| `Plant_Microbe_Ratio`    | Ratio of plant-derived to microbe-derived carbohydrates       | quantitative |                  |
| `MWHC_corrected`         | Maximum water holding capacity (corrected)                    | quantitative | g H~2~O·g-1 soil |
