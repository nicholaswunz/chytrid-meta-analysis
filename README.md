# Infection intensity predicts functional disruption
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)

This repository contains code and data needed to reproduce the article:

**Wu N. C.** (2023) Pathogen load predicts host functional disruption: A meta-analysis of an amphibian fungal panzootic. *Functional Ecology*, **37**, 900–914. DOI: [![DOI](https://zenodo.org/badge/DOI/10.1111/1365-2435.14245.svg)](https://doi.org/10.1111/1365-2435.14245)

**When using the data or code from this project, please cite it as:**

**Wu N. C.** (2022) nicholaswunz/chytrid-meta-analysis: Accepted version of paper data and code of manuscript: Pathogen load predicts host functional disruption: A meta-analysis of an amphibian fungal panzootic (Functional Ecology). *Zenodo*. DOI: [![DOI](https://zenodo.org/badge/437831315.svg)](https://zenodo.org/badge/latestdoi/437831315)

**Raw data**
- `trait_raw_data.csv` - Raw data for at the individual level used for the analysis.
- `trait_corr_data.csv` - Correlation data used for the analysis.

**Analysis workflow**
- [`supplementary_information.html`](https://nicholaswunz.github.io/chytrid-meta-analysis/supplementary_information.html) - Supplementary information which contains the *R* workflow for processing and analysing the raw data, creating figures, and supplementary material for statistical outcomes, additional figures, and descriptions from the main document.

**Files**
- `Wu_table_ref.txt` - Text file containing citation information from Table S1 in [Wu (2019)](https://espace.library.uq.edu.au/view/UQ:39c96d7).
- `naive_search_WoS_211028.txt` - Text file containing citation information from the näive Boolean search in [Web of Science](https://www.webofscience.com/wos/woscc/basic-search).

## Abstract
1. The progression of infectious disease depends on the intensity of and sensitivity to pathogen infection. Understanding commonalities in trait sensitivity to pathogen infection across studies through meta-analytic approaches can provide insight to the pathogenesis of infectious diseases. The globally devastating amphibian chytrid fungus, *Batrachochytrium dendrobatidis* (*Bd*), offers a good case system due to the widely available dataset on disruption to functional traits across species.  
2. Here, I systematically conducted a phylogenetically controlled meta-analysis to test how infection intensity affects different functional traits (e.g., behaviour, physiology, morphology, reproduction) and the survival in amphibians infected with *Bd*.  
3. There was a consistent effect of *Bd* infection on increasing energy metabolism, while other traits varied substantially. Skin integrity, hormone levels and osmoregulation were most sensitive to *Bd* infection (minimum *Bd* load ln 2.5 ZE), while higher minimum *Bd* loads were required to influence reproduction (ln 10.6 ZE). Mortality differed between life stages, where juvenile mortality was dependent on infection intensity and exposure duration, while adult mortality was depended on infection intensity only. Importantly, there were strong biases for studies on immune response, body condition, and survival, while traits such as locomotor capacity, energy metabolism, and cardiovascular traits were lacking.  
4. The influence of pathogen load on functional disruption can help inform pathogen thresholds before the onset of irreversible damage and mortality. Meta-analytic approaches can provide quantitative assessment across studies to reveal commonalities, differences, and biases of panzootic diseases, especially for understanding the ecological relevance of disease impact.  

**Keywords:** anuran, *Batrachochytrium dendrobatidis*, chytridiomycosis, emerging infectious diseases, meta-analysis, pathogen

## License
This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).
