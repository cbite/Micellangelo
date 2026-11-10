# Micellangelo

This repository contains the code, trained checkpoints, evaluation notebooks and experimental protocols associated with [*Micellangelo: a Generative Model of Cell-Topography Interactions* (bioRxiv, v2)](https://www.biorxiv.org/content/10.1101/2025.11.24.690184v2). The preprint describes generation of fluorescence images of human dermal fibroblasts conditioned on surface topography, and completion of partially observed fluorescence images.

## Preprint materials

| To find | Start here |
| --- | --- |
| Train the models or generate images with released weights | [Model setup, checkpoints and commands](Model/README.md) |
| Recreate the quantitative plots from archived measurements | [Evaluation instructions](Model/README.md#evaluation) and [CellProfiler measurement archives](Model/evaluation/cp_output_archive/) |
| Revisit the qualitative and partial observation examples | [Qualitative notebook](Model/evaluation/paper_plots_qualitative.ipynb) and [partial observation notebook](Model/evaluation/paper_plots_partial_observation.ipynb) |
| Inspect the microscopy and preprocessing methods | [Wet lab protocols](Raw_data_and_preprocessing/Wet_lab_work_protocols/), [preprocessing SOP](Raw_data_and_preprocessing/SOPs/) and [per-chip CellProfiler protocols](Raw_data_and_preprocessing/cellprofiler_protocols_for_all_the_chips.zip) |
| Find the alpha interface described in the supplement | [GUI instructions](GUI/Test%20GUI%20README.txt) |

## Data access

- [Processed model inputs (Zenodo: 17598780)](https://doi.org/10.5281/zenodo.17598780): the single-cell MATLAB `.mat` files used with the model code. See [where to extract them](Model/README.md#data).
- [Raw microscopy screen and surface metadata (Zenodo: 18631579)](https://doi.org/10.5281/zenodo.18631579): TIFFs from the screen, released after preprint v2, with metadata for ten topographies. The record describes the channel names and insert identifiers. For raw data not in this deposit, the preprint gives the corresponding author as the point of contact.

Both Zenodo records list [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) for their datasets. See [LICENSE](LICENSE) for the repository's software license.

## Citation

Konshin, N., Minartz, K., de Boer, J., & Menkovski, V. (2026). *Micellangelo: a Generative Model of Cell-Topography Interactions*. bioRxiv, version 2. [doi:10.1101/2025.11.24.690184](https://doi.org/10.1101/2025.11.24.690184).

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata. When reusing the data, also cite the specific Zenodo record used. Record the Git commit, checkpoint filename and dataset DOI with your run.
