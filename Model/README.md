# Micellangelo model and evaluation

These instructions cover the released materials for the [Micellangelo preprint](https://www.biorxiv.org/content/10.1101/2025.11.24.690184v2). Start in the repository root, then use the `Model` directory for the commands below.

## Setup

```sh
cd Model
python -m pip install -r requirements.txt
python -m pip install absl-py tqdm ipywidgets
```

For the evaluation notebooks, also install `jupyter`, `seaborn` and `scikit-learn`. The supplied [requirements](requirements.txt) pin PyTorch 2.4.1 and list the CUDA 12.4 wheel index. The preprint used a GPU for training.

## Data

Download `clean_database_for_Micellangelo_model.7z` from the [processed model input record](https://doi.org/10.5281/zenodo.17598780). Extract the archive so its `Chip_*` folders are directly under `Model/data/`, preserving the subfolders and filenames. For example, `Model/data/Chip_1130_1/QualityCheckPassed/` should contain the corresponding `.mat` files. Keep [normalization_constants.pickle](normalization_constants.pickle) in `Model/`.

The [dataset loader](topodataset.py) reads the `I` array in each MATLAB file. The stored channel order is topography, DAPI, actin, YAP. The loader reorders it to topography, YAP, actin, DAPI before normalization and separates the topography image for conditioning. The [raw TIFF archive](https://doi.org/10.5281/zenodo.18631579) is a separate resource for inspecting or repeating preprocessing.

## Released checkpoints

| Use | Checkpoint | SHA-256 |
| --- | --- | --- |
| Topography conditioned generation (OTCFM) | [Step 145000](model_weights/otcfm_topo_weights_step_145000.pt) | `419eef699e2ff05f90f14d183064ae23f95341e4f97a2dc4b1a4be525519757c` |
| Completion from partial fluorescence observations (DDC) | [Step 345000](model_weights/icfm_topo_weights_step_345000.pt) | `24c38065d2b526fb74d57b0335113f699ec1d6c1eee8bcbdcbfe6d1fd992dc7c` |

To generate TIFFs with the released OTCFM EMA weights and the topographies in `data/`, run:

```sh
python generate_fake_tiffs.py --num_samples 4 --cfg_strength 1
```

This small example writes TIFFs under `generated_images_cfg_1.0/`. For DDC inference, use the [partial observation notebook](evaluation/paper_plots_partial_observation.ipynb).

## Training

[topochip-train.py](topochip-train.py) supports the OTCFM and DDC model configurations. The commands below start new training runs. For evaluation using the released models, use the checkpoints above.

```sh
python topochip-train.py --output_dir ./otcfm_results/
python topochip-train.py --output_dir ./ddc_results/ --model icfm --use_ddc --total_steps 345000
```

See the training script for flags and defaults. These commands identify the two configurations; they do not establish the exact training provenance of the released weights.

## Evaluation

To revisit the measurement-based plots without regenerating images, extract the three ZIP files in [evaluation/cp_output_archive/](evaluation/cp_output_archive/) in place. Retain their top-level folders: `cp_output_allreal/`, `cp_output_allfake/` and `cp_workdirs_cfg_strength/`. Open [paper_plots_quantitative.ipynb](evaluation/paper_plots_quantitative.ipynb) with `Model/evaluation/` as the notebook working directory. The cell-location overlays also need real background TIFFs; point `path_real_imgs` at the real TIFF export described below.

For the image examples, use [paper_plots_qualitative.ipynb](evaluation/paper_plots_qualitative.ipynb) and [paper_plots_partial_observation.ipynb](evaluation/paper_plots_partial_observation.ipynb) from the same directory. They load the released checkpoints and model inputs. Before launching Jupyter for these notebooks or for TIFF conversion, set `PYTHONPATH` to the absolute path of `Model/` so the notebooks can import `utils`, `unet` and `topodataset`.

To recompute measurements, convert the processed real inputs to TIFFs with [real_imgs_to_tiff.ipynb](evaluation/real_imgs_to_tiff.ipynb). Generate synthetic TIFFs as above, then apply [cp-pipeline-main.cpproj](evaluation/cp-pipeline-main.cpproj) in CellProfiler separately to the real and generated image sets. The preprint reports using CellProfiler 4.2.5 for segmentation. The [per-chip preprocessing protocols](../Raw_data_and_preprocessing/cellprofiler_protocols_for_all_the_chips.zip) are separate from this evaluation pipeline.
