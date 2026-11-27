# Figures

`*.pdf` are the figures of the paper (Figs 1 to 7) and of the Supplementary Information (Supplementary Figs 1 to 8), as vector files. `source_data/` holds the numbers behind every panel as comma-separated files. `scripts/` holds the scripts that draw them, and `data/` the six small result files (`*.pkl`, 45 KB in all) that some of them read.

Every figure regenerates from the repository alone. From `figures/scripts/` run `python3 fig1_departures.py`, `fig2_buildup.py`, `fig3_plates.py`, `fig4_dose.py`, `fig5_errors.py`, `fig6_measured.py`, `fig7_map.py` and `figS1_gallery.py` to `figS8_exp.py`. Each script rewrites its PDF, a 600 dpi PNG and its source-data files in `figures/`. The scripts import the model and the simulator from `code/`, and read the experimental comparison from `code/literature_checks/experimental_comparison.csv` and the mesh test from `code/literature_checks/ripple_mesh_test.json`. Set `PERC_DATA` to use result files from another folder. The configuration images (Figs 1 and 2, Supplementary Figs 1 and 2) run the simulator and take under ten seconds each. Python 3.13 with the packages in `../requirements.txt`, plus Matplotlib and Pillow.

All 28 source-data files were regenerated from a clean copy of the repository on 27 November 2026 and are byte-identical to the ones in this folder. The drawings differ from the published PDFs only where a font is missing on the machine.

The PNG files the scripts write are not kept in the repository (`.gitignore`).
