# metawatchar

Research repository for a PhD project on **electromagnetic metamaterials at RF frequencies**.

**Author:** Watchara Amasiri, Department of Electrical and Computer Engineering, Thammasat University

> Status: early stage, work in progress. Content will be added as the research progresses.

## Approach

The work combines simulation, fabrication, and measurement:

1. **Design and simulate** structures with open-source full-wave electromagnetic solvers, using an HPC cluster for parameter sweeps and optimization.
2. **Fabricate** prototypes as printed circuit boards.
3. **Measure** prototypes with RF instruments (vector network analysis, and free-space measurement in an anechoic chamber).
4. **Compare** measured and simulated results.

## Repository layout

| Path | Contents |
|---|---|
| `sim/<case>/` | Simulation scripts, job files, and summarized results (e.g., S-parameters, figures) for each case |

Large raw outputs, such as 3-D field dumps, are not stored in this repository.

## Reproducibility

Each case folder will describe the solver and version, the input parameters, and how to regenerate its results.

## License

- **Code** (scripts, simulation models, job files) is licensed under the [MIT License](LICENSE).
- **Data and results** (measured and simulated data, figures, and documentation) are licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-CC-BY-4.0).

If you use this work, please cite the related publications (to be listed here) and credit Watchara Amasiri.
