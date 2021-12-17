# General settings

This workflow needs 3 databases already installed before running it:

- Eggnog: https://github.com/eggnogdb/eggnog-mapper/wiki
- Uniclust30 (`hh-suite3`) : https://github.com/soedinglab/hh-suite
- PDB70 (`hh-suite3`) : https://github.com/soedinglab/hh-suite

Once installed, use `config/config.yaml` to set their paths and all the configuration values.

## Software Requirements

This worflow doesn't include definition of environments, but inside a env with `egnogg-mapper=2.1.6` (probably newer versions will work), install these packages if they aren't present yet:

```sh
# way faster than conda
mamba install snakemake=6.12.3 mmseqs hhsuite bioservices biopython pandas seaborn matplotlib
pip install solrq
```

Also, for the final step of this workflow, a script from *MicrobeAnnotator*, `ko_mapper.
py`, is used, which depends on some data installed by that package. **There is no need to setup the databases for this package**, but it still needs to be installed:

```
conda config --add channels cruizperez
conda create -n microbeannotator python=3.7 pip microbeannotator=2.0.5
```
See <https://github.com/cruizperez/MicrobeAnnotator> for updated instructions

Then, set the path to `ko_mapper.py` in the config file (use `find` in the `env` path).


## `rules/` and `scripts/`

This code it supposed to work as it is, but you are always free to modify it (specially if you known what you are doing.)
