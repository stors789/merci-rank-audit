# Tested reproduction environment

The accepted B1 split21 reconstruction used Python 3.12.3 and the exact package
versions in `constraints-tested-py312.txt`. From the extracted source distribution:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -c reproduce/constraints-tested-py312.txt '.[test,reproduce]'
```

These pins record that tested Python 3.12 environment. The package's supported
version ranges remain in `pyproject.toml`; the hosted CI environment is recorded
separately in its run logs.

The full B1 reconstruction used R 4.4.1 and an explicitly supplied user library:

| R package | Tested version |
|---|---|
| Matrix | 1.7.4 |
| data.table | 1.17.8 |
| jsonlite | 2.0.0 |
| e1071 | 1.7.17 |
| proxy | 0.4.29 |

Some local R packages reported being built under R 4.4.3; the documented R 4.4.1
run nevertheless produced byte-identical split21 outputs. The local Python
environment used cached dependency wheels; pypdf and ReportLab were repacked
locally from installed distributions for that isolated test and are not
redistributed. Hosted CI obtains its dependencies independently.

For the full-input run, set numerical-library threads to one and supply a new
output directory, the four original processed files and the separately obtained
pinned MERCI checkout as described in `INPUTS.md`. The recorded run used nice 10,
one CPU and a 4 GiB address-space ceiling; its peak RSS was about 2.62 GiB. These
settings and the input hashes accompany the validation record; wall time depends
on hardware and cache state. No original biological inputs or upstream source are
included in this source distribution.
