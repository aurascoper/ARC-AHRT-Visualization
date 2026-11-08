# Runtime physics data (intentionally not bundled)

Do not substitute the synthetic parser fixture under `tests/fixtures/` for a
physical XCOM table.

Obtain a compound table for `H2O` from the official NIST XCOM interface:

- https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
- https://physics.nist.gov/PhysRefData/Xcom/Text/chap4.html

Request the standard grid spanning 0.001 MeV through at least 0.6617 MeV and
include an additional energy of 0.6617 MeV if available. XCOM's UI documents
four significant figures for user-supplied energies. Retain the original raw
XCOM output in this directory, then normalize it to `xcom_water.csv` with the
exact columns below. Do not interpolate through absorption edges.

Required normalized CSV columns (units are part of the names):

```csv
energy_MeV,mu_coherent_cm2_g,mu_incoherent_cm2_g,mu_photoelectric_cm2_g,mu_pair_nuclear_cm2_g,mu_pair_electron_cm2_g,mu_total_with_coherent_cm2_g,mu_total_without_coherent_cm2_g
```

Create `xcom_water_manifest.json` beside the CSV. It must contain:

```json
{
  "source": "https://physics.nist.gov/PhysRefData/Xcom/ ...",
  "material": "H2O",
  "units": "cm2/g; energy_MeV",
  "retrieved_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "sha256": "SHA256_OF_xcom_water.csv",
  "raw_source_file": "original_xcom_output.txt",
  "notes": "NIST XCOM compound output; normalized without changing coefficient values"
}
```

The application verifies the CSV hash and partial-to-total coefficient closure.
The XCOM table is the runtime source of attenuation coefficients; the code
contains no baked-in XCOM coefficient values.

After normalization, generate the manifest hash (then inspect and correct the
material/source/raw-source fields if needed):

```bash
python scripts/write_manifest.py data/xcom_water.csv \
  --material H2O \
  --raw-source data/original_xcom_output.txt \
  --out data/xcom_water_manifest.json
```
