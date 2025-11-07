
###count number of SCOOP peptides for each MIK
#!/usr/bin/env python3
import re
import pandas as pd

# ===== Paths (hard-coded) =====
mik2_path = "/Users/vishnu/v_files/GT/sem1/Research_sem1/sd03.txt"
scoop_path = "/Users/vishnu/v_files/GT/sem1/Research_sem1/sd01.xlsx"
out_csv    = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/scoop_counts_per_species.csv"

def header_to_key(hdr: str) -> str | None:
    """
    >A_thaliana_AT4G08850                -> Athaliana
    >T_cacao_Thecc08G107800___OUTGROUP   -> Tcacao
    """
    m = re.match(r">([^\s]+)", hdr.strip())
    if not m:
        return None
    token = m.group(1)
    parts = token.split("_")
    if len(parts) < 2:
        return None
    genus, species = parts[0], parts[1]
    return genus[0].upper() + species.lower()   # <-- critical fix

# ---- collect normalized MIK2 species keys from sd03 ----
mik2_species = set()
with open(mik2_path) as f:
    for line in f:
        if line.startswith(">"):
            k = header_to_key(line)
            if k:
                mik2_species.add(k)

print(f"Total MIK2 species (normalized): {len(mik2_species)}")

# ---- load sd01, skip the title row ----
scoop_df = pd.read_excel(scoop_path, skiprows=1)

# detect species column
species_col = next((c for c in scoop_df.columns if "species" in str(c).lower()), None)
if species_col is None:
    raise SystemExit("Could not find a 'species' column in sd01.xlsx")

# normalize species tokens the same way (strip spaces)
scoop_df[species_col] = scoop_df[species_col].astype(str).str.strip()

# ---- count matches ----
matches = (
    scoop_df[scoop_df[species_col].isin(mik2_species)]
    .groupby(species_col)
    .size()
    .reset_index(name="scoop_count")
    .sort_values("scoop_count", ascending=False)
)

missing = sorted(mik2_species - set(matches[species_col]))

print("\n=== SCOOP Counts per MIK2 species (overlap) ===")
print("None" if matches.empty else matches.to_string(index=False))

print("\n=== MIK2 species with NO corresponding SCOOPs in sd01 ===")
print("None" if not missing else missing)

matches.to_csv(out_csv, index=False)
print(f"\nWrote: {out_csv}")