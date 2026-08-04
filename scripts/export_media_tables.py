"""Export the media definitions to the formats other tools consume.

Generates, from the definitions in ``tools/media.py``:

* ``data/media/kbase_tsvs/<name>_media.tsv`` -- one per medium, for upload to
  KBase
* ``data/media/no_c_media_database.tsv`` -- the CarveMe media database for the
  two carbon-free base media (MBM and L1)

This was the second half of ``test/test_files/media/define_cobrapy_media.py``.
It lives in ``scripts/`` because it generates artifacts rather than changing the
model, and it is separate from the definitions because it needs a local clone of
the ModelSEED database while the definitions need nothing at all.

Requires ModelSEED's ``compounds.json``::

    export MODELSEED_COMPOUNDS=/path/to/ModelSEEDDatabase/Biochemistry/compounds.json
    python scripts/export_media_tables.py

Note this is *not* run in CI, both because of that external dependency and
because its outputs are committed rather than regenerated per build.
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.media import (  # noqa: E402
    MEDIA,
    convert_aliases_to_dict,
    load_modelseed_compounds,
    modelseed_id_from_exchange,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MEDIA_DATA_DIR = os.path.join(REPO_ROOT, "data", "media")
KBASE_TSV_DIR = os.path.join(MEDIA_DATA_DIR, "kbase_tsvs")

#: The carbon-free base media exported to the CarveMe database, with the
#: descriptions CarveMe shows.
CARVEME_BASE_MEDIA = {
    "mbm": "Minimal Basal Medium (Moran Lab)",
    "l1": "L1 Minimal Medium",
}


def write_media_tsv(media_dict, media_name, modelseed_db, out_dir=KBASE_TSV_DIR):
    """Write one medium as a KBase-style TSV."""
    media_df = pd.DataFrame.from_dict(media_dict, orient="index", columns=["minFlux"])
    # Fix the names of the compounds
    media_df.index = media_df.index.str.replace("EX_", "").str.replace("_e0", "")
    media_df.index.name = "compounds"
    # Make the min flux negative
    media_df["minFlux"] = -1 * media_df["minFlux"]
    # Set the max flux for everything to be 1000
    media_df["maxFlux"] = 1000
    # Add the names of the compounds
    media_df["name"] = media_df.index.map(lambda x: modelseed_db[x]["name"])
    # Add the formula of the compounds
    media_df["formula"] = media_df.index.map(lambda x: modelseed_db[x]["formula"])
    # Set the concentration to be something?
    # TODO: Does the concentration matter?
    media_df["concentration"] = 1
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, media_name + "_media.tsv")
    media_df.to_csv(path, sep="\t")
    return path


def make_carveme_media(media_dict, media_id, media_name, modelseed_db, media_db):
    """Append one medium to the CarveMe media database frame."""
    for ex_rxn, _min_flux in media_dict.items():
        met = modelseed_id_from_exchange(ex_rxn)
        name = modelseed_db[met]["name"]
        aliases = convert_aliases_to_dict(modelseed_db[met]["aliases"])
        if "BiGG" not in aliases:
            print(f"No BiGG ID for {name}")
            bigg_to_use = met
        else:
            bigg_id = aliases["BiGG"]
            if len(bigg_id) == 0:
                print(f"No BiGG ID for {name}")
                bigg_to_use = met
            if len(bigg_id) > 1:
                print(f"Multiple BiGG IDs for {name}: {bigg_id}")
            bigg_to_use = bigg_id[0]
        media_db = pd.concat(
            [
                media_db,
                pd.DataFrame(
                    {
                        "medium": media_id,
                        "description": media_name,
                        "compound": bigg_to_use,
                        "name": modelseed_db[met]["name"],
                    },
                    index=[0],
                ),
            ],
            ignore_index=True,
        )
    return media_db


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--modelseed-compounds",
        default=None,
        help="path to ModelSEED Biochemistry/compounds.json "
        "(defaults to $MODELSEED_COMPOUNDS)",
    )
    args = parser.parse_args(argv)

    modelseed_db = load_modelseed_compounds(args.modelseed_compounds)

    # Write the media TSV files
    for media_name, media_dict in MEDIA.items():
        write_media_tsv(media_dict, media_name, modelseed_db)
    print(f"Wrote {len(MEDIA)} KBase media TSVs to {os.path.relpath(KBASE_TSV_DIR, REPO_ROOT)}")

    # Make the CarveMe media database for the carbon-free base media
    media_db = pd.DataFrame(columns=["medium", "description", "compound", "name"])
    for media_id, description in CARVEME_BASE_MEDIA.items():
        media_db = make_carveme_media(
            MEDIA[media_id], media_id, description, modelseed_db, media_db
        )

    out = os.path.join(MEDIA_DATA_DIR, "no_c_media_database.tsv")
    media_db.to_csv(out, sep="\t", index=False)
    print(f"Wrote {os.path.relpath(out, REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
