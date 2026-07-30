"""Remove entities from the MIT1002 model and record them, in one step.

The point of this module is that removing a reaction from the model and logging
it in ``data/deprecatedIdentifiers/`` are the same operation. Doing them
separately is how a deprecated-identifier list rots.

Typical use, from the repo root::

    from tools.deprecate import deprecate_reactions

    deprecate_reactions(
        ["rxn00196_c0"],
        reason="no_genomic_evidence",
        notes="no candidate gene in the MIT1002 genome",
    )

There is deliberately no ``pr`` argument in that example. You remove things on
your branch before the pull request exists, so leave it blank and let CI stamp
it in -- see :func:`stamp_pr_number`.

or from the command line::

    python -m tools.deprecate reaction rxn00196_c0 \
        --reason no_genomic_evidence --dry-run

Both entry points cascade: metabolites and genes left with no reactions after a
removal are cleaned up too, and the orphaned metabolites are logged with
``reason="orphaned"``. That cascade is not a nicety -- ``test/test_sbml.py``
fails on isolated metabolites and genes, so a removal that skips it breaks CI.

The model is edited with libSBML rather than COBRApy so that formatting,
annotations and anything else COBRApy does not round-trip are preserved.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from typing import Iterable, Sequence

try:
    import libsbml
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "python-libsbml is required. It ships with cobra; try `pip install cobra`."
    ) from exc


# --------------------------------------------------------------------------
# Layout and schema
# --------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(REPO_ROOT, "model.xml")
DEPRECATED_DIR = os.path.join(REPO_ROOT, "data", "deprecatedIdentifiers")
REACTIONS_TSV = os.path.join(DEPRECATED_DIR, "deprecatedReactions.tsv")
METABOLITES_TSV = os.path.join(DEPRECATED_DIR, "deprecatedMetabolites.tsv")

COLUMNS = ["id", "name", "reason", "replaced_by", "pr", "date", "notes"]

#: Closed vocabulary for the ``reason`` column. Keep in sync with the table in
#: ``data/deprecatedIdentifiers/README.md``; ``test_deprecated.py`` enforces it.
REASONS = (
    "duplicate",
    "id_changed",
    "no_genomic_evidence",
    "no_experimental_evidence",
    "mass_imbalance",
    "infeasible_cycle",
    "dead_end",
    "orphaned",
    "erroneous_annotation",
    "out_of_scope",
    "unknown",
)

#: ``reason`` values that require ``replaced_by`` to be filled in, because the
#: claim they make is inherently relative to something else.
REASONS_REQUIRING_REPLACEMENT = ("duplicate", "id_changed")


class DeprecationError(RuntimeError):
    """Raised when a removal would leave the model or the TSVs inconsistent."""


@dataclass
class DeprecationRecord:
    """One row of a deprecated-identifier TSV."""

    id: str
    name: str = ""
    reason: str = "unknown"
    replaced_by: str = ""
    pr: str = ""
    date: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        self.id = strip_sbml_prefix(self.id)
        if self.reason not in REASONS:
            raise DeprecationError(
                f"reason {self.reason!r} is not in the controlled vocabulary. "
                f"Allowed: {', '.join(REASONS)}"
            )
        if not self.date:
            self.date = _dt.date.today().isoformat()
        if self.pr and not self.pr.startswith("#"):
            self.pr = "#" + self.pr.lstrip("#")
        # A single line only; the PR is the place for prose.
        self.notes = " ".join(self.notes.split())

    def as_row(self) -> dict:
        return {k: v for k, v in asdict(self).items() if k in COLUMNS}


# --------------------------------------------------------------------------
# Identifier helpers
# --------------------------------------------------------------------------

_PREFIX_RE = re.compile(r"^(R_|M_|G_)")


def strip_sbml_prefix(identifier: str) -> str:
    """``R_rxn00196_c0`` -> ``rxn00196_c0``. Idempotent.

    The TSVs store bare identifiers because that is what COBRApy reports, and
    matching what a user sees in ``model.reactions`` is the whole point.
    """
    return _PREFIX_RE.sub("", identifier)


def unescape_sbml_id(identifier: str) -> str:
    """Undo SBML character escaping, e.g. ``WP_012__46__1`` -> ``WP_012.1``.

    SBML ids may only contain word characters, so COBRApy encodes anything else
    as ``__<codepoint>__``. Gene locus tags contain dots, so they come back
    mangled; decode them before showing them to a human.
    """
    return re.sub(r"__(\d+)__", lambda m: chr(int(m.group(1))), identifier)


def _sbml_reaction_id(identifier: str) -> str:
    bare = strip_sbml_prefix(identifier)
    return bare if bare.startswith("R_") else "R_" + bare


def _sbml_species_id(identifier: str) -> str:
    bare = strip_sbml_prefix(identifier)
    return bare if bare.startswith("M_") else "M_" + bare


# --------------------------------------------------------------------------
# TSV I/O
# --------------------------------------------------------------------------


def read_records(path: str) -> list[DeprecationRecord]:
    """Read a deprecated-identifier TSV. Missing file reads as empty."""
    if not os.path.exists(path):
        return []
    out = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row.get("id"):
                continue
            out.append(
                DeprecationRecord(**{k: (row.get(k) or "") for k in COLUMNS})
            )
    return out


def write_records(path: str, records: Iterable[DeprecationRecord]) -> None:
    """Write records sorted by id, so the file diffs cleanly in review."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows = sorted((r.as_row() for r in records), key=lambda r: r["id"])
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def append_records(path: str, new: Sequence[DeprecationRecord]) -> list[str]:
    """Merge ``new`` into the TSV at ``path``. Returns the ids actually added.

    Re-deprecating an id that is already listed is a no-op rather than an
    error, so the helper stays safe to re-run.
    """
    existing = read_records(path)
    known = {r.id for r in existing}
    added = [r for r in new if r.id not in known]
    if added:
        write_records(path, existing + added)
    return [r.id for r in added]


# --------------------------------------------------------------------------
# Model surgery
# --------------------------------------------------------------------------


def _load(model_path: str):
    doc = libsbml.readSBMLFromFile(model_path)
    if doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR):
        raise DeprecationError(
            f"libSBML could not read {model_path}: {doc.getErrorLog().toString()}"
        )
    model = doc.getModel()
    if model is None:
        raise DeprecationError(f"no <model> element in {model_path}")
    return doc, model


def _reaction_species(reaction) -> set[str]:
    ids = set()
    for getter, count in (
        (reaction.getReactant, reaction.getNumReactants()),
        (reaction.getProduct, reaction.getNumProducts()),
        (reaction.getModifier, reaction.getNumModifiers()),
    ):
        for i in range(count):
            ids.add(getter(i).getSpecies())
    return ids


def _gene_products_of(reaction) -> set[str]:
    """Gene product ids referenced by a reaction's FBC gene association."""
    plugin = reaction.getPlugin("fbc")
    if plugin is None:
        return set()
    assoc = plugin.getGeneProductAssociation()
    if assoc is None or assoc.getAssociation() is None:
        return set()
    found = set()

    def walk(node):
        if node is None:
            return
        if isinstance(node, libsbml.GeneProductRef):
            found.add(node.getGeneProduct())
            return
        for i in range(getattr(node, "getNumAssociations", lambda: 0)()):
            walk(node.getAssociation(i))

    walk(assoc.getAssociation())
    return found


def _still_used_species(model, skip: set[str]) -> set[str]:
    used = set()
    for i in range(model.getNumReactions()):
        rxn = model.getReaction(i)
        if rxn.getId() in skip:
            continue
        used |= _reaction_species(rxn)
    return used


def _still_used_genes(model, skip: set[str]) -> set[str]:
    used = set()
    for i in range(model.getNumReactions()):
        rxn = model.getReaction(i)
        if rxn.getId() in skip:
            continue
        used |= _gene_products_of(rxn)
    return used


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------


@dataclass
class DeprecationResult:
    """What a deprecation run did, or would do under ``dry_run``."""

    reactions: list[str] = field(default_factory=list)
    metabolites: list[str] = field(default_factory=list)
    genes: list[str] = field(default_factory=list)
    dry_run: bool = False

    def summary(self) -> str:
        verb = "Would remove" if self.dry_run else "Removed"
        lines = [
            f"{verb} {len(self.reactions)} reaction(s), "
            f"{len(self.metabolites)} orphaned metabolite(s), "
            f"{len(self.genes)} orphaned gene(s)."
        ]
        for label, items in (
            ("reactions", self.reactions),
            ("metabolites", self.metabolites),
            ("genes", self.genes),
        ):
            if items:
                lines.append(f"  {label}: {', '.join(sorted(items))}")
        return "\n".join(lines)


def deprecate_reactions(
    ids: Sequence[str],
    reason: str,
    replaced_by: str = "",
    pr: str = "",
    notes: str = "",
    date: str = "",
    model_path: str = MODEL_PATH,
    reactions_tsv: str = REACTIONS_TSV,
    metabolites_tsv: str = METABOLITES_TSV,
    cascade: bool = True,
    dry_run: bool = False,
) -> DeprecationResult:
    """Remove reactions from the model and log them as deprecated.

    Parameters
    ----------
    ids
        Reaction identifiers, with or without the SBML ``R_`` prefix.
    reason
        One of :data:`REASONS`.
    replaced_by
        Identifier(s) that now serve this function, semicolon-separated.
        Required when ``reason`` is ``duplicate`` or ``id_changed``.
    cascade
        Also remove metabolites and genes left with no reactions, logging the
        metabolites with ``reason="orphaned"``. Leave this on unless you have a
        specific reason not to; ``test_sbml.py`` fails on isolated entities.
    dry_run
        Report what would happen and write nothing.
    """
    if reason in REASONS_REQUIRING_REPLACEMENT and not replaced_by:
        raise DeprecationError(
            f"reason={reason!r} requires replaced_by: say what supersedes these "
            f"reactions, otherwise the entry cannot be acted on by a reader."
        )

    doc, model = _load(model_path)
    targets, missing = [], []
    for raw in ids:
        sid = _sbml_reaction_id(raw)
        if model.getReaction(sid) is None:
            missing.append(raw)
        else:
            targets.append(sid)
    if missing:
        raise DeprecationError(
            "not in the model (already removed?): " + ", ".join(missing)
        )
    if not targets:
        return DeprecationResult(dry_run=dry_run)

    target_set = set(targets)
    doomed_species: set[str] = set()
    doomed_genes: set[str] = set()
    if cascade:
        touched_species: set[str] = set()
        touched_genes: set[str] = set()
        for sid in targets:
            rxn = model.getReaction(sid)
            touched_species |= _reaction_species(rxn)
            touched_genes |= _gene_products_of(rxn)
        doomed_species = touched_species - _still_used_species(model, target_set)
        doomed_genes = touched_genes - _still_used_genes(model, target_set)

    rxn_records = [
        DeprecationRecord(
            id=sid,
            name=model.getReaction(sid).getName() or "",
            reason=reason,
            replaced_by=replaced_by,
            pr=pr,
            date=date,
            notes=notes,
        )
        for sid in targets
    ]
    met_records = [
        DeprecationRecord(
            id=sp,
            name=(model.getSpecies(sp).getName() or "") if model.getSpecies(sp) else "",
            reason="orphaned",
            pr=pr,
            date=date,
            notes=f"orphaned by removal of {', '.join(sorted(strip_sbml_prefix(t) for t in targets))}"[
                :140
            ],
        )
        for sp in sorted(doomed_species)
    ]

    result = DeprecationResult(
        reactions=[strip_sbml_prefix(t) for t in targets],
        metabolites=[strip_sbml_prefix(s) for s in sorted(doomed_species)],
        genes=[unescape_sbml_id(strip_sbml_prefix(g)) for g in sorted(doomed_genes)],
        dry_run=dry_run,
    )
    if dry_run:
        return result

    for sid in targets:
        model.removeReaction(sid)
    for sp in doomed_species:
        model.removeSpecies(sp)
    if doomed_genes:
        plugin = model.getPlugin("fbc")
        if plugin is not None:
            for gid in doomed_genes:
                plugin.removeGeneProduct(gid)

    if not libsbml.writeSBMLToFile(doc, model_path):
        raise DeprecationError(f"failed to write {model_path}")

    append_records(reactions_tsv, rxn_records)
    if met_records:
        append_records(metabolites_tsv, met_records)
    return result


def deprecate_metabolites(
    ids: Sequence[str],
    reason: str,
    replaced_by: str = "",
    pr: str = "",
    notes: str = "",
    date: str = "",
    model_path: str = MODEL_PATH,
    metabolites_tsv: str = METABOLITES_TSV,
    dry_run: bool = False,
) -> DeprecationResult:
    """Remove metabolites from the model and log them as deprecated.

    Refuses to remove a metabolite that any reaction still uses -- that would
    leave the model invalid. Remove the reactions first with
    :func:`deprecate_reactions`, which cascades to the metabolites anyway.
    """
    if reason in REASONS_REQUIRING_REPLACEMENT and not replaced_by:
        raise DeprecationError(f"reason={reason!r} requires replaced_by")

    doc, model = _load(model_path)
    targets, missing = [], []
    for raw in ids:
        sid = _sbml_species_id(raw)
        if model.getSpecies(sid) is None:
            missing.append(raw)
        else:
            targets.append(sid)
    if missing:
        raise DeprecationError(
            "not in the model (already removed?): " + ", ".join(missing)
        )

    in_use = _still_used_species(model, skip=set())
    blocked = sorted(set(targets) & in_use)
    if blocked:
        users = {}
        for i in range(model.getNumReactions()):
            rxn = model.getReaction(i)
            for sp in _reaction_species(rxn) & set(blocked):
                users.setdefault(sp, []).append(rxn.getId())
        detail = "; ".join(
            f"{strip_sbml_prefix(sp)} used by "
            f"{', '.join(strip_sbml_prefix(r) for r in sorted(users[sp])[:5])}"
            for sp in blocked
        )
        raise DeprecationError(
            "these metabolites are still used by reactions, remove those first: " + detail
        )

    records = [
        DeprecationRecord(
            id=sid,
            name=model.getSpecies(sid).getName() or "",
            reason=reason,
            replaced_by=replaced_by,
            pr=pr,
            date=date,
            notes=notes,
        )
        for sid in targets
    ]
    result = DeprecationResult(
        metabolites=[strip_sbml_prefix(t) for t in targets], dry_run=dry_run
    )
    if dry_run:
        return result

    for sid in targets:
        model.removeSpecies(sid)
    if not libsbml.writeSBMLToFile(doc, model_path):
        raise DeprecationError(f"failed to write {model_path}")
    append_records(metabolites_tsv, records)
    return result


# --------------------------------------------------------------------------
# Stamping the PR number after the fact
# --------------------------------------------------------------------------


def stamp_pr_number(
    number: str,
    reactions_tsv: str = REACTIONS_TSV,
    metabolites_tsv: str = METABOLITES_TSV,
) -> dict[str, list[str]]:
    """Fill in the ``pr`` column for rows that do not have one yet.

    You cannot know your pull request number before you open the pull request,
    so ``--pr`` is optional when you deprecate something. CI closes the loop:
    the ``Custom-CI`` workflow runs this on every pull request and commits the
    result, the same way it already stamps the PR number into
    ``scripts/results/README.md``.

    Only blank cells are filled, so re-running is safe and a row that was
    deliberately attributed to a different PR or issue is never overwritten.

    Returns a mapping of file path to the ids that were stamped.
    """
    if not number:
        raise DeprecationError("no PR number given")
    pr = "#" + str(number).lstrip("#")
    if not re.fullmatch(r"#\d+", pr):
        raise DeprecationError(f"{number!r} does not look like a PR number")

    stamped: dict[str, list[str]] = {}
    for path in (reactions_tsv, metabolites_tsv):
        records = read_records(path)
        touched = [r for r in records if not r.pr]
        if not touched:
            continue
        for record in touched:
            record.pr = pr
        write_records(path, records)
        stamped[path] = [r.id for r in touched]
    return stamped


# --------------------------------------------------------------------------
# Notes mirror (shared with scripts/export_model.py)
# --------------------------------------------------------------------------

NOTES_REACTION_KEY = "DEPRECATED_REACTIONS"
NOTES_METABOLITE_KEY = "DEPRECATED_METABOLITES"
NOTES_INFO_KEY = "DEPRECATED_INFO"
NOTES_INFO_VALUE = (
    "reactions and metabolites listed above were deliberately removed during "
    "curation; full table with reasons at "
    "https://github.com/C-CoMP-STC/GEM-mit1002/tree/main/data/deprecatedIdentifiers"
)


def build_notes_dict(
    reactions_tsv: str = REACTIONS_TSV, metabolites_tsv: str = METABOLITES_TSV
) -> dict:
    """The ``model.notes`` mapping that mirrors the deprecated-identifier TSVs.

    Only the identifier lists and a pointer -- ``<notes>`` is a flat
    string-to-string map, so the structured detail stays in the TSVs.
    """
    rxns = sorted(r.id for r in read_records(reactions_tsv))
    mets = sorted(r.id for r in read_records(metabolites_tsv))
    notes = {}
    if rxns:
        notes[NOTES_REACTION_KEY] = "; ".join(rxns)
    if mets:
        notes[NOTES_METABOLITE_KEY] = "; ".join(mets)
    if rxns or mets:
        notes[NOTES_INFO_KEY] = NOTES_INFO_VALUE
    return notes


def sync_model_notes(
    model_path: str = MODEL_PATH,
    reactions_tsv: str = REACTIONS_TSV,
    metabolites_tsv: str = METABOLITES_TSV,
) -> dict:
    """Write the deprecated-identifier mirror into the model's ``<notes>``.

    Preserves any other ``<p>KEY: VALUE</p>`` notes already present, and only
    touches the file if the result differs, so it is safe in CI.
    """
    doc, model = _load(model_path)
    wanted = build_notes_dict(reactions_tsv, metabolites_tsv)

    existing = {}
    if model.isSetNotes():
        text = model.getNotesString()
        for key, value in re.findall(r"<p[^>]*>\s*([A-Za-z0-9_]+)\s*:\s*(.*?)\s*</p>", text, re.S):
            existing[key] = " ".join(value.split())

    merged = dict(existing)
    merged.update(wanted)
    for key in (NOTES_REACTION_KEY, NOTES_METABOLITE_KEY, NOTES_INFO_KEY):
        if key not in wanted:
            merged.pop(key, None)

    if merged == existing:
        return merged

    body = "".join(f"<p>{k}: {_escape(v)}</p>" for k, v in merged.items())
    model.setNotes(f'<body xmlns="http://www.w3.org/1999/xhtml">{body}</body>')
    if not libsbml.writeSBMLToFile(doc, model_path):
        raise DeprecationError(f"failed to write {model_path}")
    return merged


def _escape(text: str) -> str:
    return (
        text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tools.deprecate",
        description=(
            "Remove reactions or metabolites from model.xml and record them in "
            "data/deprecatedIdentifiers/ in one step."
        ),
    )
    sub = parser.add_subparsers(dest="kind", required=True)

    for kind, helptext in (
        ("reaction", "remove reactions (cascades to orphaned mets and genes)"),
        ("metabolite", "remove metabolites that no reaction uses"),
    ):
        p = sub.add_parser(kind, help=helptext)
        p.add_argument("ids", nargs="+", help="identifiers, R_/M_ prefix optional")
        p.add_argument("--reason", required=True, choices=REASONS)
        p.add_argument(
            "--replaced-by",
            default="",
            help="identifier(s) that supersede these, semicolon-separated",
        )
        p.add_argument(
            "--pr",
            default="",
            help=(
                "pull request number, e.g. '#317'. Usually omit this: CI stamps "
                "it in once the PR exists. Pass it only if you already know the "
                "relevant number, e.g. an issue this closes."
            ),
        )
        p.add_argument("--notes", default="", help="one short line; prose goes in the PR")
        p.add_argument("--date", default="", help="YYYY-MM-DD, defaults to today")
        p.add_argument("--model", default=MODEL_PATH)
        p.add_argument(
            "--dry-run",
            action="store_true",
            help="report what would happen and write nothing",
        )
        if kind == "reaction":
            p.add_argument(
                "--no-cascade",
                action="store_true",
                help="do not remove metabolites/genes orphaned by this removal",
            )

    sync = sub.add_parser(
        "sync-notes", help="regenerate the model's <notes> mirror from the TSVs"
    )
    sync.add_argument("--model", default=MODEL_PATH)

    stamp = sub.add_parser(
        "stamp-pr",
        help="fill in blank pr cells with a PR number (CI runs this for you)",
    )
    stamp.add_argument("number", help="pull request number, with or without '#'")

    args = parser.parse_args(argv)

    try:
        if args.kind == "sync-notes":
            notes = sync_model_notes(model_path=args.model)
            print(f"Model notes now carry {len(notes)} key(s): {', '.join(notes)}")
            return 0
        if args.kind == "stamp-pr":
            stamped = stamp_pr_number(args.number)
            if not stamped:
                print("No blank pr cells to fill.")
            for path, ids in stamped.items():
                print(f"Stamped {len(ids)} row(s) in {os.path.basename(path)}: {', '.join(ids)}")
            return 0
        if args.kind == "reaction":
            result = deprecate_reactions(
                args.ids,
                reason=args.reason,
                replaced_by=args.replaced_by,
                pr=args.pr,
                notes=args.notes,
                date=args.date,
                model_path=args.model,
                cascade=not args.no_cascade,
                dry_run=args.dry_run,
            )
        else:
            result = deprecate_metabolites(
                args.ids,
                reason=args.reason,
                replaced_by=args.replaced_by,
                pr=args.pr,
                notes=args.notes,
                date=args.date,
                model_path=args.model,
                dry_run=args.dry_run,
            )
    except DeprecationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(result.summary())
    if not result.dry_run:
        print("Remember to run the test suite before opening your PR: pytest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
