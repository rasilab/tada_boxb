"""Shared biological annotations, read parsing, and count pooling for the λN DMS reanalysis.

Biological identities and coordinates come from annotation CSVs. Shared processing choices
are declared below. Run all scripts from the repository root.
"""

import gzip
from contextlib import ExitStack
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, t as student_t

ANALYSIS_DIRECTORY = "analysis/ln_dms_sequencing"
ANNOTATION_DIRECTORY = f"{ANALYSIS_DIRECTORY}/annotations"
LINKAGE_DIRECTORY = f"{ANALYSIS_DIRECTORY}/data/barcode_linkage"
COUNT_DIRECTORY = f"{ANALYSIS_DIRECTORY}/data/recorder_counts"
# Shared read-processing choices and figure inclusion thresholds.
MIN_BASE_QUALITY = 30
MIN_LINKAGE_SUPPORT = 3
MIN_LINKAGE_PURITY = .90
MIN_REPORTED_SUPPORT = 2
MAX_OFF_DESIGN_SHARE = 1.00
EXCLUDE_NEIGHBOURS = True
FIGURE_MIN_LINKAGE_SUPPORT = 5
MIN_CALLABLE_READS = 1
# Barcodes required at both sites before a codon variant enters the class distribution and its
# source-data table. Shared so the panel and its exported values cannot drift apart.
CLASS_MIN_BARCODES = 5
POOLED_WT_LABEL = "WT"
EXPECTED_POSITIONS = 22
EXPECTED_INSERTS = 1387
EXPECTED_DISPLAY_CELLS = 461
# Barcodes of one variant disagree far more than their pooled read counts imply, so every test
# against wild type takes the barcode as its unit and reads that dispersion off the wild-type
# barcodes themselves. The read-level test is kept beside it for the record, not for the figures.
SIGNIFICANCE_LEVEL = .05
PHRED_OFFSET = 33
REFERENCE_BASE = "T"
EDITED_BASE = "C"
MASK = "N"
BASES = ["A", "C", "G", "T", "N"]
STOP_SYMBOL = "*"
WILD_TYPE_CODON_COLUMN = "wt_codon"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def annotation(name):
    return pd.read_csv(f"{ANNOTATION_DIRECTORY}/{name}.csv", comment="#", dtype={"expected_sequence": str})


def manifest(role=None, dataset=None):
    """The provenance manifest, optionally restricted to one role or dataset."""
    rows = annotation("dataset_manifest")
    require(rows.sample_id.is_unique, "Duplicate sample in the dataset manifest")
    if role is not None:
        rows = rows.loc[rows.role.eq(role)]
    if dataset is not None:
        rows = rows.loc[rows.dataset.eq(dataset)]
    require(len(rows) > 0, f"No manifest rows for role={role} dataset={dataset}")
    return rows.reset_index(drop=True)


def figure_sample():
    """The one RNA sample selected in the biological provenance annotation."""
    samples = manifest(role="rna")
    require(samples.use_for_figures.isin(["yes", "no"]).all(), "Invalid figure sample selection")
    selected = samples.loc[samples.use_for_figures.eq("yes"), "sample_id"]
    require(len(selected) == 1, "Select exactly one RNA sample for the retained figures")
    return selected.item()


def comparison_sites():
    """The two annotated recorder sites, ordered along the RNA."""
    sites = measured_sites(primary_only=True)
    require(sites.compare_sites.isin(["yes", "no"]).all(), "Invalid comparison-site selection")
    selected = sites.loc[sites.compare_sites.eq("yes")].sort_values("rna_position")
    require(len(selected) == 2, "Select exactly two recorder sites for comparison")
    return list(selected.site_id)


def read_files(sample_id):
    """Declared FASTQ paths for one sample, keyed by read name."""
    rows = annotation("read_files")
    rows = rows.loc[rows.sample_id.eq(sample_id)]
    require(len(rows) > 0, f"No declared read files for {sample_id}")
    require(rows.read.is_unique, f"Duplicate read declaration for {sample_id}")
    return dict(zip(rows.read, rows.path))


def layout(sample_id):
    """Read layout of one sample as {feature: row}, with 0-based `offset` added."""
    rows = annotation("read_layouts")
    rows = rows.loc[rows.sample_id.eq(sample_id)].copy()
    require(len(rows) > 0, f"No declared read layout for {sample_id}")
    require(rows.feature.is_unique, f"Duplicate feature in the layout of {sample_id}")
    rows["offset"] = rows.start - 1
    for row in rows.itertuples():
        if isinstance(row.expected_sequence, str) and row.expected_sequence != "NA":
            require(len(row.expected_sequence) == row.length,
                    f"Declared length disagrees with the expected sequence of {sample_id} {row.feature}")
    return {row.feature: row for row in rows.itertuples()}


def measured_sites(feature=None, primary_only=False):
    """Annotated adenosines, in the read coordinates the counter measures them at."""
    sites = annotation("measured_adenosines")
    require(sites.site_id.is_unique, "Duplicate measured adenosine")
    if feature is not None:
        sites = sites.loc[sites.feature.eq(feature)]
    if primary_only:
        sites = sites.loc[sites.in_primary_endpoint.eq("yes")]
    return sites.reset_index(drop=True)


def verify_sites(sample_id):
    """Check every annotated adenosine against the expected sequence of its feature.

    An editable site must carry the reference base in the construct the layout declares, its read
    position must follow from the feature's start and the site's index, and its RNA position must
    follow from the feature's declared `rna_offset`. This is what keeps a coordinate change from
    silently measuring a neighbouring base. The offset is declared per feature because the read
    regions do not all begin at RNA position one: the 22-nt boxB region carries the 19-nt hairpin
    plus flanking bases, so its offset is 21 while the recorder's is its own length.
    """
    features = layout(sample_id)
    for site in measured_sites().itertuples():
        if site.feature not in features:
            continue
        row = features[site.feature]
        require(row.expected_sequence[site.feature_index] == REFERENCE_BASE,
                f"{site.site_id} is not {REFERENCE_BASE} in the {site.feature} of {sample_id}")
        require(row.start + site.feature_index == site.read_position,
                f"{site.site_id} read position disagrees with the {site.feature} start of {sample_id}")
        require(int(row.rna_offset) - site.feature_index == site.rna_position,
                f"{site.site_id} RNA position disagrees with the declared offset of the {site.feature} of {sample_id}")
    return True


def inserts():
    """Designed λN inserts with their codon and protein consequences.

    The identical wild-type DNA design is retained once. Mapping is by complete DNA sequence;
    the positional labels are display coordinates. Validate the supplied codon annotations.
    """
    designs = pd.read_csv(f"{ANNOTATION_DIRECTORY}/insert_annotations.csv")
    require(len(designs) == EXPECTED_INSERTS, "Unexpected number of designed inserts")
    require(designs.row_num.is_unique and designs.insert_seq.is_unique, "Designed inserts are not unique")
    require(designs.codon_position.between(1, EXPECTED_POSITIONS).all(), "Codon position outside the λN segment")
    require(len(designs[["codon_position", "mut_aa"]].drop_duplicates()) == EXPECTED_DISPLAY_CELLS,
            "Unexpected number of position and amino-acid display cells")
    residues = designs[["codon_position", "wt_aa", "wt_codon"]].drop_duplicates().sort_values("codon_position")
    require(residues.codon_position.is_unique and len(residues) == EXPECTED_POSITIONS,
            "Missing or conflicting wild-type residue annotations")
    require(str(Seq("".join(residues.wt_codon)).translate()) == "".join(residues.wt_aa),
            "Wild-type codon and amino-acid annotations disagree")
    designs["variant_class"] = np.where(
        designs.mut_aa.eq(STOP_SYMBOL), "stop",
        np.where(designs.wt_aa.eq(designs.mut_aa), "synonymous", "missense"))
    wild_type = designs.loc[designs.mut_codon.eq(designs[WILD_TYPE_CODON_COLUMN])]
    require(len(wild_type) == 1, "The design does not carry exactly one canonical wild-type row")
    require(wild_type.insert_seq.iloc[0] == "".join(residues.wt_codon),
            "Canonical wild-type insert disagrees with its codon annotations")
    designs["is_wild_type"] = designs.insert_seq.eq(wild_type.insert_seq.iloc[0])
    return designs


def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def synchronized_fastq_records(paths):
    """Stream aligned read files, rejecting malformed records or unequal/mismatched pairs."""
    with ExitStack() as stack:
        handles = {name: stack.enter_context(open_fastq(path)) for name, path in paths.items()}
        while True:
            titles = {name: handle.readline() for name, handle in handles.items()}
            if not any(titles.values()):
                return
            require(all(titles.values()), "FASTQ read files end at different records")
            require(len({title.split()[0] for title in titles.values()}) == 1,
                    "FASTQ read identifiers are not synchronized")
            record = {}
            for name, handle in handles.items():
                sequence = handle.readline().rstrip("\n")
                separator = handle.readline()
                qualities = handle.readline().rstrip("\n")
                require(titles[name].startswith("@") and separator.startswith("+")
                        and len(sequence) > 0 and len(sequence) == len(qualities),
                        f"Malformed or truncated FASTQ record in {paths[name]}")
                record[name] = sequence, qualities
            yield record


def count_columns(site_ids):
    """Per-base counts retain each site's denominator: reads minus that site's N count."""
    return ["reads"] + [f"pos_{site}_{base}" for site in site_ids for base in BASES]


def site_estimates(counts, site_ids):
    """Edited and callable counts per site, taking each site's own denominator."""
    frames = []
    for site in site_ids:
        callable_count = counts.reads - counts[f"pos_{site}_{MASK}"]
        frames.append(pd.DataFrame({
            "site_id": site,
            "edited": counts[f"pos_{site}_{EDITED_BASE}"],
            "callable": callable_count,
        }, index=counts.index))
    return pd.concat(frames)


def filter_linkage_support(counts, sample_id, minimum):
    """Apply a stricter linkage-read floor to accepted RNA barcodes without recounting FASTQs.

    The source counts retain their production linkage policy. This selection can only remove
    barcodes; their assignments and all other production filters remain unchanged.
    """
    require(np.isfinite(minimum) and minimum == int(minimum) and minimum >= MIN_LINKAGE_SUPPORT,
            "The plot linkage floor must be an integer at least as strict as the source counts")
    row = manifest(role="rna").set_index("sample_id").loc[sample_id]
    linkage = pd.read_csv(f"{LINKAGE_DIRECTORY}/{row.linkage_dataset}.csv").set_index("barcode", verify_integrity=True)
    require(counts.barcode.is_unique and counts.barcode.isin(linkage.index).all(),
            "An RNA barcode is repeated or missing from the linkage table")
    evidence = linkage.loc[counts.barcode]
    require(evidence.status.eq("assigned").all()
            and np.array_equal(evidence.insert_num.to_numpy(), counts.insert_num.to_numpy()),
            "RNA barcode assignments disagree with the production linkage")
    keep = evidence.designed_support.ge(minimum).to_numpy()
    filtered = counts.loc[keep].copy()
    require(len(filtered) > 0, "No RNA barcode reaches the selected linkage support")
    print(f"{sample_id}: linkage support >= {minimum} designed reads, retaining "
          f"{len(filtered):,}/{len(counts):,} barcodes and {filtered.reads.sum():,}/{counts.reads.sum():,} RNA reads; "
          f"production purity >= {MIN_LINKAGE_PURITY:g} and other filters unchanged.")
    return filtered


def variant_barcode_sites(sample_id, site_ids, synonymous_group, min_linkage_support,
                          variant_resolution="amino_acid"):
    """Per-barcode edited and callable counts at each site, labelled by designed variant.

    This is the frame `pooled_variant_sites` sums. A test against wild type reads it unpooled,
    because the barcode, not the read, is the unit those counts actually vary over.
    """
    verify_sites(sample_id)
    require(len(site_ids) == len(set(site_ids)) and set(site_ids) <= set(measured_sites().site_id),
            "Unknown or repeated editing site")
    designs = inserts()
    require(variant_resolution in {"amino_acid", "codon"}, "Unknown variant pooling resolution")
    if variant_resolution == "codon":
        designs["variant"] = designs.wt_codon + designs.codon_position.astype(str) + designs.mut_codon
        require(synonymous_group not in set(designs.variant), "Wild-type label duplicates a codon mutation")
        designs.loc[designs.is_wild_type, "variant"] = synonymous_group
        require(designs.variant.is_unique, "Distinct designed codons were merged")
    else:
        designs["variant"] = designs.wt_aa + designs.codon_position.astype(str) + designs.mut_aa
        require(synonymous_group not in set(designs.variant), "Synonymous group duplicates a mutation")
        synonymous = designs.variant_class.eq("synonymous")
        designs.loc[synonymous, "variant"] = synonymous_group
    groups = designs.groupby(["variant", "variant_class"]).agg(designed_codons=("row_num", "size"))
    require(groups.reset_index().variant.is_unique, "A variant has conflicting classes")
    columns = ["reads"] + [f"pos_{site}_{base}" for site in site_ids for base in BASES]
    counts = pd.read_csv(f"{COUNT_DIRECTORY}/{sample_id}_reads.csv",
                         usecols=["barcode", "insert_num", *columns])
    require(len(counts) > 0 and counts.notna().all().all() and counts.barcode.is_unique,
            "Empty, missing, or duplicated barcode counts")
    require(np.isfinite(counts[columns]).all().all() and counts[columns].ge(0).all().all()
            and counts[columns].eq(np.floor(counts[columns])).all().all(), "Invalid read counts")
    counts = filter_linkage_support(counts, sample_id, min_linkage_support)
    counts = counts.merge(designs[["row_num", "variant"]], left_on="insert_num", right_on="row_num",
                          how="left", validate="many_to_one")
    require(counts.variant.notna().all(), "A counted insert has no designed variant")
    for site in site_ids:
        require(counts[[f"pos_{site}_{base}" for base in BASES]].sum(axis=1).eq(counts.reads).all(),
                f"Base counts at {site} do not sum to accepted reads")
    sites = site_estimates(counts, site_ids).join(counts[["variant", "reads", "barcode", "insert_num"]])
    require(sites.edited.between(0, sites["callable"]).all(), "Edits exceed site-callable reads")
    return sites, groups, counts


def pooled_variant_sites(sample_id, site_ids, synonymous_group, min_linkage_support,
                         variant_resolution="amino_acid"):
    """Pool barcode counts at the selected variant resolution before taking site fractions.

    Codon resolution retains every distinct designed DNA insert in all classes, including one
    canonical WT. Amino-acid resolution pools missense/stop codons at each position; synonymous
    designs form one protein group.
    The selected linkage floor applies before pooling; no RNA coverage or effect filter applies.
    Uncovered designed variants retain zero counts and undefined fractions, not zero editing.
    """
    sites, groups, counts = variant_barcode_sites(sample_id, site_ids, synonymous_group,
                                                  min_linkage_support, variant_resolution)
    summary = sites.groupby(["variant", "site_id"]).agg(
        reads=("reads", "sum"), edited=("edited", "sum"), callable=("callable", "sum"),
        barcodes=("barcode", "nunique"), observed_codons=("insert_num", "nunique"))
    grid = groups.reset_index().merge(pd.DataFrame({"site_id": site_ids}), how="cross")
    pooled = grid.merge(summary, on=["variant", "site_id"], how="left", validate="one_to_one")
    totals = ["reads", "edited", "callable", "barcodes", "observed_codons"]
    pooled[totals] = pooled[totals].fillna(0).astype("int64")
    pooled["fraction"] = pooled.edited / pooled["callable"].replace(0, np.nan)
    for site in site_ids:
        selected = pooled.loc[pooled.site_id.eq(site)]
        require(selected.reads.sum() == counts.reads.sum()
                and selected.edited.sum() == counts[f"pos_{site}_{EDITED_BASE}"].sum()
                and selected["callable"].sum() == (counts.reads - counts[f"pos_{site}_{MASK}"]).sum()
                and selected.barcodes.sum() == len(counts)
                and selected.observed_codons.sum() == counts.insert_num.nunique(),
                f"Variant pooling changed counts at {site}")
    return pooled


def hamming_one_groups(barcodes):
    """Barcodes within one substitution of another barcode in the same list.

    Two equal-length sequences at Hamming distance one agree after deleting the differing
    position, so one pass per position groups them without enumerating neighbours.
    """
    barcodes = list(barcodes)
    require(len(set(len(barcode) for barcode in barcodes)) == 1, "Barcodes of mixed length")
    neighbours = set()
    for position in range(len(barcodes[0])):
        seen = {}
        for barcode in barcodes:
            signature = barcode[:position] + barcode[position + 1:]
            if signature in seen:
                neighbours.add(barcode)
                neighbours.add(seen[signature])
            else:
                seen[signature] = barcode
    return neighbours


BASE_CODES = np.full(256, len(BASES) - 1, dtype=np.uint8)
for _index, _base in enumerate(BASES[:-1]):
    BASE_CODES[ord(_base)] = _index


def encode(text):
    return np.frombuffer(text.encode(), dtype=np.uint8)


class FeatureCheck:
    """One declared constant region, prepared for numpy comparison against a read."""

    def __init__(self, name, read, offset, expected):
        self.name = name
        self.read = read
        self.offset = offset
        self.expected = encode(expected)
        self.length = len(expected)
        # Every reference T is an adenosine on the RNA and so may legitimately read as C in an
        # edited molecule, whether or not the site annotation selects it for measurement. Deriving
        # the allowance from the construct rather than from the site table is what keeps an
        # editing-aware check from rejecting edited reads at an unannotated adenosine.
        self.edit_allowed = self.expected == ord(REFERENCE_BASE)

    def evaluate(self, sequence, quality):
        """Return (strict, editing_aware, masked_bases) or None when the read is too short."""
        if len(sequence) < self.offset + self.length:
            return None
        region = sequence[self.offset:self.offset + self.length]
        low = quality[self.offset:self.offset + self.length] < MIN_BASE_QUALITY
        mismatch = (region != self.expected) & ~low
        allowed = mismatch & self.edit_allowed & (region == ord(EDITED_BASE))
        return not mismatch.any(), not (mismatch & ~allowed).any(), int(low.sum())


def feature_checks(features, roles=("required",)):
    return [FeatureCheck(name, row.read, row.offset, row.expected_sequence)
            for name, row in features.items() if row.role in roles and isinstance(row.expected_sequence, str)]


SEQUENCE_CODES = str.maketrans("ACGT", "0123")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


DNA_CODES = np.zeros(256, dtype=bool)
for _base in "ACGT":
    DNA_CODES[ord(_base)] = True


def encode_sequence(sequence):
    """Pack a fixed-length ACGT sequence into an exact integer key.

    Base-4 packing is exact, so two different sequences never share a key, unlike a hash. The
    length is fixed by the layout, so leading zeros carry no ambiguity.
    """
    return int(sequence.translate(SEQUENCE_CODES), 4)


def packable(region):
    """Whether an encoded region holds only A, C, G, and T and can therefore be packed.

    A no-call can reach a measurement with a passing quality score, and packing it would raise, so
    the caller counts the read against its own unreadable-field tally instead.
    """
    return bool(DNA_CODES[region].all())


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def designed_insert_keys(orientation):
    """Map the packed key of every designed insert to its row number, in the read's orientation."""
    designs = inserts()
    require(designs.insert_seq.str.translate(COMPLEMENT).str[::-1].equals(designs.insert_seq_rc),
            "The stored reverse complement is not the reverse complement of the designed insert")
    column = "insert_seq_rc" if orientation == "antisense" else "insert_seq"
    require(designs[column].is_unique, f"Designed {column} values are not unique")
    keys = {encode_sequence(sequence): row for sequence, row in zip(designs[column], designs.row_num)}
    require(len(keys) == len(designs), "Designed inserts collide after packing")
    return keys


class ConstantRegions:
    """Every required constant position of one read, compared in a single vectorised pass.

    A read is checked feature by feature only when the combined comparison fails, which is what
    keeps the attribution of a rejection reason off the path that most reads take.
    """

    def __init__(self, checks, read):
        selected = [check for check in checks if check.read == read]
        require(len(selected) > 0, f"No constant region declared on read {read}")
        self.read = read
        self.checks = selected
        self.offsets = np.concatenate([np.arange(check.offset, check.offset + check.length)
                                       for check in selected])
        self.expected = np.concatenate([check.expected for check in selected])
        self.edit_allowed = np.concatenate([check.edit_allowed for check in selected])
        self.extent = int(self.offsets.max()) + 1

    def evaluate(self, bases, quality):
        """Return (strict, editing_aware), or None when the read is too short to carry the regions."""
        if len(bases) < self.extent:
            return None
        region = bases[self.offsets]
        low = quality[self.offsets] < MIN_BASE_QUALITY
        mismatch = (region != self.expected) & ~low
        if not mismatch.any():
            return True, True
        allowed = self.edit_allowed & (region == ord(EDITED_BASE))
        return False, not (mismatch & ~allowed).any()


def constant_regions(checks):
    return [ConstantRegions(checks, read) for read in sorted({check.read for check in checks})]


def rejection_reason(checks, encoded):
    """The first declared constant region an otherwise accepted read fails, for the accounting."""
    for check in checks:
        result = check.evaluate(*encoded[check.read])
        if result is None:
            return "short_read"
        if not result[1]:
            return f"structure_{check.name}"
    return ""


def barcode_dispersion(edited, callable_counts, reference_rate):
    """Intraclass correlation of editing among barcodes measured at one rate.

    A binomial model says two barcodes of one variant agree to within their read counts. They do
    not: the wild-type barcodes disagree by tens of percentage points however deeply each is
    sequenced, and this correlation is the size of that disagreement. It is read off the wild-type
    barcodes and then applied to every variant, which is what lets a variant carried by two
    barcodes be tested at all.
    """
    edited, callable_counts = np.asarray(edited, float), np.asarray(callable_counts, float)
    covered = callable_counts > 0
    edited, callable_counts = edited[covered], callable_counts[covered]
    require(len(edited) > 1 and 0 < reference_rate < 1, "Dispersion needs covered barcodes and a rate")
    expected = callable_counts * reference_rate
    pearson = (((edited - expected) ** 2) / (expected * (1 - reference_rate))).sum() / (len(edited) - 1)
    mean_size = (callable_counts ** 2).sum() / callable_counts.sum()
    require(mean_size > 1, "Barcodes carry too few reads to estimate dispersion")
    return float(max((pearson - 1) / (mean_size - 1), 0))


def clustered_rate_error(edited, callable_counts, reference_rate, dispersion):
    """Standard error of a pooled rate once its barcodes are counted as the unit, not its reads.

    The design effect is `1 + dispersion x (mean barcode size - 1)`: deeper sequencing of the same
    barcodes buys almost nothing, so the error is set by how many barcodes carry the variant.
    """
    callable_counts = np.asarray(callable_counts, float)
    callable_counts = callable_counts[callable_counts > 0]
    total = callable_counts.sum()
    if total == 0:
        return np.nan
    mean_size = (callable_counts ** 2).sum() / total
    design_effect = 1 + dispersion * (mean_size - 1)
    return float(np.sqrt(design_effect * reference_rate * (1 - reference_rate) / total))


def clustered_rate_test(edited, callable_counts, reference_rate, reference_error, dispersion):
    """Two-sided test of one variant's pooled rate against the wild-type rate.

    Both rates are the depth-weighted pooled rates the figures draw; only their errors change.
    """
    edited, callable_counts = np.asarray(edited, float), np.asarray(callable_counts, float)
    total = callable_counts.sum()
    if total == 0:
        return np.nan, np.nan, np.nan
    rate = edited.sum() / total
    error = np.sqrt(clustered_rate_error(edited, callable_counts, reference_rate, dispersion) ** 2
                    + reference_error ** 2)
    statistic = (rate - reference_rate) / error
    # The error is carried by however many barcodes hold the variant, so the tail is Student's with
    # that many degrees of freedom; a normal tail would read a two-barcode variant as if it were many.
    degrees = max(len(callable_counts[callable_counts > 0]) - 1, 1)
    return float(rate), float(statistic), float(2 * student_t.sf(abs(statistic), degrees))


def read_rate_test(edited, callable_count, reference_edited, reference_callable):
    """Fisher exact test on pooled reads, which treats PCR duplicates as independent."""
    table = [[edited, callable_count - edited],
             [reference_edited, reference_callable - reference_edited]]
    return float(fisher_exact(table).pvalue)


def adjust_bh(p_values):
    """Benjamini-Hochberg false discovery rates, leaving undefined tests undefined."""
    values = pd.Series(p_values, dtype=float)
    defined = values.dropna().sort_values()
    ranks = np.arange(1, len(defined) + 1)
    adjusted = np.minimum.accumulate((defined.to_numpy() * len(defined) / ranks)[::-1])[::-1]
    return values.where(values.isna(), pd.Series(np.minimum(adjusted, 1), index=defined.index))


def position_mean_test(rates, reference_rate, reference_error):
    """Two-sided test of the equal-weight mean of a position's variant rates against wild type.

    The figure averages the substitutions at a position with equal weight, so the substitutions
    are the units here; their spread already carries the barcode noise inside each of them.
    """
    rates = np.asarray(pd.Series(rates).dropna(), float)
    require(len(rates) > 1, "A position mean needs at least two covered variants")
    error = np.sqrt(rates.var(ddof=1) / len(rates) + reference_error ** 2)
    statistic = (rates.mean() - reference_rate) / error
    return float(rates.mean()), float(statistic), float(2 * student_t.sf(abs(statistic), len(rates) - 1))


def wild_type_reference(barcode_sites, synonymous_group):
    """The wild-type rate, its barcode-level dispersion, and its clustered error, at each site."""
    references = []
    for site, site_frame in barcode_sites.groupby("site_id", sort=False):
        wild_type = site_frame.loc[site_frame.variant.eq(synonymous_group) & site_frame["callable"].gt(0)]
        require(len(wild_type) > 1 and wild_type["callable"].sum() > 0,
                "The wild-type group carries too few barcodes to calibrate a test")
        rate = wild_type.edited.sum() / wild_type["callable"].sum()
        dispersion = barcode_dispersion(wild_type.edited, wild_type["callable"], rate)
        references.append(dict(site_id=site, reference_rate=rate, dispersion=dispersion,
                               reference_error=clustered_rate_error(wild_type.edited, wild_type["callable"],
                                                                    rate, dispersion),
                               reference_barcodes=len(wild_type),
                               reference_edited=int(wild_type.edited.sum()),
                               reference_callable=int(wild_type["callable"].sum())))
    return pd.DataFrame(references)


def wild_type_comparisons(sample_id, site_ids, synonymous_group, min_linkage_support,
                          variant_resolution="amino_acid", combine_sites=False, min_barcodes=1):
    """Test every designed variant against the pooled synonymous wild type.

    `rate` and `reference_rate` are the pooled, depth-weighted percentages the figures draw; only
    their errors differ between the two tests reported. The barcode test counts the barcodes
    carrying a variant as its sample size, through the wild-type dispersion; the read test counts
    its reads, which is what makes nearly every variant significant. Figures mark the first.
    `min_barcodes` drops variants carried by fewer barcodes before testing, and is 1 when a figure
    applies no such floor.
    """
    sites, _, _ = variant_barcode_sites(sample_id, site_ids, synonymous_group, min_linkage_support,
                                        variant_resolution)
    frame = sites.reset_index(drop=True)[["variant", "barcode", "site_id", "edited", "callable"]]
    if combine_sites:
        frame = frame.groupby(["variant", "barcode"], as_index=False)[["edited", "callable"]].sum()
        frame["site_id"] = "+".join(site_ids)
    references = wild_type_reference(frame, synonymous_group)
    results = []
    for site, site_frame in frame.groupby("site_id", sort=False):
        wild_type = site_frame.loc[site_frame.variant.eq(synonymous_group)]
        row = references.set_index("site_id").loc[site]
        reference_rate, dispersion, reference_error = row.reference_rate, row.dispersion, row.reference_error
        for variant, group in site_frame.groupby("variant", sort=False):
            covered = group.loc[group["callable"].gt(0)]
            if len(covered) < min_barcodes or covered["callable"].sum() == 0:
                continue
            rate, statistic, p_barcode = clustered_rate_test(
                covered.edited, covered["callable"], reference_rate, reference_error, dispersion)
            results.append(dict(
                variant=variant, site_id=site, barcodes=len(covered),
                edited=int(covered.edited.sum()), callable=int(covered["callable"].sum()),
                rate=rate, statistic=statistic, p_value_barcode=p_barcode,
                p_value_read=read_rate_test(int(covered.edited.sum()), int(covered["callable"].sum()),
                                            int(wild_type.edited.sum()), int(wild_type["callable"].sum()))))
    results = pd.DataFrame(results)
    require(len(results) > 0, "No variant reached the barcode floor")
    for source, target in [("p_value_barcode", "q_value_barcode"), ("p_value_read", "q_value_read")]:
        tested = results.variant.ne(synonymous_group)
        results[target] = np.nan
        results.loc[tested, target] = adjust_bh(results.loc[tested, source]).to_numpy()
    results["significant"] = results.q_value_barcode.lt(SIGNIFICANCE_LEVEL)
    results["direction"] = np.sign(results.rate - results.site_id.map(
        references.set_index("site_id").reference_rate))
    return results, references
