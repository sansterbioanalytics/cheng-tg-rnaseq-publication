#!/usr/bin/env python3
"""Curate existing notebook figures into a publication-facing folder."""

from __future__ import annotations

import csv
import shutil
import struct
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "publication_figures"


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError(f"{path} is not a valid PNG")
    return struct.unpack(">II", header[16:24])


def copy_panel(entry: dict[str, str]) -> dict[str, str]:
    source = ROOT / entry["source"]
    destination = OUT_DIR / entry["tier_dir"] / entry["group_dir"] / entry["curated_name"]
    destination.parent.mkdir(parents=True, exist_ok=True)

    if source.exists():
        shutil.copy2(source, destination)
        entry["status"] = "copied"
        entry["curated_path"] = str(destination.relative_to(ROOT))
    else:
        entry["status"] = "missing"
        entry["curated_path"] = ""

    return entry


def layout_scale(
    src_width: int,
    src_height: int,
    box_x: int,
    box_y: int,
    box_width: int,
    box_height: int,
) -> tuple[int, int, int, int]:
    scale = min(box_width / src_width, box_height / src_height)
    width = int(src_width * scale)
    height = int(src_height * scale)
    x = box_x + (box_width - width) // 2
    y = box_y + (box_height - height) // 2
    return x, y, width, height


def svg_header(width: int = 2550, height: int = 3300) -> list[str]:
    return [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="8.5in" height="11in" viewBox="0 0 {width} {height}">'
        ),
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="white"/>',
    ]


def svg_footer() -> list[str]:
    return ["</svg>"]


def add_text(lines: list[str], x: int, y: int, text: str, size: int, weight: str = "normal") -> None:
    safe = (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
    lines.append(
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="#111111">{safe}</text>'
    )


def add_panel(
    lines: list[str],
    svg_path: Path,
    relative_image: str,
    label: str,
    title: str,
    box_x: int,
    box_y: int,
    box_width: int,
    box_height: int,
) -> None:
    image_path = svg_path.parent / relative_image
    src_width, src_height = png_size(image_path)
    add_text(lines, box_x, box_y - 28, label, 34, "bold")
    add_text(lines, box_x + 55, box_y - 28, title, 24, "bold")
    lines.append(
        f'<rect x="{box_x}" y="{box_y}" width="{box_width}" height="{box_height}" '
        f'fill="#fafafa" stroke="#cfcfcf" stroke-width="2"/>'
    )
    img_x, img_y, img_width, img_height = layout_scale(
        src_width, src_height, box_x + 20, box_y + 20, box_width - 40, box_height - 40
    )
    lines.append(
        f'<image href="{relative_image}" x="{img_x}" y="{img_y}" '
        f'width="{img_width}" height="{img_height}" preserveAspectRatio="xMidYMid meet"/>'
    )


def build_sheet(output_path: Path, title: str, subtitle: str, panels: Iterable[dict[str, object]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = svg_header()
    add_text(lines, 120, 115, title, 40, "bold")
    add_text(lines, 120, 155, subtitle, 20)
    for panel in panels:
        add_panel(
            lines,
            output_path,
            str(panel["file"]),
            str(panel["label"]),
            str(panel["title"]),
            int(panel["x"]),
            int(panel["y"]),
            int(panel["w"]),
            int(panel["h"]),
        )
    lines.extend(svg_footer())
    output_path.write_text("\n".join(lines), encoding="utf-8")


impactful = [
    {
        "tier": "impactful_dense",
        "tier_dir": "01_impactful_dense",
        "group": "main_maturation_summary",
        "group_dir": "",
        "notebook": "3_TG_incoming_sample_QC_maturation.qmd",
        "source": "3_TG_incoming_sample_QC_maturation_files/figure-html/fig-composite-signature-trajectory-1.png",
        "curated_name": "01_signature_program_trajectory.png",
        "title": "Composite signature program trajectory",
        "use": "Main figure panel candidate",
        "notes": "Strong, information-dense weekly maturation summary across programs.",
    },
    {
        "tier": "impactful_dense",
        "tier_dir": "01_impactful_dense",
        "group": "main_maturation_summary",
        "group_dir": "",
        "notebook": "3_TG_incoming_sample_QC_maturation.qmd",
        "source": "3_TG_incoming_sample_QC_maturation_files/figure-html/fig-pca-projection-trajectory-1.png",
        "curated_name": "02_pca_identity_trajectory.png",
        "title": "PCA identity projection trajectory",
        "use": "Main figure panel candidate",
        "notes": "Adds geometric maturation trajectory in fixed TG reference space.",
    },
]


supplemental = [
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "reference_validation",
        "group_dir": "01_reference_validation",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-sample-correlation-heatmap-1.png",
        "curated_name": "01_sample_correlation_heatmap.png",
        "title": "Pseudobulk sample correlation heatmap",
        "use": "Supplement 1",
        "notes": "Reference design validation; strong reviewer-facing support.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "reference_validation",
        "group_dir": "01_reference_validation",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-variance-partition-1.png",
        "curated_name": "02_variance_partition.png",
        "title": "Variance explained by donor vs identity",
        "use": "Supplement 1",
        "notes": "Shows identity dominates donor effects in the reference.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "reference_validation",
        "group_dir": "01_reference_validation",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-library-complexity-saturation-1.png",
        "curated_name": "03_library_complexity_saturation.png",
        "title": "Library complexity saturation",
        "use": "Supplement 1",
        "notes": "Supports the 20-nuclei floor and pseudobulk stability claims.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "reference_validation",
        "group_dir": "01_reference_validation",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-bulk-reference-projection-1.png",
        "curated_name": "04_bulk_reference_projection.png",
        "title": "Bulk projection onto the TG reference",
        "use": "Supplement 1",
        "notes": "Connects incoming bulk samples to the reference manifold.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "reference_validation",
        "group_dir": "01_reference_validation",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-bulk-to-reference-correlation-1.png",
        "curated_name": "05_bulk_reference_correlation.png",
        "title": "Bulk similarity to TG reference identities",
        "use": "Supplement 1",
        "notes": "Good supporting view for identity assignment consistency.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "signature_diagnostics",
        "group_dir": "02_signature_diagnostics",
        "notebook": "4_TG_deconvolution_diagnostics.qmd",
        "source": "4_TG_deconvolution_diagnostics_files/figure-html/fig-signature-collinearity-1.png",
        "curated_name": "01_signature_collinearity.png",
        "title": "Signature matrix collinearity audit",
        "use": "Supplement 2",
        "notes": "Core mechanistic explanation for the initial NNLS collapse.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "signature_diagnostics",
        "group_dir": "02_signature_diagnostics",
        "notebook": "4_TG_deconvolution_diagnostics.qmd",
        "source": "4_TG_deconvolution_diagnostics_files/figure-html/synthetic-recovery-test-1.png",
        "curated_name": "02_synthetic_recovery.png",
        "title": "Synthetic mixture recovery test",
        "use": "Supplement 2",
        "notes": "Demonstrates identifiability failure in a controlled benchmark.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "signature_diagnostics",
        "group_dir": "02_signature_diagnostics",
        "notebook": "4_TG_deconvolution_diagnostics.qmd",
        "source": "4_TG_deconvolution_diagnostics_files/figure-html/ridge-deconvolution-comparison-1.png",
        "curated_name": "03_ridge_solver_comparison.png",
        "title": "Alternative solver cross-check",
        "use": "Supplement 2",
        "notes": "Shows the issue is not specific to the NNLS solver.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "signature_diagnostics",
        "group_dir": "02_signature_diagnostics",
        "notebook": "4_TG_deconvolution_diagnostics.qmd",
        "source": "4_TG_deconvolution_diagnostics_files/figure-html/fig-marker-specificity-rescreen-1.png",
        "curated_name": "04_marker_specificity_rescreen.png",
        "title": "Marker specificity rescreen",
        "use": "Supplement 2",
        "notes": "Dense but useful support for the decision to expand the marker panel.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "signature_diagnostics",
        "group_dir": "02_signature_diagnostics",
        "notebook": "5_TG_single_cell_marker_refinement.qmd",
        "source": "5_TG_single_cell_marker_refinement_files/figure-html/fig-refit-vif-comparison-1.png",
        "curated_name": "05_refit_vif_comparison.png",
        "title": "VIF before and after marker expansion",
        "use": "Supplement 2",
        "notes": "Bridges diagnostics to the corrected expanded manifest.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "weekly_qc_rerun",
        "group_dir": "03_weekly_qc_rerun",
        "notebook": "3_TG_incoming_sample_QC_maturation.qmd",
        "source": "3_TG_incoming_sample_QC_maturation_files/figure-html/fig-aggregate-neuronal-index-1.png",
        "curated_name": "01_aggregate_neuronal_index.png",
        "title": "Aggregate neuronal identity index",
        "use": "Supplement 3",
        "notes": "Compact QC summary that complements the main figure without replacing it.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "weekly_qc_rerun",
        "group_dir": "03_weekly_qc_rerun",
        "notebook": "3_TG_incoming_sample_QC_maturation.qmd",
        "source": "3_TG_incoming_sample_QC_maturation_files/figure-html/fig-deconvolution-stacked-bar-1.png",
        "curated_name": "02_deconvolution_stacked_bar.png",
        "title": "Initial deconvolution stacked composition",
        "use": "Supplement 3",
        "notes": "Documents the pre-expansion collapse state for comparison.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "weekly_qc_rerun",
        "group_dir": "03_weekly_qc_rerun",
        "notebook": "6_TG_deconvolution_rerun_expanded_markers.qmd",
        "source": "6_TG_deconvolution_rerun_expanded_markers_files/figure-html/fig-vif-recheck-1.png",
        "curated_name": "04_vif_recheck.png",
        "title": "VIF recheck in the rerun joint space",
        "use": "Supplement 3",
        "notes": "Shows the expanded panel materially improves identifiability in the final analysis space.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "weekly_qc_rerun",
        "group_dir": "03_weekly_qc_rerun",
        "notebook": "6_TG_deconvolution_rerun_expanded_markers.qmd",
        "source": "6_TG_deconvolution_rerun_expanded_markers_files/figure-html/fig-composition-before-after-1.png",
        "curated_name": "05_composition_before_after.png",
        "title": "Composition before vs after expansion",
        "use": "Supplement 3",
        "notes": "Most direct side-by-side evidence that the rerun resolved subtype collapse.",
    },
    {
        "tier": "supplemental",
        "tier_dir": "02_supplemental",
        "group": "weekly_qc_rerun",
        "group_dir": "03_weekly_qc_rerun",
        "notebook": "6_TG_deconvolution_rerun_expanded_markers.qmd",
        "source": "6_TG_deconvolution_rerun_expanded_markers_files/figure-html/fig-composition-ci-trajectory-1.png",
        "curated_name": "06_composition_ci_trajectory.png",
        "title": "Composition trajectories with bootstrap CIs",
        "use": "Supplement 3",
        "notes": "Good future supplemental candidate, but current render is missing.",
    },
]


excluded = [
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_qc",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/qc-distributions-1.png",
        "curated_name": "",
        "title": "QC distributions",
        "use": "Notebook only",
        "notes": "Methods-heavy QC context; too dense for the polished set.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_qc",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/qc-feature-count-relationship-1.png",
        "curated_name": "",
        "title": "QC feature-count relationship",
        "use": "Notebook only",
        "notes": "Useful for auditability, not a strong publication-facing figure.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/neuron-sct-umap-1.png",
        "curated_name": "",
        "title": "Neuron SCT UMAP",
        "use": "Notebook only",
        "notes": "Helpful exploration figure, but not central to the weekly QC story.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/neuron-umap-cluster-1.png",
        "curated_name": "",
        "title": "Neuron UMAP clusters",
        "use": "Notebook only",
        "notes": "Reference-construction detail rather than a polished summary figure.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/neuron-marker-plots-1.png",
        "curated_name": "",
        "title": "Neuron marker plots",
        "use": "Notebook only",
        "notes": "Dense exploratory marker panel; better left in the notebook.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/marker-violin-plots-1.png",
        "curated_name": "",
        "title": "Marker violin plots",
        "use": "Notebook only",
        "notes": "Dense panel better used as supporting exploration, not final curation.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/pseudobulk-marker-heatmap-1.png",
        "curated_name": "",
        "title": "Neuronal pseudobulk marker heatmap",
        "use": "Notebook only",
        "notes": "Scientifically useful but too large and dense for the polished subset.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/tg-broad-purity-marker-heatmap-1.png",
        "curated_name": "",
        "title": "Broad TG purity marker heatmap",
        "use": "Notebook only",
        "notes": "Broad reference support that is redundant with cleaner reviewer-facing figures.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "single_cell_reference_build",
        "group_dir": "",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/combined-pseudobulk-marker-heatmap-1.png",
        "curated_name": "",
        "title": "Integrated marker heatmap",
        "use": "Notebook only",
        "notes": "Too information-dense to fit cleanly into the current polished package.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "weekly_qc_methods",
        "group_dir": "",
        "notebook": "3_TG_incoming_sample_QC_maturation.qmd",
        "source": "3_TG_incoming_sample_QC_maturation_files/figure-html/incoming-qc-thresholds-1.png",
        "curated_name": "",
        "title": "Incoming QC thresholds",
        "use": "Notebook only",
        "notes": "Important operational QC, but not part of the stronger narrative figure set.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "weekly_qc_methods",
        "group_dir": "",
        "notebook": "3_TG_incoming_sample_QC_maturation.qmd",
        "source": "3_TG_incoming_sample_QC_maturation_files/figure-html/incoming-reference-correlation-1.png",
        "curated_name": "",
        "title": "Incoming sample similarity to reference identities",
        "use": "Notebook only",
        "notes": "Good audit figure, but secondary to the curated supplemental selection.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "diagnostics",
        "group_dir": "",
        "notebook": "5_TG_single_cell_marker_refinement.qmd",
        "source": "5_TG_single_cell_marker_refinement_files/figure-html/fig-cluster-purity-1.png",
        "curated_name": "",
        "title": "Cluster purity",
        "use": "Notebook only",
        "notes": "Useful for internal confidence, but not necessary in the polished set.",
    },
    {
        "tier": "exclude_notebook_only",
        "tier_dir": "03_exclude_notebook_only",
        "group": "diagnostics",
        "group_dir": "",
        "notebook": "5_TG_single_cell_marker_refinement.qmd",
        "source": "5_TG_single_cell_marker_refinement_files/figure-html/fig-effective-rank-1.png",
        "curated_name": "",
        "title": "Effective rank of neuronal block",
        "use": "Notebook only",
        "notes": "Strong analytic support, but less immediately interpretable than the chosen supplemental figures.",
    },
]


all_entries = [copy_panel(entry.copy()) for entry in impactful + supplemental if entry["curated_name"]]
all_entries.extend(excluded)


bulk_first_story = [
    {
        "story_order": 1,
        "stage": "bulk_to_reference",
        "tier_dir": "04_bulk_first_story",
        "group_dir": "01_bulk_to_reference",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-bulk-reference-projection-1.png",
        "curated_name": "01_bulk_reference_projection.png",
        "title": "Bulk TG samples projected onto the integrated TG reference",
        "role": "Main flow",
        "notes": "Lead with bulk RNA-seq; this is the clearest entry point into the reference comparison.",
    },
    {
        "story_order": 2,
        "stage": "bulk_to_reference",
        "tier_dir": "04_bulk_first_story",
        "group_dir": "01_bulk_to_reference",
        "notebook": "2_TG_reviewer_figures.qmd",
        "source": "2_TG_reviewer_figures_files/figure-html/fig-bulk-to-reference-correlation-1.png",
        "curated_name": "02_bulk_reference_correlation.png",
        "title": "Bulk sample similarity to single-cell-derived reference identities",
        "role": "Main flow",
        "notes": "Makes the bulk-to-reference identity matches explicit before moving into the reference build.",
    },
    {
        "story_order": 3,
        "stage": "processed_single_cell_reference",
        "tier_dir": "04_bulk_first_story",
        "group_dir": "02_processed_single_cell_reference",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/combined-pseudobulk-pca-1.png",
        "curated_name": "01_integrated_reference_pca.png",
        "title": "Processed single-cell pseudobulk reference structure",
        "role": "Main flow",
        "notes": "Shows how the processed single-cell data were condensed into the integrated TG reference used for bulk interpretation.",
    },
    {
        "story_order": 4,
        "stage": "signature_programs",
        "tier_dir": "04_bulk_first_story",
        "group_dir": "03_signature_programs",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/tg-broad-purity-marker-heatmap-1.png",
        "curated_name": "01_broad_type_marker_programs.png",
        "title": "Broad single-cell marker programs with neuron and satellite glia contrast",
        "role": "Main flow",
        "notes": "Best broad signature view for the two dominant interpretive compartments: neuron and satellite glia.",
    },
    {
        "story_order": 5,
        "stage": "signature_programs",
        "tier_dir": "04_bulk_first_story",
        "group_dir": "03_signature_programs",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/pseudobulk-marker-heatmap-1.png",
        "curated_name": "02_neuronal_subtype_marker_heatmap.png",
        "title": "Neuronal subtype signatures in the processed single-cell reference",
        "role": "Supplemental flow",
        "notes": "Deeper signature detail within the neuronal compartment after the broad neuron-vs-satellite-glia split is established.",
    },
    {
        "story_order": 6,
        "stage": "signature_programs",
        "tier_dir": "04_bulk_first_story",
        "group_dir": "03_signature_programs",
        "notebook": "1_TG_sc_reprocessing.qmd",
        "source": "1_TG_sc_reprocessing_files/figure-html/combined-pseudobulk-marker-heatmap-1.png",
        "curated_name": "03_integrated_marker_heatmap.png",
        "title": "Integrated neuronal and non-neuronal marker heatmap",
        "role": "Supplemental flow",
        "notes": "Most complete reference signature panel, but denser than the preferred broad marker program view.",
    },
]

story_entries = [copy_panel(entry.copy()) for entry in bulk_first_story]


OUT_DIR.mkdir(parents=True, exist_ok=True)
(OUT_DIR / "03_exclude_notebook_only").mkdir(parents=True, exist_ok=True)


with (OUT_DIR / "figure_tier_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "tier",
            "group",
            "notebook",
            "title",
            "source",
            "curated_path",
            "use",
            "status",
            "notes",
        ],
    )
    writer.writeheader()
    for entry in all_entries:
        writer.writerow(
            {
                "tier": entry["tier"],
                "group": entry["group"],
                "notebook": entry["notebook"],
                "title": entry["title"],
                "source": entry["source"],
                "curated_path": entry.get("curated_path", ""),
                "use": entry["use"],
                "status": entry.get("status", "excluded"),
                "notes": entry["notes"],
            }
        )


build_sheet(
    OUT_DIR / "01_impactful_dense" / "Figure_01_impactful_dense_draft.svg",
    "Figure 1 draft - impactful and dense",
    "8.5x11 target layout centered on the two strongest weekly maturation panels.",
    [
        {
            "file": "01_signature_program_trajectory.png",
            "label": "A",
            "title": "Composite signature program trajectory",
            "x": 140,
            "y": 250,
            "w": 2270,
            "h": 1080,
        },
        {
            "file": "02_pca_identity_trajectory.png",
            "label": "B",
            "title": "Neuron-only PCA identity projection trajectory",
            "x": 140,
            "y": 1520,
            "w": 2270,
            "h": 1380,
        },
    ],
)


build_sheet(
    OUT_DIR
    / "02_supplemental"
    / "01_reference_validation"
    / "Supplement_01_reference_validation_draft.svg",
    "Supplement 1 draft - reference validation",
    "Reviewer-facing support for the TG pseudobulk reference and bulk-to-reference mapping.",
    [
        {
            "file": "01_sample_correlation_heatmap.png",
            "label": "A",
            "title": "Sample correlation heatmap",
            "x": 120,
            "y": 250,
            "w": 1100,
            "h": 1100,
        },
        {
            "file": "02_variance_partition.png",
            "label": "B",
            "title": "Variance partition",
            "x": 1320,
            "y": 250,
            "w": 1100,
            "h": 720,
        },
        {
            "file": "03_library_complexity_saturation.png",
            "label": "C",
            "title": "Library complexity saturation",
            "x": 1320,
            "y": 860,
            "w": 1100,
            "h": 490,
        },
        {
            "file": "04_bulk_reference_projection.png",
            "label": "D",
            "title": "Bulk reference projection",
            "x": 120,
            "y": 1520,
            "w": 1100,
            "h": 900,
        },
        {
            "file": "05_bulk_reference_correlation.png",
            "label": "E",
            "title": "Bulk-to-reference correlation",
            "x": 1320,
            "y": 1520,
            "w": 1100,
            "h": 900,
        },
    ],
)


with (OUT_DIR / "04_bulk_first_story" / "narrative_flow_manifest.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "story_order",
            "stage",
            "notebook",
            "title",
            "source",
            "curated_path",
            "role",
            "status",
            "notes",
        ],
    )
    writer.writeheader()
    for entry in story_entries:
        writer.writerow(
            {
                "story_order": entry["story_order"],
                "stage": entry["stage"],
                "notebook": entry["notebook"],
                "title": entry["title"],
                "source": entry["source"],
                "curated_path": entry.get("curated_path", ""),
                "role": entry["role"],
                "status": entry.get("status", "missing"),
                "notes": entry["notes"],
            }
        )


build_sheet(
    OUT_DIR / "04_bulk_first_story" / "Figure_02_bulk_first_story_draft.svg",
    "Figure 2 draft - bulk RNA-seq to single-cell signature flow",
    "Start with bulk TG RNA-seq, then reveal the processed single-cell reference and the neuron-vs-satellite-glia signature split.",
    [
        {
            "file": "01_bulk_to_reference/01_bulk_reference_projection.png",
            "label": "A",
            "title": "Bulk TG projection onto the reference",
            "x": 120,
            "y": 250,
            "w": 2300,
            "h": 820,
        },
        {
            "file": "01_bulk_to_reference/02_bulk_reference_correlation.png",
            "label": "B",
            "title": "Bulk-to-reference identity correlation",
            "x": 120,
            "y": 1260,
            "w": 1100,
            "h": 980,
        },
        {
            "file": "02_processed_single_cell_reference/01_integrated_reference_pca.png",
            "label": "C",
            "title": "Processed single-cell reference PCA",
            "x": 1320,
            "y": 1260,
            "w": 1100,
            "h": 980,
        },
        {
            "file": "03_signature_programs/01_broad_type_marker_programs.png",
            "label": "D",
            "title": "Broad marker programs: neuron vs satellite glia",
            "x": 120,
            "y": 2440,
            "w": 2300,
            "h": 700,
        },
    ],
)


build_sheet(
    OUT_DIR / "04_bulk_first_story" / "Supplement_bulk_first_signature_detail_draft.svg",
    "Bulk-first supplement draft - signature detail",
    "Deeper single-cell signature detail after the main bulk-to-reference narrative is established.",
    [
        {
            "file": "03_signature_programs/02_neuronal_subtype_marker_heatmap.png",
            "label": "A",
            "title": "Neuronal subtype signature heatmap",
            "x": 120,
            "y": 250,
            "w": 2300,
            "h": 1220,
        },
        {
            "file": "03_signature_programs/03_integrated_marker_heatmap.png",
            "label": "B",
            "title": "Integrated neuronal and non-neuronal marker heatmap",
            "x": 120,
            "y": 1650,
            "w": 2300,
            "h": 1400,
        },
    ],
)


build_sheet(
    OUT_DIR
    / "02_supplemental"
    / "02_signature_diagnostics"
    / "Supplement_02_signature_diagnostics_draft.svg",
    "Supplement 2 draft - signature diagnostics",
    "Mechanistic support for why the original signature matrix collapsed and why the expanded panel helps.",
    [
        {
            "file": "01_signature_collinearity.png",
            "label": "A",
            "title": "Signature collinearity audit",
            "x": 120,
            "y": 250,
            "w": 1100,
            "h": 950,
        },
        {
            "file": "02_synthetic_recovery.png",
            "label": "B",
            "title": "Synthetic mixture recovery",
            "x": 1320,
            "y": 250,
            "w": 1100,
            "h": 950,
        },
        {
            "file": "03_ridge_solver_comparison.png",
            "label": "C",
            "title": "Alternative solver cross-check",
            "x": 120,
            "y": 1420,
            "w": 1100,
            "h": 700,
        },
        {
            "file": "05_refit_vif_comparison.png",
            "label": "D",
            "title": "VIF before vs after expansion",
            "x": 1320,
            "y": 1420,
            "w": 1100,
            "h": 700,
        },
        {
            "file": "04_marker_specificity_rescreen.png",
            "label": "E",
            "title": "Marker specificity rescreen",
            "x": 120,
            "y": 2350,
            "w": 2300,
            "h": 700,
        },
    ],
)


build_sheet(
    OUT_DIR
    / "02_supplemental"
    / "03_weekly_qc_rerun"
    / "Supplement_03_weekly_qc_rerun_draft.svg",
    "Supplement 3 draft - weekly QC and rerun detail",
    "Context panels showing aggregate tracking, initial composition, and the refined corrected rerun.",
    [
        {
            "file": "01_aggregate_neuronal_index.png",
            "label": "A",
            "title": "Aggregate neuronal index",
            "x": 120,
            "y": 250,
            "w": 1100,
            "h": 640,
        },
        {
            "file": "04_vif_recheck.png",
            "label": "B",
            "title": "VIF recheck after expansion",
            "x": 1320,
            "y": 250,
            "w": 1100,
            "h": 640,
        },
        {
            "file": "02_deconvolution_stacked_bar.png",
            "label": "C",
            "title": "Initial deconvolution composition",
            "x": 120,
            "y": 1080,
            "w": 2300,
            "h": 720,
        },
        {
            "file": "05_composition_before_after.png",
            "label": "D",
            "title": "Composition before vs after expansion",
            "x": 120,
            "y": 1990,
            "w": 2300,
            "h": 850,
        },
    ],
)
