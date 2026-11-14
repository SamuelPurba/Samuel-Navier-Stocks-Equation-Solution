# Scientific & Academic Documentation Standards (IEEE / Scopus Q1)

When authoring or updating scientific documentation, mathematical monographs, or `README.md` files in this workspace, strictly enforce the following four pillars:

## 1. Human-Centric Readability & Pedagogical Intuition ("100% Mudah Dibaca Manusia")
- **Intuition Before Formalism:** Always precede dense mathematical derivations with a plain-language executive summary (*Inti Sari*), physical/engineering intuition, and clean ASCII/Mermaid visual schematics.
- **Readable Typography for Compound Fractions:** Avoid 3- or 4-story vertically stacked `\frac` expressions that shrink inner text to microscopic sizes; write compound ratios as side-by-side products (`\frac{A}{B} \times \left(\frac{C}{D}\right)`) while preserving 100% of the manuscript's algebraic steps and variants.
- **Explicit Gap Tagging:** Never invent unstated physical definitions or empirical datasets; tag any undefined manuscript symbol or unit explicitly with `[Repository Gap — Author Clarification Required]`.

## 2. 100% Markdown, LaTeX, & Diagram Rendering Hygiene ("100% Benar")
- **Display Math (`$$`):** Place opening and closing `$$` delimiters on their own dedicated lines.
- **Safe Delimiters (`\left` / `\right`):** Never use `\left\{`, `\right\}`, `\left\lbrace`, or `\right\rbrace` in Markdown math, and never place `\left` / `\right` inside integral limits (`_{...}^{...}`). Use `\left[` / `\right]` and `\left(` / `\right)` exclusively.
- **Plain-Text Headings & TOC Slugs:** Never place `$ ... $` math delimiters inside `#`, `##`, or `###` headings or Table of Contents links; use clean Unicode symbols (`π`, `σ`, `ΔP`, `ΔL`) in headings and verify that every `#anchor` link resolves to an exact GitHub heading slug.
- **Markdown Tables:** Never use raw pipe characters (`|` or `\|`) inside `$ ... $` within Markdown table cells; always use `\Vert`.
- **Mermaid Diagrams:** Use `flowchart TD` or `flowchart LR` with double-quoted labels free of `<br>`, `{`, `}`, `|`, or `>`.

## 3. Zero-Dependency Reproducible Verification ("100% Works")
- **Standard-Library Python Solvers:** Embedded or standalone numerical verification scripts (`solver.py`, `test_solver.py`) must use 100% Python 3 standard library (`math`, `unittest`, `json`, `argparse`) with zero external dependencies.
- **Automated Verification Gate:** Before concluding any documentation or solver update, run `python -m unittest discover -v` and `python solver.py --verify` to verify `0.0` residual error, bitwise code parity, and zero rendering/anchor defects.

## 4. Canonical Author Identity, Core Equations & Architecture Invariants
- **Author & Affiliation:** Always credit **Samuel Hasiholan Omega Purba, S. Tr. T.** (*Department of Electrical Engineering / Teknik Elektro, Study Program in Robotics and Artificial Intelligence Engineering / Prodi Teknik Robotika dan Kecerdasan buatan, Politeknik Negeri Batam*).
- **Required IEEE Section Hierarchy:**
  1. Title, Academic Header & Verified Badges
  2. Abstract (150–250 words) & *Index Terms*
  3. Nomenclature Table (with explicit `[Repository Gap — Author Clarification Required]` tags for unstated units/definitions)
  4. *I. Introduction* (Motivation, Geometric Schematic, Spectral Cutoff Intuition, Scope)
  5. *II. Mathematical Framework and Step-by-Step Derivation* (Cross-check of $\pi = 4 - \frac{\Delta}{a}$ vs. $\pi = \frac{4 - \Delta}{a}$ and $\frac{2}{30} \to \frac{1}{15}$ vs. $\frac{8}{30} \to \frac{4}{15}$, Core 6-Step Sequence, and Complete `.docx`/`.pdf` Derivations)
  6. *III. Interactive Command Protocols* (`/Team-Viewer`, `/learn`, `/humanize` specification table, Mermaid workflow, and example prompts)
  7. *IV. Repository Structure and Reproducibility* (Directory tree, zero-dependency `solver.py`, and verified console output)
  8. *V. Limitations, Open Questions and Future Work*
  9. *VI. Citation and References* (IEEE reference entry & BibTeX block) + *Repository Audit Log*
