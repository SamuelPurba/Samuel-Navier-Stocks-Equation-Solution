# Scientific & Academic Documentation Standards (IEEE / Scopus Q1)

When authoring or updating scientific documentation, mathematical monographs, or `README.md` files in this workspace, strictly enforce the following four pillars:

## 1. Human-Centric Readability & Pedagogical Intuition ("100% Mudah Dibaca Manusia")
- **Intuition Before Formalism:** Always precede dense mathematical derivations with a plain-language executive summary (*Inti Sari*), physical/engineering intuition, and clean ASCII/Mermaid visual schematics.
- **Dual-Layer Notation (Shorthand + Verbatim Manuscript Form):** Always provide a Quick Notation Guide ($\Delta P \equiv \Delta_{(P_{(Square)} - P_{(Circle)})}$, $I_\zeta \equiv \int_{0}^{\zeta(s) = -\frac{1}{12}} r^4 f(r)\,dr$) and an at-a-glance summary table before multi-step derivations, and pair each step with:
  1. **Penjelasan Mudah** (1-sentence plain-language explanation),
  2. **Bentuk Ringkas** (compact shorthand formula for instant visual comprehension), and
  3. **Bentuk Lengkap Manuskrip** (full verbatim equation in display math `$$`).
- **Readable Typography for Compound Fractions:** Avoid 3- or 4-story vertically stacked `\frac` expressions that shrink inner text to microscopic sizes; write compound ratios as side-by-side products (`\frac{A}{B} \times \left(\frac{C}{D}\right)`) while preserving 100% of the manuscript's algebraic steps and variants.
- **Explicit Gap Tagging:** Never invent unstated physical definitions or empirical datasets; tag any undefined manuscript symbol or unit explicitly with `[Repository Gap — Author Clarification Required]`.

## 2. 100% Markdown, LaTeX, & Diagram Rendering Hygiene ("100% Benar")
- **Display Math (`$$`):** Place opening and closing `$$` delimiters on their own dedicated lines.
- **Safe Delimiters (`\left` / `\right`):** Never use `\left\{`, `\right\}`, `\left\lbrace`, or `\right\rbrace` in Markdown math, and never place `\left` / `\right` inside integral limits (`_{...}^{...}`). Use `\left[` / `\right]` and `\left(` / `\right)` exclusively.
- **Plain-Text Headings & TOC Slugs:** Never place `$ ... $` math delimiters inside `#`, `##`, or `###` headings or Table of Contents links; use clean Unicode symbols (`π`, `σ`, `ΔP`, `ΔL`) in headings and verify that every `#anchor` link resolves to an exact GitHub heading slug.
- **Markdown Tables:** Never use raw pipe characters (`|` or `\|`) inside `$ ... $` within Markdown table cells; always use `\Vert`.
- **Mermaid Diagrams:** Use `flowchart TD` or `flowchart LR` with double-quoted labels free of `<br>`, `{`, `}`, `|`, or `>`.
- **Inline Math & Emphasis Isolation:** Never wrap inline `$ ... $` expressions inside Markdown emphasis delimiters (`*...$...$...*` or `**...$...$...**`), and avoid `\text{...}` containing spaces inside inline `$ ... $` within Markdown tables (write e.g. `$a = 101325.0$ Pa` or `$m^2$` instead of `*(... $a = 101325.0\text{ Pa}$)*`). Keep punctuation, prose units, and bold/italic markers outside `$ ... $` so CommonMark/KaTeX never leak literal `$` delimiters.

## 3. Zero-Dependency Reproducible Verification ("100% Works")
- **Standard-Library Python Solvers vs. Optional Visual Generators:** Core numerical verification scripts (`solver.py`, `test_solver.py`) must remain 100% Python 3 standard library (`math`, `unittest`, `json`, `argparse`) with zero external dependencies. Optional media/simulation scripts (`generate_simulation.py` producing `Samuel-Navier_Stokes_Simulation.mp4` via Pillow + `ffmpeg`) may read from `solver.py` but must never introduce third-party dependencies into `solver.py` or `test_solver.py`.
- **Automated Verification Gate:** Before concluding any documentation or solver update, run `python -m unittest discover -v`, `python solver.py --verify`, and an inline-math/emphasis scan confirming `0.0` residual error, bitwise code parity, `0` emphasis-wrapped or `\text`-containing inline `$ ... $` spans, and zero anchor/Mermaid defects.

## 4. Canonical Author Identity, Core Equations & Architecture Invariants
- **Author & Affiliation:** Always credit **Samuel Hasiholan Omega Purba, S. Tr. T.** (*Department of Electrical Engineering / Teknik Elektro, Study Program in Robotics and Artificial Intelligence Engineering / Prodi Teknik Robotika dan Kecerdasan buatan, Politeknik Negeri Batam*).
- **Formal Title, Repo Alignment & Minimal Badge Strip:** Standardize the formal H1 repository title as **Navier–Stokes Equations Solution** while preserving both the canonical manuscript title (*"Samuel – Navier Stocks Equation Solution"*) and active repository title (*"Samuel–Navier Stocks Equation Solution"*) in the bibliographic subtitle/notes. Keep the header badge strip minimal (`Standard`, `Manuscript`, `Python 3 Audit`), omitting noisy `Unit Tests` and `Repository Reconciled` badges unless explicitly requested.
- **Tri-Domain Manuscript Formulation & Integrity Audit (`P01–P52`):**
  - Document **Branch A** ($\Delta_{(P_{(Square)} - P_{(Circle)})} = (4 - \pi) \times a$, Steps 1–7, `P07–P29`), **Branch B** ($\Delta_{(L_{(Square)} - L_{(Circle)})} = (1 - \pi \times 4^{-1}) \times m^2$, Steps 1–7, `P31–P45`), and the newly added **Einstein Tensor Extension (`P46–P52`)** ($G_{\mu\nu} = 8\pi G T_{\mu\nu} \mapsto G_{\mu\nu} = 30 \times \sigma \times I_\zeta^{-1} G T_{\mu\nu}$ via $\pi = (4^{-1} \times 15 \times \sigma) \times I_\zeta^{-1}$) cleanly in display math (`$$`), cross-checked against the live 52-paragraph `.docx` and 5-page `.pdf`.
  - **PDF/OCR Sanitization:** Exclude PDF text-extraction artifacts (stray $\theta$ symbols, malformed angle-bracket $\pi$ fractions, `f_0` instead of `\int_0`, and garbled subscripts) in favor of clean `.docx` OMML notation.
  - **Factor-of-$4^{-1}$ & Factor-of-$4$ Erratum Audit (Branch B) & Tensor Coupling Gap (`P46–P52`):** Explicitly examine Branch B Step 1 ($\Delta L = (1 - \pi \times 4^{-1}) \times m^2$) against Step 3 ($\pi = \frac{m^2 - \Delta L}{m^2}$ vs. strict $4 \times \frac{m^2 - \Delta L}{m^2}$) and Step 6–7 ($\pi$ vs. $\pi \times 4^{-1}$), and transparently flag `[Repository Gap — Author Clarification Required]` on unstated physical details in `P46–P52` (e.g., natural units $c=1$ vs. omitted $c^{-4}$, and unstated Navier–Stokes-to-Einstein PDE coupling) without altering source equations.
- **Required IEEE Section Hierarchy & Closing Statements:**
  1. Title, Academic Header & Minimal Verified Badges
  2. Abstract (150–250 words) & *Index Terms*
  3. Nomenclature Table (with explicit `[Repository Gap — Author Clarification Required]` tags for unstated units/definitions)
  4. *I. Introduction* (Motivation, Geometric Schematic, Spectral Cutoff Intuition, Scope)
  5. *II. Mathematical Framework and Step-by-Step Derivation* (Cross-check of $\pi = 4 - \frac{\Delta}{a}$ vs. $\pi = \frac{4 - \Delta}{a}$ and $\frac{2}{30} \to \frac{1}{15}$ vs. $\frac{8}{30} \to \frac{4}{15}$, Core Sequence, Branch A, Branch B, Einstein Tensor Extension `P46–P52`, and Mathematical Audit & Verification)
  6. *III. Interactive Command Protocols* (`/Team-Viewer`, `/learn`, `/humanize` specification table, Mermaid workflow, and example prompts)
  7. *IV. Repository Structure and Reproducibility* (Directory tree covering all live workspace files including `generate_simulation.py` and `Samuel-Navier_Stokes_Simulation.mp4`, zero-dependency `solver.py`, and verified console output)
  8. *V. Limitations, Open Questions and Future Work*
  9. *VI. Citation and References* (IEEE reference entry & BibTeX block) + preserve the author's **`## 📜 Lisensi & Hak Cipta Statement`** (MIT License, Hak Cipta © 2026 Samuel Hasiholan Omega Purba / Samuel Hasiholan Omega, S. Tr. T.). When a user prompt includes an explicit Phase 5 output contract (single ` ```markdown ... ``` ` block + `(a)(b)(c)` Grounding & Audit Log), strictly satisfy both the repository `README.md` update and the requested response output contract.
