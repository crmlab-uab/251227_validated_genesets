# Copilot-instructions systematic review — 2025-12-28

Repos reviewed:
- /data/bRNA3F: `.github/copilot-instructions.md` (template, authoritative)
- /data/251207_mNSC: `.github/copilot-instructions.md` (copied from template)
- /data/251227_validated_genesets: `.github/copilot-instructions.md` (copied; contains nested fence artifacts)
- /data/250829_rstudio: NO `.github/copilot-instructions.md` found

High-level findings
- Core rules are present and detailed in the `bRNA3F` template (multi-repo safety, git pattern, config architecture, render caveats, timestamps). This is the authoritative source.
- `251207_mNSC` contains a full copy of the template (consistent with the sync pattern).
- `251227_validated_genesets` contains a copy but it has an extra wrapper fence and is abridged (`... (copied from template repo)`), indicating an incomplete copy-paste.
- `250829_rstudio` is missing a copilot-instructions file entirely.

Inconsistencies & issues
- Missing file in `250829_rstudio` — repository lacks agent guidance.
- `251227_validated_genesets` copy contains duplicated/backtick wrapper markers (`````instructions) and an ellipsis placeholder; this can confuse automated readers and is non-authoritative.
- Slightly different scopes: the template contains many project-specific rules (render behavior, chunk numbering). Analysis repos should inherit rules but not duplicate large sections that become stale.
- No automated verification that the analysis repos keep the template file in sync (recommend md5 checksum workflow described in template be applied to these files as well).

Recommended harmonization
1. Keep `/data/bRNA3F/.github/copilot-instructions.md` as the canonical source. Do not edit analysis copies directly — instead, update the template and propagate.
2. Replace analysis repo copies with a small header that points to the canonical template and includes a one-line local override section (if needed). Example header:

   **Header (analysis repo)**
   - Purpose: short repo-specific guidance (1–3 bullets)
   - Canonical instructions: reference to template at `/data/bRNA3F/.github/copilot-instructions.md`

3. For `251227_validated_genesets`: replace existing `.github/copilot-instructions.md` with a clean header (no duplicated fences) and a single-line pointer to the canonical file.
4. For `250829_rstudio`: add a minimal `.github/copilot-instructions.md` header that points to the canonical template.
5. Add a lightweight validation script (in `scripts/` or CI) that asserts:
   - Each analysis repo has a `.github/copilot-instructions.md` file
   - The file contains a canonical-hash line matching the template (use `md5sum`), or simply contains the expected pointer phrase.

Proposed files to change
- `/data/251227_validated_genesets/.github/copilot-instructions.md` — replace with header + pointer
- `/data/250829_rstudio/.github/copilot-instructions.md` — add header + pointer
- Optionally add `/data/bRNA3F/scripts/verify_copilot_instructions.sh` to check repository coverage and write checksums into `docs/sessions/`.

Suggested minimal header template (for analysis repos)
```
# Copilot instructions — analysis repo
This repository uses the central copilot instructions maintained in the bRNA3F template repo.
See: /data/bRNA3F/.github/copilot-instructions.md

Local overrides: (none)
```

Next steps I can take
- (A) Patch `251227_validated_genesets` and add the header file.
- (B) Add the header file to `250829_rstudio`.
- (C) Add a verification script and run it across the four repos, write an audit report into `docs/`.

Which actions should I perform now? (I can do A+B+C in sequence.)
