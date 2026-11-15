# Presentation update plan (after report 13, "labels as measurements")

Status: PLAN v1, 2026-11-15. Source of numbers: reports/13_labels_as_measures.md and the CSVs it cites (data/derived). Deck source: docs/presentation_deck/project/
(64 slides, plain white content slides, navy dividers, Arial, tags such as "Q5 · METHOD"). The style, CSS classes and slide templates stay unchanged.
Rules: the deck stays self-contained (no references to earlier versions, reviews or corrections); one message per slide; method before result; short bullets; figures over text;
legends outside the plot; every number must match reports/13 and the tables; claims no stronger than the report.

## A. Edits to existing slides (ids unchanged)
| Slide | Edit |
|---|---|
| `contents` | Q1 answer -> "110 genes: a lower bound; no database holds them all"; Q5 answer -> "Only where a neighbour has at least 30% identity" |
| `q1` (divider) | short answer: "110 genes, a lower bound: no single database holds them all, and 11 more reviewed candidates exist" |
| `q5` (divider) | short answer: "Only where a close neighbour exists: below about 30% identity, transfer is no better than guessing" |
| `q5m1` (method) | card 2 "Exclude family" stays; card 3 baselines stays; add note: "Then repeat with neighbours limited to a maximum identity (100 to 20%)" |
| `function` (result) | title -> "Across families, neighbours beat guessing by a few points that the data cannot confirm"; lead -> "Family-excluded gain: +9 to +11 points, intervals include zero"; small -> "Next: the gain depends on neighbour identity" |
| `hypotheses` | add to third card text: "most need a neighbour of at least 30% identity (see previous results)" ONLY IF the count is computed; plan: compute the number of the 30 genes whose nearest labelled neighbour has identity >= 30% and state it |
| `review` | Weakened column: "Function transfer: 85% → no better than guessing below 30% identity" |
| `next` | Limits: add "110 is a lower bound (11 further reviewed candidates)"; Next: "Add the 11 candidates; experimental-only function labels" |
| `takehome` | add a fifth statement: "Function transfer needs a neighbour with at least about 30% identity." (font 48px to keep five lines in the box) |
| `followups` | unchanged |

## B. New slides (ids, section, type, content, figure)
Q1 (after `sources`/`attrition`/`gaps`, before `models`):
1. `q1m2` METHOD "Add a third source and list what the first two miss": cards: Collect (reviewed human UniProt entries with InterPro IPR029058 or 12 Pfam alpha/beta-hydrolase families); Collapse (to genes as before); Compare (overlap with the 110; candidates are listed, not analysed).
2. `cands` RESULT "The third source holds all 110 genes and 11 more": cards "121 / 78" (IPR029058 / Pfam entries), "11" candidates, "0 of 110" missing; list of 11 gene names as note (ABHD12B, AADACL2, AADACL3, AADACL4, TEX30, FAM135A, FAM135B, TMCO4, TMEM53, ARB2A, TG); caveat: near-superset, so no recall estimate; candidates need structure QC.
Q2 (after `triad`):
3. `q2m3` METHOD "Compare the elbow with a motif and with structural transfer": cards: Elbow (exact call); Motif (legacy pattern at the annotated residue); Transfer (residue aligned to the annotated nucleophile of the nearest other gold protein, leave one out); stratified by evidence of the annotation.
4. `meas` RESULT "Three measures agree on the 31 annotated nucleophiles": cards "31 / 31" elbow, "26 / 31" motif, "31 / 31" transfer; note: experimental subset (n = 7) 7 / 6 / 7; ordering does not change. Caveat: all gold genes are positives and mostly paralogs.
Q4 (after `majority`/`notcompact`, i.e. before the Q5 divider; plus the multiview result):
5. `q4m4` METHOD "Rebuild the domain from 3D distance instead of CATH ranges": cards: Core (200 confident residues nearest the elbow C-alpha, original order); Compare (TM-align and the ABHD tests on cores, K = 150, 180, 200); Test (paired bootstrap of core minus CATH, equivalence margin 0.10).
6. `cores` RESULT "The ABHD structure signal does not depend on CATH trimming": cards "0.78" core vs CATH distance correlation, "0.40 vs 0.40" R-precision core vs CATH, "0.85 vs 0.80" nearest-neighbour rate (+0.05 [0.00, 0.15]); note: cores are centred on the elbow; sequence on cores is fragmentary (median 20 segments).
7. `multiview` RESULT "Across families, structure and sequence agree; EC agrees with neither": figure (bars with CIs, cross-family pairs): structure-sequence 0.50, sequence-FunFam 0.54, structure-FunFam 0.36, structure-core 0.77, EC vs structure -0.13, EC vs sequence 0.03, EC vs FunFam 0.07; side text: EC is coarse (41 of 80 are 3.1.1); the dependent pairs (sequence-FunFam) are not independent confirmation.
Q5 (after `q5m1`/`function`/`hypotheses`):
8. `q5m3` METHOD "Check how well each EC label is supported": cards: Parse (each UniProt reaction block: EC, Rhea, evidence codes); Classify (experimental, curator, similarity, rule, none); Use (restrict queries and neighbours to experimentally supported labels).
9. `prov` RESULT "Most EC labels have some experimental support": cards "58 of 80" some experimental, "51" only experimental, "8" no record; note: 41 of 80 are 3.1.1; restricting to experimental labels changes accuracy by a few points (structure 0.88 vs 0.84 any neighbour; 0.66 vs 0.68 family-excluded).
10. `q5m4` METHOD "Limit how similar the neighbour may be": cards: Identity (global alignment on trimmed domains); Cutoff (maximum identity 100, 40, 35, 30, 25, 20%); Statistic (accuracy minus majority baseline, stem-cluster bootstrap CI).
11. `decay` RESULT "Below 30% identity, no space beats guessing the most common class": figure (difference to baseline vs maximum identity, 4 spaces, CI bars); side: "Structure drops to the random-neighbour level (0.39 to 0.41)"; "Pools are nested: read as a curve".
12. `nnbins` RESULT "Transfer works when the nearest neighbour has at least 30% identity": figure or cards: accuracy by identity of the actual nearest labelled neighbour (structure: <0.25: 0.38, n=13; 0.25-0.30: 1.00, n=5; 0.30-0.40: 0.90, n=10; >0.40: 0.92, n=52; same for sequence/FunFam/elbow window in a small table) as a table (4 spaces x 4 bins) with n.

Total new slides: 12 (deck 76 slides). `deck.json` order: insert as above; sections unchanged.

## C. New figures (scripts/23_presentation_figures.py additions; same palette, Arial, white, legends outside)
- `decay.png` (slide 11): lines with CI bars from `function_decay.csv` mode=all, spaces tmalign/seqid/funfam/anchor_elbow, x = maximum identity (reversed), y = difference to the majority baseline; legend above.
- `multiview.png` (slide 7): horizontal bars with CI, cross-family pairs from `multiview_association.csv`.
- (Slide 12 is an HTML table; no figure.)

## D. Speaker notes and sources
Each new slide has notes with the detail and caveats from reports/13; source lines name the supplementary tables (no repository paths except method slides, as in the current deck).

## E. Acceptance checks
1. Every number on new or edited slides appears in reports/13 or its CSVs (reviewer recomputes at least 15).
2. No new text refers to earlier versions, corrections, reviews or "now".
3. Tags correct (section, METHOD/RESULT); each new result has a preceding method slide in its section; ids unique; deck.json order and sections consistent; 76 slide files.
4. Layout limits (title <= 2 lines at 56 px; figure heights as in existing slides; cards with numbers use xs/xxs classes for 3/4-card rows); figures have legends outside the axes.
5. Claims: no stronger than the report (for example "no space beats guessing below 30% identity" with the quoted intervals; "110 is a lower bound" not "the family has N members").


---
# v2 amendments after independent review (supersede conflicting parts above)
- Numbers: family-excluded gain is "+8 to +11 points" (sequence +0.075); the 20% cutoff is a subset (38 of 80) and is left out of the figure and noted; strata comparisons say "within the intervals".
- Wording: decay title "When neighbours are limited to 30% identity or less, no space beats guessing"; structure falls to the random-neighbour level; gain appears around 35-40% (intervals include zero at 35%).
  function title "Across families the gain over guessing is about ten points, but its interval includes zero". multiview title "Structure and sequence agree moderately; EC distance is not detectably related to any view".
  cores title "The structure-space ABHD signal is unchanged when CATH ranges are replaced by cores" (R-precision 0.40 vs 0.40; nearest-neighbour equivalence not shown; sequence on cores weaker; cores centred on the elbow).
  meas title "The elbow and structural transfer recover all 31; the motif misses five" (sensitivity only; all gold genes are positives). cands: "11 reviewed candidates, not all hydrolases (for example TG)".
- Structure: drop `nnbins` (its table goes into the `decay` notes with n); merge the Q5 methods into one (`q5m3`: label support and neighbour similarity); no `q5m1` note; add `q4m5` (multi-view method); keep `q2m3` + `meas`.
  New slides: q1m2, cands, q2m3, meas, q4m4, cores, q4m5, multiview, q5m3, prov, decay = 11 (deck 75). Q5 order: q5, q5m1, function, q5m3, prov, decay, hypotheses, q5m2, anchored.
- Other edits: question card ("110 reviewed human genes, a lower bound"), gaps card 3 and sources notes (ABHD12B appears in the third source), next notes, takehome notes (five messages), review wording ("roughly 30%"), contents/q1/q4/q5 dividers, hypotheses (7 of the 18 eligible predictions have a nearest labelled neighbour with identity >= 30%: MEST 0.32, PRSS16 0.32, DPP6 0.38, ABHD14A 0.47, EPHX4 0.48, PNLIPRP1 0.67, CES1P1 0.91).
- Text: "bootstrap over paralog groups" replaces "stem-cluster bootstrap". New figures are uploaded as blobs; the PDF is regenerated by the owner after publishing.
- Acceptance additions: each "no better than guessing" statement carries n or the interval; reconcile every hit of "ten points", "110", "ABHD12B", "four messages", "65 percent", "60-70" across slides and notes; no table cell with n < 10 carries a headline.
