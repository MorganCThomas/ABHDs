# ABHD2 vs MGLL

- Structure: TM-score 0.63 (normalised by ABHD2) / 0.62 (by MGLL); RMSD 3.6 Å over 225 aligned residues (294 and 303 in the domains)
- Sequence identity (global alignment of the trimmed domains): 21%
- Same paralog group or FunFam: no
- Catalytic categories: ABHD2 Ser-His-Asp; MGLL Ser-His-Asp
- EC3: ABHD2 3.1.1; MGLL 3.1.1
- Second oxyanion donor: ABHD2 I134; MGLL A51
- Active-site cavity (median grid units, 3.2 Å, like-for-like): ABHD2 283.0; MGLL 264.0 (model-dependent; see report 12)

## Distance in each view

| view | distance | percentile_of_all_pairs | rank_of_b_for_a | rank_of_a_for_b |
|---|---|---|---|---|
| structure | 0.374 | 24.454 | 23 | 56 |
| structure_short | 0.366 | 50.142 | 45 | 71 |
| structure_core | 0.310 | 11.860 | 11 | 31 |
| sequence | 0.972 | 26.422 | 18 | 52 |
| funfam | 296.447 | 31.363 | 22 | 60 |
| elbow_window | 0.762 | 5.363 | 4 | 11 |
| catalytic_windows | 0.902 | 25.456 | 37 | 62 |
| local_site | 0.004 | 23.987 | 59 | 25 |
| triad | 0.333 | 27.707 | 38 | 49 |
| language_model | 0.069 | 9.741 | 25 | 22 |

## Catalytic residues (aligned by TM-align)

| role | a_residue | b_aligned_residue | b_own_residue | ca_displacement |
|---|---|---|---|---|
| nucleophile position | S207 | S122 | S122 | 1.1 Å |
| His | H376 | H269 | H269 | 1.0 Å |
| acid | D345 | D239 | D239 | 1.6 Å |

## Pocket-lining residues of ABHD2 mapped onto MGLL

4 of 24 identical, 12 of the same physicochemical class (random core positions of this protein pair would be identical in roughly 15–25%).

| a | b_aligned | identical | b_also_lines_pocket |
|---|---|---|---|
| I134 | A51 |  | yes |
| F177 | S91 |  |  |
| T178 | D92 |  |  |
| Y179 | F93 |  |  |
| S207 | S122 | yes | yes |
| L208 | M123 |  | yes |
| C234 | S146 |  |  |
| Q235 | P147 |  |  |
| G236 | L148 |  | yes |
| A242 | A151 | yes | yes |
| Q243 | N152 |  | yes |
| F246 | S155 |  | yes |
| M247 | A156 |  | yes |
| N257 | F159 |  | yes |
| L300 | P178 |  | yes |
| I303 | L205 |  | yes |
| D304 | K206 |  |  |
| M308 | G210 |  | yes |
| Y320 | N215 |  |  |
| Y321 | A216 |  |  |
| L347 | L241 | yes | yes |
| V348 | C242 |  | yes |
| L352 | G246 |  |  |
| H376 | H269 | yes |  |
