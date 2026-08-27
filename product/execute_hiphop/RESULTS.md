# Product steps 2+3 results

Executed 2026-08-27 under product/transfer_route_b/DESIGN.md and
product/execute_hiphop/DESIGN.md (configs: product_transfer.json,
product_execute.json).

## Step 2 (transfer, route B verbatim)

1,177/1,177 universe proteins through the trained injection projection
-> 1,177 yeast task rankings (results/, 0 embedding NaN).

## Step 3 (HIP/HOP execution)

- 3,825,250 (target x InChIKey) best-dose pairs evaluated.
- **7,100 pairs q<0.1 over 1,038 targets and 484 compounds**; rho
  0.014-0.273.
- Top targets by max rho (with scF direction score): STING1 0.273 (0.41),
  GPR149 0.254 (0.37), ADGRG4 0.248 (0.38), GPR75 0.246 (0.38), GPR160
  0.244 (0.39), ADGRV1 0.239 (0.39), RYR3 0.230 (0.41), RYR2 0.229
  (0.40), KCNH7 0.213 (0.38), PIEZO2 0.187 (0.35), CATSPER3 0.186 (0.37).
- Compound breadth: median 4 significant targets per compound; tail up
  to 317 (see DESIGN's registered observation - specificity treatment
  required in step 4).

## Next

Step 4 (convergence): B-reversed back-mapping precision check, compound
known-target intersection, specificity treatment, final target selection.
