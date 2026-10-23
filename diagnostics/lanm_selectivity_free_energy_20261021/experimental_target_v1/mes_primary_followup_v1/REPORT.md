# MES primary-source follow-up

## Result

**No numerical La–MES constant or matched La/Dy MES correction is newly admitted.** The existing observable contract remains unchanged. This bounded follow-up recovered a useful primary methodological warning and a lower-confidence conflicting observation, but neither supplies the missing assay conversion.

No molecular calls, submissions, fitting, paid purchase or contact occurred. Reserved validation outcomes were not opened. Search results mentioning unrelated biological systems were not followed or used as new validation data.

## Evidence recovered

| Source | Direct evidence | Admission |
|---|---|---|
| [Anwar and Azab 2001](https://pubs.acs.org/doi/10.1021/je000187%2B), JCED 46, 613–618 | Publisher abstract confirms La/Ce/Pr/Eu,25.0±0.1°C and0.1M KNO3; binary hydroxo/dimeric species are included in the study | Abstract-level identity and medium only. La numerical values, MES protonation and equilibrium definitions remain unavailable |
| [Azab et al.2011 primary PDF](https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Chemical%20%26%20Engineering%20Data/56/je200099n.pdf), Table 13/p2620 | Recovered indexed primary table again: 1:1 Dy–MES logK4.19±0.02; Eu4.27±0.03;25.0±0.1°C,0.1M KNO3. Its uncertainty footnote labels3s | Confirms the existing table record, not a new correction. Complete free-ligand/protonation definitions were not recovered; no La row appears in this table |
| [Beaubien 2015 thesis](https://espace.inrs.ca/id/eprint/3351/1/T00812.pdf), printed p23/PDF p41 | Reports that ion-exchange experiments by Sébastien Leguay found no more than10% La/Ce–MES association and could not determine reliable constants | This is a thesis report of another investigator's tests, not recovered underlying primary data. No matching test conditions, raw results or Dy measurement are supplied in that passage; **not an empirical assay bound** |
| [Leguay, Campbell and Fortin 2016](https://doi.org/10.1071/EN15136), [institutional author manuscript](https://espace.inrs.ca/4699/1/P2923_PP.pdf), methods | The primary experimental study explicitly omits pH buffers to avoid rare-earth buffer complexation, citing Anwar 2001 | A methodological concern, not a MES constant or a verification of the thesis's weak-binding claim. This manuscript contains no MES numerical determination |

The nominal temperature and nitrate medium of the 2001/2011 potentiometric papers match one another, but this alone does not establish compatible equilibrium conventions. They differ from the10°C,50mM MES/100mM NaCl 2020 peptide assay. The 2025 assay temperature remains separately unspecified.

## Access and stopping record

After reading the earlier search inventory, this follow-up used the specific DOI, publisher, DOI registry/repository locations, the indexed primary 2011 table and one new laboratory-method lead. The direct 2001 publisher PDF returned HTTP403; its accessible publisher page is an abstract/purchase page. Crossref supplies only that same PDF link. OpenAlex currently reports closed access, no repository full text and no alternative location. The indexed 2011 PDF could be searched, but direct download failed with no route to host; the guessed 2001 mirror URL was also inaccessible. No access restriction was bypassed, and inaccessible requests were not repeatedly retried.

The thesis and 2016 manuscript were retrieved from INRS and retained with hashes. The 2016 follow-up closes the thesis lead for this bounded search: its buffer-free methods cannot provide a numerical La–MES value. Further broad searches or reusing a review's La value would not resolve the primary-data gap.

## Consequence for the observable

Keep `alpha_M = total unbound soluble metal / free M3+` and the conditional exchange correction `RT ln(alpha_Dy/alpha_La)` separate from the intrinsic exchange calculation. The old excluded secondary La value remains excluded. Neither a numeric correction nor a finite empirical interval is justified by this follow-up; do not assume MES is inert, and do not transfer the thesis's≤10% observation to 50 mM MES/pH6.

To admit a future correction, the required evidence is the actual 2001 La binary/protonation tables and definitions, compatible Dy definitions including relevant stoichiometries, and a defensible transfer between temperature/media—or a paired free-metal/buffer-concentration measurement under the chosen assay. This is an evidence requirement, not a request to purchase/contact anyone or perform new experiments.

Machine-readable provenance, retrieval outcomes and source hashes are in `EVIDENCE.json` and the matching workspace. Original contract, mapping inputs and results are preserved unchanged.
