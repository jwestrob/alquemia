# Why 1KB0 still needs a new full-system preparation

The pinned raw mmCIF resolves the first rejected connection: TRO512 is present
and covalently joins LEU511 to GLU513. The old fixed-core-oriented standard-residue
normalization removed it; subsequent whole-chain parsing joined511directlyto513.
This is omitted modified chemistry, not a missing crystallographic residue.
A separate five-residue loop gap573to579 is genuinely absent in the archived
chain and spans18.969263Angstrom in the normalized preparation.

The raw structure also contains heme c802, covalently joined to CYS604/CYS607;
its iron contacts HIS608/MET647. The original normalization deliberately excluded
noncore heterogens. That is a recorded fixed-core model choice; it does not
supply a chemically complete full-system input. The deposited structure and
primary structural paper identify its cytochrome domain and heme chemistry.
[Primary structure and citation](https://www.rcsb.org/structure/1KB0).

A simple terminal-cap repair cannot resolve all these problems. Do not replace
TRO with TRP, omit heme, guess iron spin/charge, bridge the gap, choose another
site, or count a changed preparation as the original transfer passing. Keep
1KB0 unsupported in the completed canonical panel. A future explicit full-system
preparation needs all three issues resolved under a new recorded policy.

Local source pins and exact deposited covalent connections are retained in
ONE_KB0_RAW_CHEMISTRY.json and workspaces/mace_omol_20260917/
one_kb0_raw_chemistry_v1/result.json. No preparation, score, new inference or
baseline change was performed in this read-only diagnosis.
