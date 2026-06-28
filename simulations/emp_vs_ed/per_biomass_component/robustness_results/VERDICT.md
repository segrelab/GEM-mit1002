# Is the per-biomass-component ED signal real? — No.

## Test
For the amino acids that showed ED flux when their production was maximized
(cysteine, leucine, glutamate, glutamine, alanine, valine, proline), we asked
whether ED is *required* for maximal production or just one of many equally
optimal solutions that pFBA happens to select. Method: FVA on the ED reaction
(`rxn01477`) at 99.9% of maximal production, plus loopless FVA. See
`robustness_check.py`.

## Result: ED is an alternative optimum, not a requirement
| amino acid | ED (pFBA) | ED min | ED max | ED min (loopless) | required? |
|---|---|---|---|---|---|
| Cysteine | 0.0 | 0.0 | 378 | 0.0 | No |
| Leucine | 10.0 | 0.0 | 401 | 0.0 | No |
| Glutamate | 10.0 | 0.0 | 402 | 0.0 | No |
| Glutamine | 10.0 | 0.0 | 397 | 0.0 | No |
| Alanine | 10.0 | 0.0 | 384 | 0.0 | No |
| Valine | 10.0 | 0.0 | 385 | 0.0 | No |
| Proline | 10.0 | 0.0 | 388 | 0.0 | No |

For every amino acid the ED flux can be **0** while production stays at 99.9% of
its maximum — so ED is never required. pFBA reports a non-zero ED flux (10) only
because, among many equally parsimonious solutions, the solver's tie-break lands
on ED. The result holds under loopless FVA, so it is not a loop artifact; it is a
genuine alternative-optimum (degeneracy) artifact. EMP is equally unconstrained.

## Implication for the paper
The per-biomass-component maximization (and the ±1-log stoichiometric
perturbation, which shares the same degeneracy) do not provide evidence that
biomass composition drives ED use. They should not be included as ED-rationale
results. This is consistent with the NADPH titration finding that ED is never
selected as a redox source.

## Robustly supported story
ED is growth-suboptimal on glucose (forced_ed_tradeoffs) and confers no
stoichiometric advantage for amino-acid production or NADPH supply. The one
robust driver of ED use is substrate access: uronic acids that enter at KDPG
obligately require it (tradeoff_on_different_substrates).
