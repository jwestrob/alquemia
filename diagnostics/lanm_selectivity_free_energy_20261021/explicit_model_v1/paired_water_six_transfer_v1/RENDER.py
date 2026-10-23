"""Render only declared water-pair changes on the existing consumed32 inputs."""
import math


def require(condition, message):
    if not condition:
        raise ValueError(message)


def render(adapter, cfg, row, values):
    """The caller supplies the pinned existing transfer DATA adapter."""
    require(len(values) == 6 and all(
        isinstance(v, (int, float)) and not isinstance(v, bool)
        and math.isfinite(v) for v in values), "Six finite parameters required")
    for j, value in enumerate(values[:4]):
        require((.05 <= value <= .5) if j % 2 else
                (.85 * adapter.e.INITIAL[j] <= value <= 1.15 * adapter.e.INITIAL[j]),
                "Original radius/damping bounds")
    for j, metal in enumerate(("La", "Dy")):
        require(.5 * adapter.e.EPSILON[metal] <= values[4 + j]
                <= 1.5 * adapter.e.EPSILON[metal], "Declared epsilon bounds")
    resolver = adapter.module(cfg["resolver"], "six_transfer_water_resolver")
    prm = adapter.verify(row["input"]["parameters"]).read_text()
    original_key = adapter.verify(row["input"]["key"]).read_text()
    before = resolver.parse_parameters(prm, original_key)
    metal = row["geometry"]["metal"]
    require(metal in ("La", "Dy"), "Metal outside declared pair")
    mt, offset, epsilon_index = (500, 0, 4) if metal == "La" else (501, 2, 5)
    epsilon = values[epsilon_index]
    require([t for t, atom in before["atom"].items() if atom["class"] == 34]
            == [36], "Unexpected water class alias")
    recipient = 36 in row["input"]["types"]
    require(recipient == (row["input"]["id"] not in adapter.DRY),
            "Dry/wet recipient mismatch")
    unchanged = not recipient or (
        list(values[offset:offset + 2]) == list(adapter.e.INITIAL[offset:offset + 2])
        and epsilon == adapter.e.EPSILON[metal])
    key = original_key if unchanged else original_key.rstrip() + (
        f"\n\n# Paired water development override: {metal} only\n"
        f"vdwpr 34 {mt} {values[offset]:.17g} {epsilon:.17g}\n"
        f"polpair 36 {mt} {values[offset + 1]:.17g}\n")
    after = resolver.parse_parameters(prm, key)
    for field in ("settings", "atom", "vdw", "polarize", "metal_multipoles"):
        require(before[field] == after[field], "Unrelated parameters changed")
    for field, allowed in (("vdwpr", (34, mt)), ("polpair", (36, mt))):
        require({k: v for k, v in before[field].items() if k != allowed}
                == {k: v for k, v in after[field].items() if k != allowed},
                "Non-water override changed")
    for atom_type in set(row["input"]["types"]) - {mt, 36} | {37}:
        require(resolver.resolve_pair(before, mt, atom_type)
                == resolver.resolve_pair(after, mt, atom_type),
                "Water-H or donor pair changed")
    pair = resolver.resolve_pair(after, mt, 36)
    if recipient:
        require(pair["Rmin_pair_A"] == values[offset]
                and pair["thole"] == values[offset + 1]
                and pair["epsilon_kcal_mol"] == epsilon, "Rendered pair differs")
    return key, {"water_recipient": recipient, "original_key_reused": unchanged,
                 "resolved_water_pair": pair,
                 "effective_values": list(values[offset:offset + 2]) + [epsilon]
                 if recipient else None}
