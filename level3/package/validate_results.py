"""Format check only. Compares nothing with any reference. Usage: python3 validate_results.py results.json"""
import sys, json
r = json.load(open(sys.argv[1])); cases = json.load(open("cases.json"))["cases"]; ok = True
def need(c, msg):
    global ok
    if not c: ok = False; print("FORMAT PROBLEM:", msg)
need(len(r["pairs"]["disc"]) == 4000 and len(r["pairs"]["cyl"]) == 4000, "pairs must have 4000 entries each")
for k in ("m0_t0.0", "m1_t0.25", "m1_t0.05", "sq_L8", "tri_L8", "sph_L8"): need(len(r["configs"].get(k, [])) == 20, f"configs[{k}] must have 20 entries")
for k in ("square", "triangle", "spheroid"): need(len(r["pairs"].get(k, [])) == 4000, f"pairs[{k}] must have 4000 entries")
for c, spec in cases.items():
    s = r["stats"].get(c); need(s is not None, f"stats missing {c}")
    if s: need(all(k in s for k in ("mean", "se", "nreal")) and s["nreal"] > 0 and s["mean"] > 0, f"stats[{c}] incomplete")
    if spec.get("record_raw"):
        raw = r["raw"].get(c); need(raw is not None and all(k in raw and len(raw[k]) == (s or {}).get("nreal", -1) for k in ("N_any", "N_x", "N_y", "N_z")), f"raw[{c}] arrays must have nreal entries")
        if raw and all(len(raw[k]) == len(raw["N_any"]) for k in raw):
            bad = sum(1 for i in range(len(raw["N_any"])) if raw["N_any"][i] != min(raw["N_x"][i], raw["N_y"][i], raw["N_z"][i]))
            need(bad == 0, f"raw[{c}]: N_any must equal min(N_x, N_y, N_z); {bad} realisations violate it")
print("FORMAT OK" if ok else "FORMAT ERRORS FOUND")
