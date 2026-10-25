"""Deterministic inputs (step 3) and the cases in cases.json (step 4).

Run from anywhere; the script changes to the workspace root.
"""
from __future__ import annotations

import json
import math
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "work"))

import core

MASTER_SEED = 20261025
CHECKPOINT = ROOT / "work" / "checkpoints"
PROGRESS = ROOT / "work" / "PROGRESS.txt"


def log(msg):
    line = msg if msg.endswith("\n") else msg + "\n"
    print(line, end="", flush=True)
    PROGRESS.parent.mkdir(parents=True, exist_ok=True)
    with PROGRESS.open("a", encoding="utf-8") as f:
        f.write(line)


def parse_n_insert(expr, L):
    """Upper limit written as 'factor*L^3'. Use the greatest integer not above it."""
    factor = float(str(expr).split("*")[0])
    return int(math.floor(factor * (L ** 3) + 1e-9))


def pick_ncell(L, typical_range):
    """Cell edge a little above the typical interaction range, at least ~4 cells when L allows."""
    target = max(float(typical_range) * 1.05, 1.0)
    cap = max(L / 4.0, 1.0)
    target = min(target, cap)
    ncell = int(math.floor(L / target))
    if ncell < 1:
        ncell = 1
    if ncell > 48:
        ncell = 48
    return ncell


def rng_for(case_index, real_index):
    ss = np.random.SeedSequence([MASTER_SEED, int(case_index), int(real_index)])
    return np.random.Generator(np.random.PCG64(ss))


def _normals_from_zphi(z, phi):
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=-1)


def _frames_from_normals(nrm, spin):
    """Vectorised frame_from_normal. spin shape (n,)."""
    n = nrm.shape[0]
    e1 = np.zeros((n, 3))
    e2 = np.zeros((n, 3))
    ref = np.zeros((n, 3))
    use_z = np.abs(nrm[:, 2]) < 0.9
    ref[use_z, 2] = 1.0
    ref[~use_z, 0] = 1.0
    a = np.cross(ref, nrm)
    a /= np.linalg.norm(a, axis=1)[:, None]
    b = np.cross(nrm, a)
    c = np.cos(spin)[:, None]
    s = np.sin(spin)[:, None]
    e1r = c * a + s * b
    e2r = np.cross(nrm, e1r)
    return e1r, e2r


def _fisher_about(directors, u, phi, kappa):
    ek = math.exp(-2.0 * kappa)
    arg = u + (1.0 - u) * ek
    arg = np.maximum(arg, 1e-300)
    c = 1.0 + np.log(arg) / kappa
    c = np.clip(c, -1.0, 1.0)
    s = np.sqrt(np.maximum(0.0, 1.0 - c * c))
    # Same frame convention as core.frame_from_normal(..., spin=0).
    ref = np.zeros_like(directors)
    use_z = np.abs(directors[:, 2]) < 0.9
    ref[use_z, 2] = 1.0
    ref[~use_z, 0] = 1.0
    e1 = np.cross(ref, directors)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(directors, e1)
    return (s * np.cos(phi))[:, None] * e1 + (s * np.sin(phi))[:, None] * e2 + c[:, None] * directors


def _nearest_directors(pos, seeds, directors, L):
    raw = seeds[None, :, :] - pos[:, None, :]
    shift = -np.floor(raw / L + 0.5)
    d = raw + shift * L
    dist2 = np.einsum("ijk,ijk->ij", d, d)
    k = np.argmin(dist2, axis=1)
    return directors[k]


def generate(kind, L, n, rng, spec):
    """Draw one realisation. Random-number order is documented in NOTES.md."""
    pos = np.empty((n, 3))
    nrm = np.empty((n, 3))
    e1 = np.empty((n, 3))
    e2 = np.empty((n, 3))
    rad = np.empty(n)
    aux = np.empty(n)
    rb = np.empty(n)

    if kind == "thin_disc":
        u = rng.random((n, 5))
        pos[:] = u[:, 0:3] * L
        nrm[:] = _normals_from_zphi(u[:, 3] * 2 - 1, u[:, 4] * 2 * math.pi)
        e1[:] = 0.0
        e2[:] = 0.0
        rad[:] = 1.0
        aux[:] = 0.0
        rb[:] = 1.0
        sim_kind = 0
        typical = 2.0
    elif kind in ("cylinder", "cylinder_aligned", "cylinder_lognormal"):
        R = float(spec.get("R", 1.0))
        t = float(spec["t"])
        h = 0.5 * t
        if kind == "cylinder":
            u = rng.random((n, 5))
            pos[:] = u[:, 0:3] * L
            nrm[:] = _normals_from_zphi(u[:, 3] * 2 - 1, u[:, 4] * 2 * math.pi)
            rad[:] = R
            aux[:] = h
            rb[:] = math.sqrt(R * R + h * h)
        elif kind == "cylinder_lognormal":
            mu = float(spec["mean_lnR"])
            sigma = float(spec["sigma_lnR"])
            cap = float(spec["R_cap"])
            for i in range(n):
                pos[i] = rng.random(3) * L
                while True:
                    Ri = math.exp(float(rng.normal(mu, sigma)))
                    if Ri <= cap:
                        break
                rad[i] = Ri
                z = float(rng.random()) * 2.0 - 1.0
                phi = float(rng.random()) * 2.0 * math.pi
                rr = math.sqrt(max(0.0, 1.0 - z * z))
                nrm[i, 0] = rr * math.cos(phi)
                nrm[i, 1] = rr * math.sin(phi)
                nrm[i, 2] = z
                aux[i] = h
                rb[i] = math.sqrt(Ri * Ri + h * h)
        else:
            n_seeds = int(round((L / float(spec["domain_size"])) ** 3))
            su = rng.random((n_seeds, 5))
            seeds = su[:, 0:3] * L
            directors = _normals_from_zphi(su[:, 3] * 2 - 1, su[:, 4] * 2 * math.pi)
            u = rng.random((n, 5))
            pos[:] = u[:, 0:3] * L
            assigned = _nearest_directors(pos, seeds, directors, L)
            nrm[:] = _fisher_about(assigned, u[:, 3], u[:, 4] * 2 * math.pi, float(spec["kappa"]))
            rad[:] = R
            aux[:] = h
            rb[:] = math.sqrt(R * R + h * h)
        e1[:] = 0.0
        e2[:] = 0.0
        sim_kind = 1
        typical = 2.0 * float(np.median(rb))
    elif kind in ("square_plate", "triangle_plate"):
        u = rng.random((n, 6))
        pos[:] = u[:, 0:3] * L
        nrm[:] = _normals_from_zphi(u[:, 3] * 2 - 1, u[:, 4] * 2 * math.pi)
        e1[:], e2[:] = _frames_from_normals(nrm, u[:, 5] * 2 * math.pi)
        rad[:] = 1.0
        aux[:] = 0.0
        if kind == "square_plate":
            rb[:] = core.RB_SQUARE
            sim_kind = 2
            typical = 2.0 * core.RB_SQUARE
        else:
            rb[:] = core.RB_TRIANGLE
            sim_kind = 3
            typical = 2.0 * core.RB_TRIANGLE
    elif kind == "oblate_spheroid":
        R = float(spec["R"])
        c = float(spec["c"])
        u = rng.random((n, 5))
        pos[:] = u[:, 0:3] * L
        nrm[:] = _normals_from_zphi(u[:, 3] * 2 - 1, u[:, 4] * 2 * math.pi)
        # Rotationally symmetric: in-plane spin fixed at 0.
        e1[:], e2[:] = _frames_from_normals(nrm, np.zeros(n))
        rad[:] = R
        aux[:] = c
        rb[:] = max(R, c)
        sim_kind = 4
        typical = 2.0 * max(R, c)
    else:
        raise ValueError(kind)
    return sim_kind, pos, nrm, e1, e2, rad, aux, rb, typical


def observable(kind, N, L, rad, spec):
    if N < 0:
        return None
    if kind in ("thin_disc", "square_plate", "triangle_plate", "oblate_spheroid"):
        return (4.0 * math.pi / 3.0) * N / (L ** 3)
    t = float(spec["t"])
    if kind == "cylinder_lognormal":
        return math.pi * t * float(np.sum(rad[:N] ** 2)) / (L ** 3)
    R = float(spec.get("R", 1.0))
    return N * (math.pi * R * R * t) / (L ** 3)


def _simulate_arrays(sim_kind, pos, nrm, e1, e2, rad, aux, rb, L, typical, stop_all):
    ncell = pick_ncell(L, typical)
    return core.simulate(sim_kind, pos, nrm, e1, e2, rad, aux, rb, L, ncell, 1 if stop_all else 0), ncell


def run_realisation(kind, spec, L, n_cap, case_index, real_index, stop_all):
    """Return N tuple, observable, flags. May extend raw realisations that have wrapped."""
    caps = [n_cap]
    if stop_all:
        caps = [n_cap, n_cap * 2, n_cap * 3]
    extended = 0
    censored = 0
    out = None
    rad = None
    used = n_cap
    for attempt, n_use in enumerate(caps):
        rng = rng_for(case_index, real_index)
        sim_kind, pos, nrm, e1, e2, rad, aux, rb, typical = generate(kind, L, n_use, rng, spec)
        out, _ = _simulate_arrays(sim_kind, pos, nrm, e1, e2, rad, aux, rb, L, typical, stop_all)
        used = n_use
        if attempt > 0:
            extended = 1
        Nany, Nx, Ny, Nz = (int(out[0]), int(out[1]), int(out[2]), int(out[3]))
        if not stop_all:
            break
        if Nany < 0:
            break
        if Nx > 0 and Ny > 0 and Nz > 0:
            break
    Nany, Nx, Ny, Nz = (int(out[0]), int(out[1]), int(out[2]), int(out[3]))
    if stop_all and Nany > 0 and (Nx < 0 or Ny < 0 or Nz < 0):
        # Censored directional time: the search stopped at `used` insertions.
        # The stored value is that cap so min(N_x,N_y,N_z) still equals N_any.
        censored = 1
        if Nx < 0:
            Nx = used
        if Ny < 0:
            Ny = used
        if Nz < 0:
            Nz = used
    if Nany > 0 and Nx > 0 and Ny > 0 and Nz > 0:
        if Nany != min(Nx, Ny, Nz):
            raise RuntimeError(f"N_any mismatch {Nany, Nx, Ny, Nz}")
    value = observable(kind, Nany, L, rad, spec)
    return {
        "N_any": Nany, "N_x": Nx, "N_y": Ny, "N_z": Nz,
        "value": value, "extended": extended, "censored": censored, "used": used,
    }


# ---------------------------------------------------------------------------
# Deterministic inputs
# ---------------------------------------------------------------------------

def _in_cyl(x, p, n, R, h, tol=1e-5):
    d = x - p
    axial = float(np.dot(d, n))
    if abs(axial) > h + tol:
        return False
    return float(np.dot(d, d) - axial * axial) <= (R + tol) ** 2


def _ext_cyl(n, R, h, u):
    nu = abs(float(np.dot(n, u)))
    return h * nu + R * math.sqrt(max(0.0, 1.0 - nu * nu))


def _in_ell(x, p, e1, e2, n, R, c, tol=1e-5):
    d = x - p
    a1 = float(np.dot(d, e1))
    a2 = float(np.dot(d, e2))
    an = float(np.dot(d, n))
    return (a1 / (R + tol)) ** 2 + (a2 / (R + tol)) ** 2 + (an / (c + tol)) ** 2 <= 1.0


def _ext_ell(e1, e2, n, R, c, u):
    return math.sqrt((R * float(np.dot(e1, u))) ** 2 + (R * float(np.dot(e2, u))) ** 2
                     + (c * float(np.dot(n, u))) ** 2)


def certify_pairs(pairs):
    """Check witnesses and separating directions on the deterministic pairs. Returns a dict of counts."""
    z = np.load(ROOT / "inputs" / "pairs_inputs.npz")
    n1, n2, d = z["n1"], z["n2"], z["d"]
    R1, R2, h1, h2 = z["R1"], z["R2"], z["h1"], z["h2"]
    disc_bad = cyl_hit_bad = cyl_sep_bad = 0
    for i in range(n1.shape[0]):
        hit, wx, wy, wz = core.disc_analyze(
            0, 0, 0, n1[i, 0], n1[i, 1], n1[i, 2], R1[i],
            d[i, 0], d[i, 1], d[i, 2], n2[i, 0], n2[i, 1], n2[i, 2], R2[i])
        if bool(hit) != pairs["disc"][i]:
            disc_bad += 1
        if hit:
            w = np.array([wx, wy, wz])
            # Plane and radius membership.
            if abs(float(np.dot(w, n1[i]))) > 1e-6 or float(np.dot(w, w)) > (R1[i] + 1e-6) ** 2:
                disc_bad += 1
            r2 = w - d[i]
            if abs(float(np.dot(r2, n2[i]))) > 1e-6 or float(np.dot(r2, r2)) > (R2[i] + 1e-6) ** 2:
                disc_bad += 1
        ch, wx, wy, wz, vx, vy, vz = core.cylinder_analyze(
            0, 0, 0, n1[i, 0], n1[i, 1], n1[i, 2], R1[i], h1[i],
            d[i, 0], d[i, 1], d[i, 2], n2[i, 0], n2[i, 1], n2[i, 2], R2[i], h2[i])
        if bool(ch) != pairs["cyl"][i]:
            cyl_hit_bad += 1
        if ch:
            w = np.array([wx, wy, wz])
            if not _in_cyl(w, np.zeros(3), n1[i], R1[i], h1[i]) or not _in_cyl(w, d[i], n2[i], R2[i], h2[i]):
                cyl_hit_bad += 1
        else:
            v = np.array([vx, vy, vz], dtype=np.float64)
            dist = float(np.linalg.norm(v))
            rb = math.sqrt(R1[i] ** 2 + h1[i] ** 2) + math.sqrt(R2[i] ** 2 + h2[i] ** 2)
            if dist >= 1e-8:
                u = v / dist
                gap = float(np.dot(-d[i], u)) - _ext_cyl(n1[i], R1[i], h1[i], u) - _ext_cyl(n2[i], R2[i], h2[i], u)
                if abs(gap - dist) > 1e-4 or gap < -1e-5:
                    cyl_sep_bad += 1
            elif float(np.linalg.norm(d[i])) <= rb + 1e-6:
                # Miss with a zero separator that is not a bounding-sphere reject.
                cyl_sep_bad += 1
    # Shapes
    zs = np.load(ROOT / "inputs" / "shapes_pairs_inputs.npz")
    plate_bad = {"square": 0, "triangle": 0}
    for name, nv in (("square", 4), ("triangle", 3)):
        for i in range(4000):
            hit, wx, wy, wz = core.plate_analyze(
                nv, 0, 0, 0,
                zs[f"{name}_e1_1"][i, 0], zs[f"{name}_e1_1"][i, 1], zs[f"{name}_e1_1"][i, 2],
                zs[f"{name}_e2_1"][i, 0], zs[f"{name}_e2_1"][i, 1], zs[f"{name}_e2_1"][i, 2],
                zs[f"{name}_n1"][i, 0], zs[f"{name}_n1"][i, 1], zs[f"{name}_n1"][i, 2],
                zs[f"{name}_d"][i, 0], zs[f"{name}_d"][i, 1], zs[f"{name}_d"][i, 2],
                zs[f"{name}_e1_2"][i, 0], zs[f"{name}_e1_2"][i, 1], zs[f"{name}_e1_2"][i, 2],
                zs[f"{name}_e2_2"][i, 0], zs[f"{name}_e2_2"][i, 1], zs[f"{name}_e2_2"][i, 2],
                zs[f"{name}_n2"][i, 0], zs[f"{name}_n2"][i, 1], zs[f"{name}_n2"][i, 2])
            if bool(hit) != pairs[name][i]:
                plate_bad[name] += 1
            if hit:
                w = np.array([wx, wy, wz])
                if not math.isfinite(wx):
                    plate_bad[name] += 1
    sph_hit_bad = sph_sep_bad = 0
    for i in range(4000):
        hit, wx, wy, wz, vx, vy, vz = core.spheroid_analyze(
            0, 0, 0,
            zs["spheroid_e1_1"][i, 0], zs["spheroid_e1_1"][i, 1], zs["spheroid_e1_1"][i, 2],
            zs["spheroid_e2_1"][i, 0], zs["spheroid_e2_1"][i, 1], zs["spheroid_e2_1"][i, 2],
            zs["spheroid_n1"][i, 0], zs["spheroid_n1"][i, 1], zs["spheroid_n1"][i, 2],
            float(zs["spheroid_R1"][i]), float(zs["spheroid_c1"][i]),
            zs["spheroid_d"][i, 0], zs["spheroid_d"][i, 1], zs["spheroid_d"][i, 2],
            zs["spheroid_e1_2"][i, 0], zs["spheroid_e1_2"][i, 1], zs["spheroid_e1_2"][i, 2],
            zs["spheroid_e2_2"][i, 0], zs["spheroid_e2_2"][i, 1], zs["spheroid_e2_2"][i, 2],
            zs["spheroid_n2"][i, 0], zs["spheroid_n2"][i, 1], zs["spheroid_n2"][i, 2],
            float(zs["spheroid_R2"][i]), float(zs["spheroid_c2"][i]))
        if bool(hit) != pairs["spheroid"][i]:
            sph_hit_bad += 1
        if hit:
            w = np.array([wx, wy, wz])
            if not (math.isfinite(wx) and _in_ell(
                    w, np.zeros(3), zs["spheroid_e1_1"][i], zs["spheroid_e2_1"][i], zs["spheroid_n1"][i],
                    float(zs["spheroid_R1"][i]), float(zs["spheroid_c1"][i]))
                    and _in_ell(
                    w, zs["spheroid_d"][i], zs["spheroid_e1_2"][i], zs["spheroid_e2_2"][i], zs["spheroid_n2"][i],
                    float(zs["spheroid_R2"][i]), float(zs["spheroid_c2"][i]))):
                sph_hit_bad += 1
        else:
            v = np.array([vx, vy, vz], dtype=np.float64)
            dist = float(np.linalg.norm(v))
            if dist >= 1e-8:
                u = v / dist
                gap = float(np.dot(-zs["spheroid_d"][i], u))
                gap -= _ext_ell(zs["spheroid_e1_1"][i], zs["spheroid_e2_1"][i], zs["spheroid_n1"][i],
                                float(zs["spheroid_R1"][i]), float(zs["spheroid_c1"][i]), u)
                gap -= _ext_ell(zs["spheroid_e1_2"][i], zs["spheroid_e2_2"][i], zs["spheroid_n2"][i],
                                float(zs["spheroid_R2"][i]), float(zs["spheroid_c2"][i]), u)
                if abs(gap - dist) > 1e-4 or gap < -1e-5:
                    sph_sep_bad += 1
    return {
        "disc_bad": disc_bad,
        "cyl_hit_bad": cyl_hit_bad,
        "cyl_sep_bad": cyl_sep_bad,
        "square_bad": plate_bad["square"],
        "triangle_bad": plate_bad["triangle"],
        "sph_hit_bad": sph_hit_bad,
        "sph_sep_bad": sph_sep_bad,
    }


def compute_pairs():
    z = np.load(ROOT / "inputs" / "pairs_inputs.npz")
    n1, n2, d = z["n1"], z["n2"], z["d"]
    R1, R2, h1, h2 = z["R1"], z["R2"], z["h1"], z["h2"]
    disc = []
    cyl = []
    for i in range(n1.shape[0]):
        disc.append(bool(core.disc_hit(
            0, 0, 0, n1[i, 0], n1[i, 1], n1[i, 2], float(R1[i]),
            d[i, 0], d[i, 1], d[i, 2], n2[i, 0], n2[i, 1], n2[i, 2], float(R2[i]))))
        cyl.append(bool(core.cylinder_hit(
            0, 0, 0, n1[i, 0], n1[i, 1], n1[i, 2], float(R1[i]), float(h1[i]),
            d[i, 0], d[i, 1], d[i, 2], n2[i, 0], n2[i, 1], n2[i, 2], float(R2[i]), float(h2[i]))))
    zs = np.load(ROOT / "inputs" / "shapes_pairs_inputs.npz")
    out = {"disc": disc, "cyl": cyl}
    for name, nv in (("square", 4), ("triangle", 3)):
        flags = []
        e11, e21, n1s = zs[f"{name}_e1_1"], zs[f"{name}_e2_1"], zs[f"{name}_n1"]
        e12, e22, n2s = zs[f"{name}_e1_2"], zs[f"{name}_e2_2"], zs[f"{name}_n2"]
        disp = zs[f"{name}_d"]
        for i in range(disp.shape[0]):
            flags.append(bool(core.plate_hit(
                nv, 0, 0, 0, e11[i, 0], e11[i, 1], e11[i, 2], e21[i, 0], e21[i, 1], e21[i, 2],
                n1s[i, 0], n1s[i, 1], n1s[i, 2],
                disp[i, 0], disp[i, 1], disp[i, 2],
                e12[i, 0], e12[i, 1], e12[i, 2], e22[i, 0], e22[i, 1], e22[i, 2],
                n2s[i, 0], n2s[i, 1], n2s[i, 2])))
        out[name] = flags
    sph = []
    for i in range(zs["spheroid_d"].shape[0]):
        sph.append(bool(core.spheroid_hit(
            0, 0, 0,
            zs["spheroid_e1_1"][i, 0], zs["spheroid_e1_1"][i, 1], zs["spheroid_e1_1"][i, 2],
            zs["spheroid_e2_1"][i, 0], zs["spheroid_e2_1"][i, 1], zs["spheroid_e2_1"][i, 2],
            zs["spheroid_n1"][i, 0], zs["spheroid_n1"][i, 1], zs["spheroid_n1"][i, 2],
            float(zs["spheroid_R1"][i]), float(zs["spheroid_c1"][i]),
            zs["spheroid_d"][i, 0], zs["spheroid_d"][i, 1], zs["spheroid_d"][i, 2],
            zs["spheroid_e1_2"][i, 0], zs["spheroid_e1_2"][i, 1], zs["spheroid_e1_2"][i, 2],
            zs["spheroid_e2_2"][i, 0], zs["spheroid_e2_2"][i, 1], zs["spheroid_e2_2"][i, 2],
            zs["spheroid_n2"][i, 0], zs["spheroid_n2"][i, 1], zs["spheroid_n2"][i, 2],
            float(zs["spheroid_R2"][i]), float(zs["spheroid_c2"][i]))))
    out["spheroid"] = sph
    return out


def _config_N(kind_sim, pos, nrm, e1, e2, rad, aux, rb, L):
    pos = np.array(pos, dtype=np.float64, copy=True)
    # Keep centres inside the cube. Inputs inspected earlier already lie in (0, L).
    pos %= L
    n = pos.shape[0]
    if e1 is None:
        e1 = np.zeros_like(pos)
        e2 = np.zeros_like(pos)
    out, ncell = _simulate_arrays(kind_sim, pos, nrm, e1, e2, rad, aux, rb, L, float(np.max(rb) * 2), True)
    # Cross-check the first call's N_any against an all-pairs grid (ncell=1) once per batch outside.
    return int(out[0]), int(out[1]), int(out[2]), int(out[3]), ncell


def compute_configs():
    z = np.load(ROOT / "inputs" / "configs_inputs.npz")
    configs = {}
    checks = {}
    specs = {
        "m0_t0.0": (0, 0.0),
        "m1_t0.25": (1, 0.25),
        "m1_t0.05": (1, 0.05),
    }
    L = 6.0
    for name, (sk, t) in specs.items():
        vals = []
        for s in range(20):
            pos = np.array(z[f"{name}_s{s}_pos"], dtype=np.float64)
            nrm = np.array(z[f"{name}_s{s}_nrm"], dtype=np.float64)
            n = pos.shape[0]
            rad = np.ones(n)
            h = 0.5 * t
            aux = np.full(n, h)
            rb = np.full(n, math.sqrt(1.0 + h * h))
            Nany, Nx, Ny, Nz, ncell = _config_N(sk, pos, nrm, None, None, rad, aux, rb, L)
            vals.append(Nany)
            if s == 0:
                # All-pairs (one cell) must agree.
                pos2 = np.array(pos, copy=True) % L
                e1 = np.zeros_like(pos2)
                e2 = np.zeros_like(pos2)
                alt = core.simulate(sk, pos2, nrm, e1, e2, rad, aux, rb, L, 1, 1)
                checks[name] = {"ncell": ncell, "match": bool(int(alt[0]) == Nany and int(alt[1]) == Nx
                                                              and int(alt[2]) == Ny and int(alt[3]) == Nz)}
        configs[name] = vals
    zs = np.load(ROOT / "inputs" / "shapes_configs_inputs.npz")
    L8 = 8.0
    shape_specs = {
        "sq_L8": (2, core.RB_SQUARE, None),
        "tri_L8": (3, core.RB_TRIANGLE, None),
        "sph_L8": (4, 1.0, 0.001),
    }
    for name, (sk, rb_or_R, cval) in shape_specs.items():
        vals = []
        for s in range(20):
            pos = np.array(zs[f"{name}_s{s}_pos"], dtype=np.float64)
            nrm = np.array(zs[f"{name}_s{s}_n"], dtype=np.float64)
            e1 = np.array(zs[f"{name}_s{s}_e1"], dtype=np.float64)
            e2 = np.array(zs[f"{name}_s{s}_e2"], dtype=np.float64)
            n = pos.shape[0]
            if sk == 4:
                rad = np.ones(n)
                aux = np.full(n, cval)
                rb = np.ones(n)  # max(R, c) = 1
            else:
                rad = np.ones(n)
                aux = np.zeros(n)
                rb = np.full(n, rb_or_R)
            Nany, Nx, Ny, Nz, ncell = _config_N(sk, pos, nrm, e1, e2, rad, aux, rb, L8)
            vals.append(Nany)
            if s == 0:
                pos2 = np.array(pos, copy=True) % L8
                alt = core.simulate(sk, pos2, nrm, e1, e2, rad, aux, rb, L8, 1, 1)
                checks[name] = {"ncell": ncell, "match": bool(int(alt[0]) == Nany)}
        configs[name] = vals
    return configs, checks


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

def summarise(vals):
    arr = np.array([v for v in vals if v is not None], dtype=np.float64)
    if arr.size == 0:
        return {"mean": 0.0, "se": 0.0, "nreal": 0, "discarded": len(vals)}
    mean = float(arr.mean())
    if arr.size > 1:
        se = float(arr.std(ddof=1) / math.sqrt(arr.size))
    else:
        se = 0.0
    return {"mean": mean, "se": se, "nreal": int(arr.size), "discarded": int(len(vals) - arr.size)}


def run_case(name, spec, case_index, n_workers):
    ck_path = CHECKPOINT / f"{name}.json"
    if ck_path.exists():
        with ck_path.open(encoding="utf-8") as f:
            saved = json.load(f)
        if saved.get("n_done") == spec["nreal"] and saved.get("complete"):
            log(f"resume {name} from checkpoint")
            return saved
    L = float(spec["L"])
    n_cap = parse_n_insert(spec["n_insert"], L)
    stop_all = 1 if spec.get("record_raw") else 0
    nreal = int(spec["nreal"])
    log(f"start {name} L={L} n_insert={n_cap} nreal={nreal} stop_all={stop_all}")
    t0 = time.perf_counter()
    results = [None] * nreal
    done = {"n": 0}
    lock = threading.Lock()

    def work(i):
        rec = run_realisation(spec["kind"], spec, L, n_cap, case_index, i, stop_all)
        with lock:
            done["n"] += 1
            if done["n"] % 10 == 0 or done["n"] == nreal:
                log(f"  {name} {done['n']}/{nreal}  elapsed={time.perf_counter() - t0:.1f}s")
        return i, rec

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        for i, rec in ex.map(work, range(nreal)):
            results[i] = rec
    values = [rec["value"] if rec["N_any"] > 0 else None for rec in results]
    stats = summarise(values)
    raw = None
    if spec.get("record_raw"):
        kept = [rec for rec in results if rec["N_any"] > 0]
        raw = {
            "N_any": [rec["N_any"] for rec in kept],
            "N_x": [rec["N_x"] for rec in kept],
            "N_y": [rec["N_y"] for rec in kept],
            "N_z": [rec["N_z"] for rec in kept],
        }
    saved = {
        "complete": True,
        "n_done": nreal,
        "n_insert": n_cap,
        "stats": stats,
        "raw": raw,
        "extended": int(sum(rec["extended"] for rec in results)),
        "censored": int(sum(rec["censored"] for rec in results)),
        "seconds": time.perf_counter() - t0,
    }
    CHECKPOINT.mkdir(parents=True, exist_ok=True)
    with ck_path.open("w", encoding="utf-8") as f:
        json.dump(saved, f)
    log(f"done {name} mean={stats['mean']:.6g} se={stats['se']:.3g} n={stats['nreal']} "
        f"discarded={stats['discarded']} extended={saved['extended']} censored={saved['censored']} "
        f"time={saved['seconds']:.1f}s")
    return saved


def warmup():
    L = 10.0
    n = 25
    rng = np.random.default_rng(0)
    pos = rng.uniform(0, L, size=(n, 3))
    nrm = np.zeros((n, 3))
    nrm[:, 2] = 1.0
    e1 = np.zeros((n, 3))
    e1[:, 0] = 1.0
    e2 = np.zeros((n, 3))
    e2[:, 1] = 1.0
    rad = np.ones(n)
    aux = np.zeros(n)
    rb = np.ones(n)
    core.simulate(0, pos, nrm, e1, e2, rad, aux, rb, L, 4, 1)
    aux[:] = 0.1
    rb[:] = math.sqrt(1.0 + 0.01)
    core.simulate(1, pos, nrm, e1, e2, rad, aux, rb, L, 4, 0)
    rb[:] = core.RB_SQUARE
    core.simulate(2, pos, nrm, e1, e2, rad, np.zeros(n), rb, L, 4, 0)
    rb[:] = core.RB_TRIANGLE
    core.simulate(3, pos, nrm, e1, e2, rad, np.zeros(n), rb, L, 4, 0)
    aux[:] = 0.001
    rb[:] = 1.0
    core.simulate(4, pos, nrm, e1, e2, rad, aux, rb, L, 4, 0)
    log("warmup compiled")


def main():
    if PROGRESS.exists():
        PROGRESS.unlink()
    t_start = time.perf_counter()
    warmup()
    log("deterministic pairs")
    pairs = compute_pairs()
    for k, v in pairs.items():
        log(f"  pairs {k}: {sum(v)} / {len(v)} intersections")
    log("certify pairs")
    cert = certify_pairs(pairs)
    log(f"  certify {cert}")
    if any(cert.values()):
        log("CERTIFY FAILURES — stopping before the cases")
        (ROOT / "work" / "certify.json").write_text(json.dumps(cert, indent=2), encoding="utf-8")
        sys.exit(1)
    log("deterministic configs")
    configs, checks = compute_configs()
    for k, v in configs.items():
        log(f"  configs {k}: {v}")
    log(f"  config grid checks {checks}")
    if not all(c["match"] for c in checks.values()):
        log("CONFIG GRID MISMATCH — stopping")
        sys.exit(1)
    cases = json.loads((ROOT / "cases.json").read_text(encoding="utf-8"))["cases"]
    n_workers = max(1, min(20, os.cpu_count() or 1))
    log(f"workers {n_workers}")
    case_out = {}
    for idx, (name, spec) in enumerate(cases.items()):
        case_out[name] = run_case(name, spec, idx, n_workers)
    runtime = time.perf_counter() - t_start
    stats = {name: rec["stats"] for name, rec in case_out.items()}
    raw = {}
    for name, spec in cases.items():
        if spec.get("record_raw"):
            raw[name] = case_out[name]["raw"]
    results = {
        "implementation": {
            "language": "Python 3.14 with numba 0.68 (nopython, nogil) for contact tests, GJK and the insertion loop",
            "rng": "numpy.random.Generator(PCG64) via numpy.random.SeedSequence",
            "seeds": "SeedSequence([20261025, case_index, realisation_index]); case_index is the order in cases.json starting at 0; realisation_index starts at 0",
            "runtime_seconds": round(runtime, 3),
            "notes": (
                "Closed contacts. Thin discs use the intersection-line chord test; "
                "cylinders and spheroids use GJK with analytic support (tolerance 1e-9); "
                "plates use polygon-plane clipping and treat parallel planes as non-intersecting. "
                "Wrapping uses a periodic union-find. "
                f"Directional extensions: { {k: case_out[k]['extended'] for k in case_out} }. "
                f"Censored directional counts: { {k: case_out[k]['censored'] for k in case_out} }."
            ),
        },
        "pairs": pairs,
        "configs": configs,
        "stats": stats,
        "raw": raw,
    }
    out_path = ROOT / "results.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(results, f)
    meta = {
        "runtime_seconds": runtime,
        "certify": cert,
        "config_checks": checks,
        "cases": {k: {"n_insert": case_out[k]["n_insert"], "extended": case_out[k]["extended"],
                      "censored": case_out[k]["censored"], "seconds": case_out[k]["seconds"],
                      "stats": case_out[k]["stats"]} for k in case_out},
    }
    (ROOT / "work" / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log(f"wrote results.json in {runtime:.1f}s")
    import subprocess
    proc = subprocess.run([sys.executable, "validate_results.py", "results.json"],
                          cwd=str(ROOT), capture_output=True, text=True)
    log(proc.stdout.strip())
    if proc.stderr.strip():
        log(proc.stderr.strip())
    log(f"validator exit {proc.returncode}")


if __name__ == "__main__":
    main()
