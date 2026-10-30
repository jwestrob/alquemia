"""Variational, real, closed-shell AO-block absolutely localized orbitals.

The energy and Fock matrix come from the unchanged full-system PySCF backend.
Only occupied-orbital AO support is restricted.  This module does not zero Fock
elements, alter integrals, assign an atomic population constraint, or run any
calculation at import time.  See THEORY.md for equations and limitations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import time
from typing import Any, Callable, Sequence

import numpy as np
from scipy import linalg


class ALMOError(ValueError):
    """An input, rank, or invariant prevents a well-defined ALMO calculation."""


@dataclass
class ALMOEvaluation:
    energy: float
    density: np.ndarray
    fock: np.ndarray
    coefficient_gradient: np.ndarray
    invariants: dict[str, Any]


@dataclass
class ALMOResult:
    coefficients: np.ndarray
    density: np.ndarray
    energy: float
    converged: bool
    status: str
    gradient_norm: float
    coefficient_gradient: np.ndarray
    history: list[dict[str, Any]]
    invariants: dict[str, Any]
    iterations: int
    evaluations: int
    initial_energy: float
    options: dict[str, Any]
    trial_history: list[dict[str, Any]] = field(default_factory=list)
    message: str = ""
    elapsed_seconds: float = 0.0

    @property
    def C(self) -> np.ndarray:
        return self.coefficients

    @property
    def delta_energy(self) -> float:
        return self.energy - self.initial_energy

    def summary(self) -> dict[str, Any]:
        """JSON-compatible scalar/record summary; save arrays separately."""
        return {
            "energy": self.energy, "initial_energy": self.initial_energy,
            "delta_energy": self.delta_energy, "converged": self.converged,
            "status": self.status, "message": self.message,
            "gradient_norm": self.gradient_norm, "iterations": self.iterations,
            "evaluations": self.evaluations, "invariants": self.invariants,
            "options": self.options, "history": self.history,
            "trial_history": self.trial_history,
            "elapsed_seconds": self.elapsed_seconds,
            "stationarity_only": True, "minimum_or_physical_state_certified": False,
        }


def _real_array(value: Any, name: str) -> np.ndarray:
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ALMOError(f"{name}: complex orbitals/backends are not implemented")
    out = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ALMOError(f"{name}: nonfinite values")
    return out


def _symmetric(value: Any, name: str, tolerance: float = 1e-10) -> np.ndarray:
    matrix = _real_array(value, name)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ALMOError(f"{name}: expected one real square restricted-spin matrix")
    if np.max(np.abs(matrix - matrix.T), initial=0.0) > tolerance:
        raise ALMOError(f"{name}: nonsymmetric matrix")
    return (matrix + matrix.T) * 0.5


def _dot(a: Sequence[np.ndarray], b: Sequence[np.ndarray]) -> float:
    return float(sum(np.einsum("ij,ij->", x, y) for x, y in zip(a, b)))


def _norm(a: Sequence[np.ndarray]) -> float:
    return math.sqrt(max(0.0, _dot(a, a)))


def _project(u: Sequence[np.ndarray], v: Sequence[np.ndarray]) -> list[np.ndarray]:
    # Horizontal (Grassmann) tangent: occupied-occupied gauge rotations vanish.
    return [z - q @ (q.T @ z) for q, z in zip(u, v)]


def _polar(block: np.ndarray, rank_tolerance: float) -> np.ndarray:
    if block.shape[1] == 0:
        return block.copy()
    eigenvalues, eigenvectors = linalg.eigh(block.T @ block)
    if eigenvalues[0] <= rank_tolerance * max(1.0, eigenvalues[-1]):
        raise ALMOError("rank-deficient fragment occupied coefficients")
    return block @ ((eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T)


class ALMOProblem:
    """Full-backend energy on an explicit fragment AO-support constraint.

    ``fragment_aos`` must partition all AO indices. Occupied columns are grouped
    in the same fragment order. Fragment occupations are numbers of doubly
    occupied spatial orbitals, not electron counts. An arbitrary (nonorthogonal)
    full-rank supported C can be passed to :meth:`evaluate` for derivative tests.
    No overlap eigenvectors are silently removed.
    """

    def __init__(
        self, mf: Any, fragment_aos: Sequence[Sequence[int]],
        fragment_nocc: Sequence[int], *, overlap_eigenvalue_tolerance: float = 1e-10,
        occupied_eigenvalue_tolerance: float = 1e-10,
        support_tolerance: float = 1e-12, invariant_tolerance: float = 1e-8,
    ) -> None:
        self.mf = mf
        self.mol = mf.mol
        if self.mol.spin != 0 or self.mol.nelectron % 2:
            raise ALMOError("only closed-shell spin=0, even-electron systems are supported")
        if getattr(mf, "with_solvent", None) is not None:
            raise ALMOError("solvent wrappers require a separately validated energy derivative")
        # ROHF/UHF/GHF, custom Fock wrappers and fractional occupations must not
        # be accidentally treated as the intended real closed-shell Hamiltonian.
        from pyscf.scf import hf
        if not isinstance(mf, hf.RHF):
            raise ALMOError("backend must be a PySCF RHF or RKS instance")
        for name, value in {
            "overlap_eigenvalue_tolerance": overlap_eigenvalue_tolerance,
            "occupied_eigenvalue_tolerance": occupied_eigenvalue_tolerance,
            "support_tolerance": support_tolerance,
            "invariant_tolerance": invariant_tolerance,
        }.items():
            if not np.isfinite(value) or value <= 0:
                raise ALMOError(f"{name} must be finite and positive")
        self.overlap_tolerance = float(overlap_eigenvalue_tolerance)
        self.occupied_tolerance = float(occupied_eigenvalue_tolerance)
        self.support_tolerance = float(support_tolerance)
        self.invariant_tolerance = float(invariant_tolerance)
        self.overlap = _symmetric(mf.get_ovlp(), "AO overlap")
        self.hcore = _symmetric(mf.get_hcore(), "core Hamiltonian")
        self.nao = self.overlap.shape[0]
        if self.hcore.shape != self.overlap.shape:
            raise ALMOError("core Hamiltonian/overlap dimensions differ")
        if len(fragment_aos) != len(fragment_nocc) or not fragment_aos:
            raise ALMOError("fragment AO and occupation lists must have equal nonzero length")
        self.fragment_aos: list[np.ndarray] = []
        self.fragment_nocc: list[int] = []
        for indices, occupation in zip(fragment_aos, fragment_nocc):
            raw = np.asarray(indices)
            if raw.ndim != 1 or not len(raw) or raw.dtype.kind not in "iu":
                raise ALMOError("each fragment needs a nonempty integer AO index list")
            if isinstance(occupation, (bool, np.bool_)) or int(occupation) != occupation:
                raise ALMOError("fragment_nocc entries must be integers")
            nocc = int(occupation)
            if not 0 <= nocc <= len(raw):
                raise ALMOError("fragment occupation lies outside its AO dimension")
            self.fragment_aos.append(raw.astype(int, copy=True))
            self.fragment_nocc.append(nocc)
        flattened = np.concatenate(self.fragment_aos)
        if sorted(flattened.tolist()) != list(range(self.nao)):
            raise ALMOError("fragment_aos must partition every AO exactly once")
        self.nocc = sum(self.fragment_nocc)
        if self.nocc == 0 or 2 * self.nocc != self.mol.nelectron:
            raise ALMOError("fragment occupied counts do not equal the molecular electron count")
        offset = 0
        self.column_slices = []
        self.support = np.zeros((self.nao, self.nocc), dtype=bool)
        self.whiteners, self.inverse_whiteners = [], []
        self.fragment_overlap_eigenvalues = []
        for indices, nocc in zip(self.fragment_aos, self.fragment_nocc):
            columns = slice(offset, offset + nocc)
            self.column_slices.append(columns)
            self.support[indices, columns] = True
            offset += nocc
            block_overlap = self.overlap[np.ix_(indices, indices)]
            values, vectors = linalg.eigh(block_overlap)
            if values[0] <= self.overlap_tolerance * max(1.0, values[-1]):
                raise ALMOError("fragment AO overlap is singular/ill-conditioned; no basis truncation performed")
            self.whiteners.append(vectors / np.sqrt(values))
            self.inverse_whiteners.append(np.sqrt(values)[:, None] * vectors.T)
            self.fragment_overlap_eigenvalues.append((float(values[0]), float(values[-1])))
        # Positive definite full overlap is required for the unmodified AO space.
        self.overlap_extrema = (
            float(linalg.eigh(self.overlap, subset_by_index=[0, 0], eigvals_only=True)[0]),
            float(linalg.eigh(self.overlap, subset_by_index=[self.nao - 1, self.nao - 1], eigvals_only=True)[0]),
        )
        if self.overlap_extrema[0] <= self.overlap_tolerance * max(1.0, self.overlap_extrema[1]):
            raise ALMOError("full AO overlap is singular/ill-conditioned; no basis truncation performed")
        # Build a fixed configured grid without density-dependent point removal.
        # An existing grid is retained exactly, which lets callers share the
        # ordinary-reference grid. Never rebuild grids during optimization.
        for name in ("grids", "nlcgrids"):
            grid = getattr(mf, name, None)
            if grid is not None and grid.coords is None:
                if name == "nlcgrids" and not (hasattr(mf, "do_nlc") and mf.do_nlc()):
                    continue
                grid.build(with_non0tab=True)
        self.backend_xc = getattr(mf, "xc", None)
        self.evaluations = 0

    def assemble_coefficients(self, blocks: Sequence[np.ndarray]) -> np.ndarray:
        """Embed local AO coefficient blocks in full AO/fragment-column order."""
        if len(blocks) != len(self.fragment_aos):
            raise ALMOError("wrong number of fragment coefficient blocks")
        coefficients = np.zeros((self.nao, self.nocc))
        for block, indices, columns, nocc in zip(
            blocks, self.fragment_aos, self.column_slices, self.fragment_nocc
        ):
            block = _real_array(block, "fragment coefficients")
            if block.shape != (len(indices), nocc):
                raise ALMOError("fragment coefficient shape does not match AO and occupied counts")
            coefficients[indices, columns] = block
        return coefficients

    def _check_coefficients(self, coefficients: Any) -> np.ndarray:
        coefficients = _real_array(coefficients, "occupied coefficients")
        if coefficients.shape != (self.nao, self.nocc):
            raise ALMOError(f"occupied coefficients need shape {(self.nao, self.nocc)}")
        if np.max(np.abs(coefficients[~self.support]), initial=0.0) > self.support_tolerance:
            raise ALMOError("occupied coefficient violates its fragment AO support")
        return coefficients

    def initial_coefficients(
        self, initial_coefficients: Any = None, fragment_mfs: Sequence[Any] | None = None,
    ) -> np.ndarray:
        """Normalize supported initial orbitals without changing their subspaces.

        Prefer converged isolated-fragment orbitals with the exact fragment AO
        ordering. If neither input is supplied, the full-system core Hamiltonian
        restricted to each fragment supplies an explicit core-only guess.
        """
        if initial_coefficients is not None and fragment_mfs is not None:
            raise ALMOError("provide initial_coefficients or fragment_mfs, not both")
        if fragment_mfs is not None:
            if len(fragment_mfs) != len(self.fragment_aos):
                raise ALMOError("wrong number of fragment SCF objects")
            blocks = []
            for fragment, indices, nocc in zip(fragment_mfs, self.fragment_aos, self.fragment_nocc):
                if not fragment.converged:
                    raise ALMOError("initial isolated-fragment SCF is not converged")
                if fragment.mol.spin != 0 or fragment.mol.nelectron != 2 * nocc:
                    raise ALMOError("initial fragment electron/spin assignment mismatch")
                occupations = _real_array(fragment.mo_occ, "fragment occupations")
                if occupations.ndim != 1 or not np.all((occupations == 0) | (occupations == 2)):
                    raise ALMOError("fractional or unrestricted fragment occupations are not supported")
                block = _real_array(fragment.mo_coeff, "fragment MO coefficients")[:, occupations == 2]
                if block.shape != (len(indices), nocc):
                    raise ALMOError("initial fragment orbital shape mismatch")
                fragment_overlap = _symmetric(fragment.get_ovlp(), "fragment overlap")
                if not np.allclose(fragment_overlap, self.overlap[np.ix_(indices, indices)], atol=1e-10, rtol=1e-10):
                    raise ALMOError("fragment and full-system AO overlap/order differ")
                blocks.append(block)
            initial_coefficients = self.assemble_coefficients(blocks)
        elif initial_coefficients is None:
            blocks = []
            for indices, nocc, whitener in zip(self.fragment_aos, self.fragment_nocc, self.whiteners):
                transformed = whitener.T @ self.hcore[np.ix_(indices, indices)] @ whitener
                _, vectors = linalg.eigh(transformed)
                blocks.append(whitener @ vectors[:, :nocc])
            initial_coefficients = self.assemble_coefficients(blocks)
        elif isinstance(initial_coefficients, (list, tuple)):
            initial_coefficients = self.assemble_coefficients(initial_coefficients)
        coefficients = self._check_coefficients(initial_coefficients)
        return self.from_whitened(self.to_whitened(coefficients, normalize=True))

    def to_whitened(self, coefficients: np.ndarray, *, normalize: bool = False) -> list[np.ndarray]:
        coefficients = self._check_coefficients(coefficients)
        blocks = [inverse @ coefficients[indices, columns]
                  for inverse, indices, columns in zip(
                      self.inverse_whiteners, self.fragment_aos, self.column_slices)]
        return [_polar(block, self.occupied_tolerance) for block in blocks] if normalize else blocks

    def from_whitened(self, blocks: Sequence[np.ndarray]) -> np.ndarray:
        return self.assemble_coefficients([x @ u for x, u in zip(self.whiteners, blocks)])

    def projected_gradient(self, coefficients: np.ndarray, evaluation: ALMOEvaluation) -> list[np.ndarray]:
        blocks = self.to_whitened(coefficients)
        gradient = [x.T @ evaluation.coefficient_gradient[indices, columns]
                    for x, indices, columns in zip(
                        self.whiteners, self.fragment_aos, self.column_slices)]
        return _project(blocks, gradient)

    def evaluate(self, coefficients: np.ndarray) -> ALMOEvaluation:
        """Evaluate the full energy and exact analytic dE/dC (Hartree).

        ``coefficient_gradient`` includes forbidden directions for derivative
        transparency. Only entries within ``problem.support`` are variational.
        Energies include the backend's nuclear repulsion. Each call builds one
        full potential; no SCF kernel or DIIS modification is used here.
        """
        coefficients = self._check_coefficients(coefficients)
        if self.mf.mol is not self.mol or getattr(self.mf, "xc", None) != self.backend_xc:
            raise ALMOError("backend molecule/functional changed after problem construction")
        sc = self.overlap @ coefficients
        metric = (coefficients.T @ sc)
        metric = (metric + metric.T) * 0.5
        eigenvalues = linalg.eigvalsh(metric)
        if eigenvalues[0] <= self.occupied_tolerance * max(1.0, eigenvalues[-1]):
            raise ALMOError("combined occupied overlap is singular/ill-conditioned")
        chol = linalg.cho_factor(metric, lower=True, check_finite=False)
        dual = linalg.cho_solve(chol, coefficients.T, check_finite=False).T
        density = 2.0 * dual @ coefficients.T
        density = (density + density.T) * 0.5
        electron_count = float(np.einsum("ij,ji->", density, self.overlap))
        # Evaluate the actual density identity with thin occupied factors.
        residual = 2.0 * (density @ sc) @ dual.T - 2.0 * density
        idempotency_error = float(linalg.norm(residual) / max(1.0, linalg.norm(2.0 * density)))
        support_error = float(np.max(np.abs(coefficients[~self.support]), initial=0.0))
        invariants = {
            "electron_count": electron_count, "expected_electron_count": int(self.mol.nelectron),
            "electron_count_error": abs(electron_count - self.mol.nelectron),
            "density_idempotency_relative_error": idempotency_error,
            "coefficient_support_max_abs": support_error,
            "occupied_overlap_min_eigenvalue": float(eigenvalues[0]),
            "occupied_overlap_max_eigenvalue": float(eigenvalues[-1]),
            "occupied_overlap_condition_number": float(eigenvalues[-1] / eigenvalues[0]),
            "ao_overlap_min_eigenvalue": self.overlap_extrema[0],
            "ao_overlap_max_eigenvalue": self.overlap_extrema[1],
            "fragment_nocc": list(self.fragment_nocc),
            "support_preserved": support_error <= self.support_tolerance,
            "rank_preserved": True,
        }
        if (invariants["electron_count_error"] > self.invariant_tolerance
                or idempotency_error > self.invariant_tolerance):
            raise ALMOError(f"density invariant failure: {invariants}")
        invariants["passed"] = True
        self.evaluations += 1
        potential = self.mf.get_veff(self.mol, density, hermi=1)
        # Preserve potential's PySCF energy tags when calling energy_tot.
        energy = float(self.mf.energy_tot(dm=density, h1e=self.hcore, vhf=potential))
        fock = _symmetric(self.hcore + np.asarray(potential), "full Fock")
        if not math.isfinite(energy):
            raise ALMOError("backend returned a nonfinite energy")
        fdual = fock @ dual
        gradient = 4.0 * (fdual - (self.overlap @ dual) @ (coefficients.T @ fdual))
        if not np.all(np.isfinite(gradient)):
            raise ALMOError("nonfinite analytic orbital gradient")
        return ALMOEvaluation(energy, density, fock, gradient, invariants)


def _preconditioner(
    problem: ALMOProblem, u: Sequence[np.ndarray], fock: np.ndarray,
    method: str, gap_floor: float,
) -> Callable[[Sequence[np.ndarray]], list[np.ndarray]]:
    if method == "identity":
        return lambda value: [v.copy() for v in value]
    spectra = []
    for indices, x, occupied in zip(problem.fragment_aos, problem.whiteners, u):
        transformed = x.T @ fock[np.ix_(indices, indices)] @ x
        transformed = (transformed + transformed.T) * 0.5
        all_energies, all_vectors = linalg.eigh(transformed)
        occupied_energies, occupied_vectors = linalg.eigh(occupied.T @ transformed @ occupied)
        denominators = 4.0 * np.maximum(
            np.abs(all_energies[:, None] - occupied_energies[None, :]), gap_floor)
        spectra.append((all_vectors, occupied_vectors, denominators))

    def apply(value: Sequence[np.ndarray]) -> list[np.ndarray]:
        result = [v @ ((v.T @ g @ r) / denominator) @ r.T
                  for g, (v, r, denominator) in zip(value, spectra)]
        return _project(u, result)

    return apply


def solve_almo(
    mf: Any, fragment_aos: Sequence[Sequence[int]], fragment_nocc: Sequence[int],
    initial_coefficients: Any = None, fragment_mfs: Sequence[Any] | None = None, *,
    max_iterations: int = 300, max_evaluations: int = 1500,
    gradient_tolerance: float = 1e-6, energy_tolerance: float = 1e-10,
    overlap_eigenvalue_tolerance: float = 1e-10,
    occupied_eigenvalue_tolerance: float = 1e-10,
    support_tolerance: float = 1e-12, invariant_tolerance: float = 1e-8,
    memory: int = 12, preconditioner: str = "orbital_gap", gap_floor: float = 0.1,
    max_step_norm: float = 0.5, armijo: float = 1e-4,
    backtrack_factor: float = 0.5, max_backtracks: int = 24,
    energy_roundoff_tolerance: float = 1e-12,
    callback: Callable[[dict[str, Any]], None] | None = None,
) -> ALMOResult:
    """Minimize the exact full-system energy over a product of Grassmann spaces.

    Polar retraction preserves each fragment's orthonormal occupied frame in its
    local AO metric. Projected-transport L-BFGS uses a positive orbital-gap
    preconditioner and bounded Armijo backtracking. A failed quasi-Newton line
    search gets one explicit preconditioned-gradient attempt with cleared
    history. Backend exceptions, exhausted budgets, or a second failed search
    return a nonconverged last-accepted state; invalid initial inputs raise.

    A converged result certifies first-order constrained stationarity only.
    It does not certify a global/local minimum or a physical electronic state.
    """
    started = time.monotonic()
    options = {
        "max_iterations": max_iterations, "max_evaluations": max_evaluations,
        "gradient_tolerance": gradient_tolerance, "energy_tolerance": energy_tolerance,
        "overlap_eigenvalue_tolerance": overlap_eigenvalue_tolerance,
        "occupied_eigenvalue_tolerance": occupied_eigenvalue_tolerance,
        "support_tolerance": support_tolerance, "invariant_tolerance": invariant_tolerance,
        "memory": memory, "preconditioner": preconditioner, "gap_floor": gap_floor,
        "max_step_norm": max_step_norm, "armijo": armijo,
        "backtrack_factor": backtrack_factor, "max_backtracks": max_backtracks,
        "energy_roundoff_tolerance": energy_roundoff_tolerance,
    }
    for name in ("max_iterations", "max_evaluations", "memory", "max_backtracks"):
        value = options[name]
        if isinstance(value, bool) or int(value) != value or value < (0 if name == "max_iterations" else 1):
            raise ALMOError(f"{name} must be a valid nonnegative/positive integer")
    for name in ("gradient_tolerance", "energy_tolerance", "gap_floor", "max_step_norm", "energy_roundoff_tolerance"):
        if not np.isfinite(options[name]) or options[name] <= 0:
            raise ALMOError(f"{name} must be finite and positive")
    if preconditioner not in ("identity", "orbital_gap"):
        raise ALMOError("preconditioner must be identity or orbital_gap")
    if not 0 < armijo < 1 or not 0 < backtrack_factor < 1:
        raise ALMOError("Armijo and backtracking factors must lie between zero and one")
    problem = ALMOProblem(
        mf, fragment_aos, fragment_nocc,
        overlap_eigenvalue_tolerance=overlap_eigenvalue_tolerance,
        occupied_eigenvalue_tolerance=occupied_eigenvalue_tolerance,
        support_tolerance=support_tolerance, invariant_tolerance=invariant_tolerance,
    )
    coefficients = problem.initial_coefficients(initial_coefficients, fragment_mfs)
    u = problem.to_whitened(coefficients)
    evaluation = problem.evaluate(coefficients)
    initial_energy = evaluation.energy
    gradient = problem.projected_gradient(coefficients, evaluation)
    history: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    # s,y pairs share the current horizontal tangent space via projection.
    pairs: list[tuple[list[np.ndarray], list[np.ndarray]]] = []
    status, message, converged = "maximum_iterations", "iteration budget exhausted", False
    last_energy_change = 0.0

    def record(iteration: int, step: float, direction_kind: str) -> None:
        row = {
            "iteration": iteration, "evaluations": problem.evaluations,
            "energy": evaluation.energy, "energy_change": last_energy_change,
            "gradient_norm": _norm(gradient), "step": step,
            "direction": direction_kind, "invariants": dict(evaluation.invariants),
            "elapsed_seconds": time.monotonic() - started,
        }
        history.append(row)
        if callback is not None:
            callback(dict(row))

    record(0, 0.0, "initial")
    for iteration in range(max_iterations + 1):
        gradient_norm = _norm(gradient)
        if gradient_norm <= gradient_tolerance and abs(last_energy_change) <= energy_tolerance:
            status, message, converged = "converged", "constrained orbital gradient and energy-change criteria passed", True
            break
        if iteration == max_iterations:
            break
        if problem.evaluations >= max_evaluations:
            status, message = "maximum_evaluations", "full-backend evaluation budget exhausted"
            break
        apply_h0 = _preconditioner(problem, u, evaluation.fock, preconditioner, gap_floor)
        q = [g.copy() for g in gradient]
        alphas = []
        for s, y in reversed(pairs):
            alpha = _dot(s, q) / _dot(s, y)
            alphas.append(alpha)
            q = [a - alpha * b for a, b in zip(q, y)]
        direction = apply_h0(q)
        for (s, y), alpha in zip(pairs, reversed(alphas)):
            beta = _dot(y, direction) / _dot(s, y)
            direction = [a + (alpha - beta) * b for a, b in zip(direction, s)]
        direction = _project(u, [-v for v in direction])
        if _dot(gradient, direction) >= -1e-14 * max(1.0, gradient_norm * _norm(direction)):
            pairs.clear()
            direction = [-v for v in apply_h0(gradient)]
        accepted = None
        fatal = False
        for attempt in range(2):
            kind = "lbfgs" if attempt == 0 and pairs else "preconditioned_gradient"
            if attempt:
                if kind == "preconditioned_gradient" and not pairs:
                    # There is no different direction to retry when the first
                    # attempt already used steepest/preconditioned descent.
                    break
                pairs.clear()
                direction = [-v for v in apply_h0(gradient)]
                kind = "preconditioned_gradient_restart"
            direction_norm = _norm(direction)
            if direction_norm > max_step_norm:
                direction = [v * (max_step_norm / direction_norm) for v in direction]
            slope = _dot(gradient, direction)
            if not math.isfinite(slope) or slope >= 0:
                status, message, fatal = "non_descent_direction", "no finite strict descent direction", True
                break
            step = 1.0
            for backtrack in range(max_backtracks):
                if problem.evaluations >= max_evaluations:
                    status, message, fatal = "maximum_evaluations", "full-backend evaluation budget exhausted", True
                    break
                trial = {
                    "iteration": iteration + 1, "attempt": attempt, "backtrack": backtrack,
                    "step": step, "direction": kind, "slope": slope, "accepted": False,
                }
                try:
                    new_u = [_polar(a + step * d, occupied_eigenvalue_tolerance)
                             for a, d in zip(u, direction)]
                    new_coefficients = problem.from_whitened(new_u)
                    new_evaluation = problem.evaluate(new_coefficients)
                    trial["energy"] = new_evaluation.energy
                    trial["evaluations"] = problem.evaluations
                except ALMOError as exc:
                    trial["error"] = str(exc)
                    trials.append(trial)
                    # An inadmissible occupied metric is a rejected trial step;
                    # all other invariant/backend defects terminate the solve.
                    if "combined occupied overlap" in str(exc) or "rank-deficient fragment" in str(exc):
                        step *= backtrack_factor
                        continue
                    status, message, fatal = "evaluation_failed", str(exc), True
                    break
                except Exception as exc:
                    trial["error"] = f"{type(exc).__name__}: {exc}"
                    trials.append(trial)
                    status, message, fatal = "backend_failed", trial["error"], True
                    break
                slack = max(energy_roundoff_tolerance,
                            8.0 * np.finfo(float).eps * max(1.0, abs(evaluation.energy)))
                trial["roundoff_allowance"] = slack
                trial["armijo_bound"] = evaluation.energy + armijo * step * slope + slack
                trial["accepted"] = bool(new_evaluation.energy <= trial["armijo_bound"])
                trials.append(trial)
                if trial["accepted"]:
                    accepted = (new_u, new_coefficients, new_evaluation, step, direction, kind)
                    break
                step *= backtrack_factor
            if accepted is not None or fatal:
                break
        if fatal:
            break
        if accepted is None:
            status, message = "line_search_failed", "bounded Armijo search failed; last accepted state retained"
            break
        new_u, new_coefficients, new_evaluation, step, direction, kind = accepted
        new_gradient = problem.projected_gradient(new_coefficients, new_evaluation)
        transported_gradient = _project(new_u, gradient)
        s = _project(new_u, [step * d for d in direction])
        y = [new - old for new, old in zip(new_gradient, transported_gradient)]
        transported_pairs = []
        for old_s, old_y in pairs:
            next_s, next_y = _project(new_u, old_s), _project(new_u, old_y)
            if _dot(next_s, next_y) > 1e-10 * _norm(next_s) * _norm(next_y):
                transported_pairs.append((next_s, next_y))
        if _dot(s, y) > 1e-10 * _norm(s) * _norm(y):
            transported_pairs.append((s, y))
        pairs = transported_pairs[-memory:]
        last_energy_change = new_evaluation.energy - evaluation.energy
        u, coefficients, evaluation, gradient = new_u, new_coefficients, new_evaluation, new_gradient
        record(iteration + 1, step, kind)
    return ALMOResult(
        coefficients=coefficients.copy(), density=evaluation.density.copy(), energy=evaluation.energy,
        converged=converged, status=status, gradient_norm=_norm(gradient),
        coefficient_gradient=evaluation.coefficient_gradient.copy(), history=history,
        invariants=dict(evaluation.invariants), iterations=len(history) - 1,
        evaluations=problem.evaluations, initial_energy=initial_energy, options=options,
        trial_history=trials, message=message, elapsed_seconds=time.monotonic() - started,
    )


def release_constraints(
    mf: Any, result: ALMOResult, *, conv_tol: float = 1e-10,
    conv_tol_grad: float = 1e-6, max_cycle: int = 200,
) -> dict[str, Any]:
    """Run ordinary SCF from the ALMO density; this explicitly mutates ``mf``.

    Use a separately prepared same-Hamiltonian backend if retaining its earlier
    SCF attributes matters. This helper does not change basis, method, grids,
    solver class, or apply a fallback. An unconverged release remains failed.
    """
    if result.density.shape != (mf.mol.nao_nr(), mf.mol.nao_nr()):
        raise ALMOError("release backend density dimension mismatch")
    mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = conv_tol, conv_tol_grad, max_cycle
    energy = float(mf.kernel(dm0=result.density.copy()))
    finite = math.isfinite(energy)
    return {
        "energy": energy, "converged": bool(mf.converged and finite),
        "status": "converged" if mf.converged and finite else "unconverged",
        "energy_change_from_almo": energy - result.energy,
        "density": np.asarray(mf.make_rdm1()).copy(),
        "coefficients": np.asarray(mf.mo_coeff).copy(),
        "occupations": np.asarray(mf.mo_occ).copy(),
        "conv_tol": conv_tol, "conv_tol_grad": conv_tol_grad, "max_cycle": max_cycle,
    }
