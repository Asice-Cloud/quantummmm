import numpy as np
import matplotlib.pyplot as plt

"""Toy honeycomb Majorana model with a Z2 gauge field and Berry holonomy.

This script is a *toy* scaffold for implementing the vison Berry experiment
on a finite Kitaev honeycomb lattice, mirroring the structure used in
run_pip_vortex_berry.py for the 2D p+ip case.

Current goals of this scaffold:

  1. Provide a concrete brick-wall representation of a finite honeycomb
     lattice with three types of bonds (x, y, z) and a static Z2 gauge
     field u_ij = ±1 on each bond.
  2. Build the free-Majorana Hamiltonian

        H = (i/4) \sum_{<ij>_\alpha} J_\alpha u_{ij} c_i c_j ,

     as an N x N Hermitian matrix (in the Majorana basis).
  3. Given a *sequence* of gauge configurations u^{(k)}, compute at each step
     the low-energy subspace and the non-Abelian Berry holonomy obtained by
     parallel transport around the loop in gauge-configuration space.

IMPORTANT: The concrete choice of gauge path here is deliberately simple and
serves only as a template. To turn this into a true "two-vison around each
other" experiment, one should implement a branch-cut construction on the
honeycomb plaquettes as outlined in kit4-3 §4.3.8, and replace
`toy_gauge_path` accordingly.
"""


def honeycomb_indices(Lx: int, Ly: int):
    """Return site index maps for a brick-wall representation of honeycomb.

    We take a rectangular array of Lx x Ly unit cells. Each unit cell has
    two sites: A (sublattice 0) and B (sublattice 1).

    We label sites by a single integer index i in [0, 2*Lx*Ly-1], with

        i = 2 * (y * Lx + x) + s,

    where x in [0, Lx-1], y in [0, Ly-1], and s in {0, 1} for A/B.
    """
    idx_of = {}
    pos_of = {}
    idx = 0
    for y in range(Ly):
        for x in range(Lx):
            for s in (0, 1):  # 0: A, 1: B
                idx_of[(x, y, s)] = idx
                pos_of[idx] = (x, y, s)
                idx += 1
    return idx_of, pos_of


def build_honeycomb_bonds(Lx: int, Ly: int):
    """Build nearest-neighbor bonds (i, j, alpha) for a brick-wall honeycomb.

    Convention (one of many equivalent brick-wall realizations):

      - z-bond: A(x, y) -- B(x, y)
      - x-bond: B(x, y) -- A(x+1, y)
      - y-bond: B(x, y) -- A(x, y+1)

    with open boundaries. alpha in {"x", "y", "z"} labels the bond type.
    """
    idx_of, _ = honeycomb_indices(Lx, Ly)
    bonds = []  # list of (i, j, alpha)

    for y in range(Ly):
        for x in range(Lx):
            iA = idx_of[(x, y, 0)]  # A sublattice
            iB = idx_of[(x, y, 1)]  # B sublattice

            # z-bond within the unit cell
            bonds.append((iA, iB, "z"))

            # x-bond: B(x, y) -- A(x+1, y)
            if x + 1 < Lx:
                iA_right = idx_of[(x + 1, y, 0)]
                bonds.append((iB, iA_right, "x"))

            # y-bond: B(x, y) -- A(x, y+1)
            if y + 1 < Ly:
                iA_up = idx_of[(x, y + 1, 0)]
                bonds.append((iB, iA_up, "y"))

    return bonds


def honeycomb_adjacency(Lx: int, Ly: int):
    """Adjacency list of the honeycomb graph (Majorana sites as nodes).

    Uses the same brick-wall bond convention as build_honeycomb_bonds.
    Returns a dict i -> sorted list of neighbor indices.
    """
    bonds = build_honeycomb_bonds(Lx, Ly)
    N = 2 * Lx * Ly
    adj = {i: set() for i in range(N)}
    for (i, j, _alpha) in bonds:
        adj[i].add(j)
        adj[j].add(i)
    # sort for deterministic behavior
    return {i: sorted(nbs) for i, nbs in adj.items()}


def _canonical_cycle6(path: list[int]) -> tuple[int, ...]:
    """Return a canonical representative for a length-6 cycle.

    We consider cycles equivalent up to rotation and reversal. This
    helper returns the lexicographically smallest tuple among all such
    representatives, to allow de-duplication in a set.
    """
    assert len(path) == 6
    n = 6
    best = None
    seq = path[:]  # copy
    for k in range(n):
        rot = seq[k:] + seq[:k]
        tup_fwd = tuple(rot)
        tup_rev = tuple(reversed(rot))
        if best is None or tup_fwd < best:
            best = tup_fwd
        if tup_rev < best:
            best = tup_rev
    assert best is not None
    return best


def honeycomb_plaquettes(Lx: int, Ly: int) -> list[list[int]]:
    """Enumerate elementary hexagonal plaquettes as 6-cycles in the graph.

    We treat the Majorana lattice as an undirected graph and search for
    simple cycles of length 6 using a depth-limited DFS, then de-duplicate
    them using a canonical representation. This is sufficient for the
    small Lx, Ly systems used in our toy experiments.
    """
    adj = honeycomb_adjacency(Lx, Ly)
    N = 2 * Lx * Ly
    cycles: set[tuple[int, ...]] = set()

    for start in range(N):
        stack: list[tuple[int, list[int]]] = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) >= 7:
                continue
            for nb in adj[node]:
                if nb == start and len(path) == 6:
                    cyc = path[:]  # length-6 simple cycle
                    canon = _canonical_cycle6(cyc)
                    cycles.add(canon)
                elif nb not in path and len(path) < 6:
                    stack.append((nb, path + [nb]))

    # Convert to list-of-lists
    return [list(c) for c in sorted(cycles)]


def build_majorana_matrix(Lx: int, Ly: int,
                          Jx: float, Jy: float, Jz: float,
                          u_ij: dict[tuple[int, int], int]):
    """Construct the antisymmetric matrix A_ij for the Majorana Hamiltonian.

    Hamiltonian (Majorana basis):

        H = (i/4) * sum_{<ij>_alpha} J_alpha * u_ij * c_i c_j.

    We encode this as an N x N Hermitian matrix H_mat via

        H_mat = 0.5j * A,

    where A is real antisymmetric with A_{ij} = J_alpha * u_ij on a bond.
    """
    bonds = build_honeycomb_bonds(Lx, Ly)
    N = 2 * Lx * Ly

    A = np.zeros((N, N), dtype=float)

    for (i, j, alpha) in bonds:
        if alpha == "x":
            J = Jx
        elif alpha == "y":
            J = Jy
        elif alpha == "z":
            J = Jz
        else:
            raise ValueError(f"Unknown bond type alpha={alpha}")

        key = (i, j) if i < j else (j, i)
        u = u_ij.get(key, 1)
        val = float(J * u)

        if i < j:
            A[i, j] += val
            A[j, i] -= val
        else:
            A[j, i] += val
            A[i, j] -= val

    H_mat = 0.5j * A
    return H_mat


def plaquette_fluxes(Lx: int, Ly: int,
                     u_ij: dict[tuple[int, int], int]) -> list[int]:
    """Compute Z2 flux W_p = prod_{(ij) in p} u_ij for each hexagonal plaquette.

    The plaquettes are represented as 6-cycles of site indices. For each
    directed edge (i_k -> i_{k+1}) we use the bond key (min(i,j), max(i,j))
    to look up u_ij, with default +1 if not explicitly present.
    """
    plaquettes = honeycomb_plaquettes(Lx, Ly)
    fluxes: list[int] = []
    for cyc in plaquettes:
        W = 1
        for k in range(6):
            i = cyc[k]
            j = cyc[(k + 1) % 6]
            key = (i, j) if i < j else (j, i)
            u = u_ij.get(key, 1)
            W *= int(u)
        fluxes.append(W)
    return fluxes


def low_energy_subspace(H: np.ndarray, dim_sub: int = 2):
    """Return eigenvalues and eigenvectors of the lowest-|E| states of H.

    H is Hermitian (N x N). We diagonalize, sort by |E|, and return the
    first dim_sub eigenvalues and eigenvectors (columns of V_sel).
    """
    w, v = np.linalg.eigh(H)
    idx = np.argsort(np.abs(w))
    sel = idx[:dim_sub]
    w_sel = w[sel]
    V_sel = v[:, sel]
    return w_sel, V_sel


def project_to_unitary(W: np.ndarray) -> np.ndarray:
    """Project an overlap matrix W to the closest unitary via polar decomposition.

    W = U S V^\dagger -> U_polar = U V^\dagger.
    """
    U, s, Vh = np.linalg.svd(W)
    return U @ Vh


def normalize_to_su2(U: np.ndarray) -> np.ndarray:
    """Strip global phase so that det(U) = 1 (SU(2) representative).

    For a 2x2 unitary U, divide by sqrt(det U); either branch differs by
    a global phase, irrelevant for gate equivalence.
    """
    det = np.linalg.det(U)
    sqrt_det = np.sqrt(det)
    return U / sqrt_det


def ising_R_and_dehn_su2() -> tuple[np.ndarray, np.ndarray]:
    """Return SU(2)-normalized Ising R^{σσ} and its square (Dehn twist)."""
    R1 = np.exp(-0.125j * np.pi)
    Rpsi = np.exp(0.375j * np.pi)
    R = np.diag([R1, Rpsi])
    R2 = np.diag([R1**2, Rpsi**2])
    R_su2 = normalize_to_su2(R)
    R2_su2 = normalize_to_su2(R2)
    return R_su2, R2_su2


def su2_fidelity(U_target: np.ndarray, U_candidate: np.ndarray) -> tuple[float, complex]:
    """Phase-agnostic overlap F = |Tr(U_target^† U_candidate)| / 2 for 2x2 SU(2)."""
    M = U_target.conj().T @ U_candidate
    tr = np.trace(M)
    F = np.abs(tr) / 2.0
    return F, tr


def initial_gauge(Lx: int, Ly: int) -> dict[tuple[int, int], int]:
    """Trivial Z2 gauge field: u_ij = +1 on all bonds.

    This corresponds to the simplest flux-free sector in this toy scaffold.
    A more faithful implementation would choose u_ij to realize W_p = +1 on
    every plaquette in a specific gauge.
    """
    bonds = build_honeycomb_bonds(Lx, Ly)
    u = {}
    for (i, j, alpha) in bonds:
        key = (i, j) if i < j else (j, i)
        u[key] = 1
    return u


def toy_gauge_path(Lx: int, Ly: int, n_steps: int = 20):
    """Construct a simple closed path in Z2 gauge-field space.

    This is *not yet* a true two-vison braiding path. Instead, we take a
    small set of bonds (a short "cut" segment) and flip their u_ij values
    gradually along a loop in parameter space. The purpose is to provide a
    working example of the Berry-holonomy computation pipeline.

    To upgrade this to a vison-braiding path, implement the branch-cut
    construction from kit4-3 §4.3.8: choose two plaquettes as vison cores
    and move one vison around the other by updating u_ij along a sequence
    of cuts.
    """
    bonds = build_honeycomb_bonds(Lx, Ly)
    if not bonds:
        raise ValueError("No bonds constructed; check Lx, Ly")

    # Choose a small subset of bonds to form a toy cut (first few bonds).
    cut_edges = []
    for (i, j, alpha) in bonds:
        key = (i, j) if i < j else (j, i)
        cut_edges.append(key)
        if len(cut_edges) >= 4:
            break

    configs = []
    u0 = initial_gauge(Lx, Ly)

    # Half loop: gradually flip all bonds in cut_edges to -1
    for k in range(n_steps):
        u = dict(u0)
        # For a toy path we simply flip the first m bonds where m grows.
        m = (k * len(cut_edges)) // n_steps
        for e in cut_edges[:m + 1]:
            u[e] = -u[e]
        configs.append(u)

    # Second half: undo the flips to return to u0, closing the loop
    for k in range(n_steps):
        u = dict(u0)
        m = (n_steps - 1 - k) * len(cut_edges) // n_steps
        for e in cut_edges[:m + 1]:
            u[e] = -u[e]
        configs.append(u)

    # Ensure loop closes exactly at u0
    configs.append(u0)
    return configs


def plaquette_centers(Lx: int, Ly: int) -> list[tuple[float, float]]:
    """Compute geometric centers (x,y) for each hexagonal plaquette.

    We average the brick-wall lattice coordinates (x,y) of the six sites
    belonging to a plaquette, ignoring the sublattice index.
    The order of returned centers matches honeycomb_plaquettes(Lx, Ly).
    """
    _, pos_of = honeycomb_indices(Lx, Ly)
    plaquettes = honeycomb_plaquettes(Lx, Ly)
    centers: list[tuple[float, float]] = []
    for cyc in plaquettes:
        xs = []
        ys = []
        for site in cyc:
            x, y, s = pos_of[site]
            xs.append(float(x))
            ys.append(float(y))
        cx = float(np.mean(xs))
        cy = float(np.mean(ys))
        centers.append((cx, cy))
    return centers


def plaquette_adjacency(Lx: int, Ly: int) -> list[list[int]]:
    """Adjacency graph between plaquettes (share a bond).

    Returns an adjacency list adj[p] = list of neighboring plaquette indices.
    Two plaquettes are adjacent if they share at least one bond (i,j).
    """
    plaquettes = honeycomb_plaquettes(Lx, Ly)
    # Map bond -> list of incident plaquettes
    bond_to_plaq: dict[tuple[int, int], list[int]] = {}
    for p_idx, cyc in enumerate(plaquettes):
        for k in range(6):
            i = cyc[k]
            j = cyc[(k + 1) % 6]
            key = (i, j) if i < j else (j, i)
            bond_to_plaq.setdefault(key, []).append(p_idx)

    n_p = len(plaquettes)
    adj: list[set[int]] = [set() for _ in range(n_p)]
    for _bond, inc in bond_to_plaq.items():
        if len(inc) == 2:
            a, b = inc
            adj[a].add(b)
            adj[b].add(a)

    return [sorted(list(s)) for s in adj]


def shortest_plaquette_path(Lx: int, Ly: int, p_start: int, p_end: int) -> list[int]:
    """Breadth-first search for a shortest path between two plaquettes."""
    if p_start == p_end:
        return [p_start]
    adj = plaquette_adjacency(Lx, Ly)
    from collections import deque

    n_p = len(adj)
    visited = [False] * n_p
    parent = [-1] * n_p
    q: deque[int] = deque()
    q.append(p_start)
    visited[p_start] = True

    found = False
    while q:
        v = q.popleft()
        if v == p_end:
            found = True
            break
        for nb in adj[v]:
            if not visited[nb]:
                visited[nb] = True
                parent[nb] = v
                q.append(nb)

    if not found:
        raise RuntimeError(f"No plaquette path between {p_start} and {p_end}")

    path: list[int] = []
    cur = p_end
    while cur != -1:
        path.append(cur)
        cur = parent[cur]
    path.reverse()
    return path


def branch_cut_bonds_for_pair(Lx: int, Ly: int, pA: int, pB: int) -> list[tuple[int, int]]:
    """Compute a set of bonds forming a Z2 branch cut from plaquette pA to pB.

    We first find a shortest path in plaquette adjacency space, then for
    each neighboring pair of plaquettes pick one of their shared bonds.
    The resulting bond set, when flipped (u_ij -> -u_ij), creates visons
    (W_p = -1) exactly at pA and pB.
    """
    plaquettes = honeycomb_plaquettes(Lx, Ly)
    path = shortest_plaquette_path(Lx, Ly, pA, pB)

    # Precompute bond sets for each plaquette
    plaq_bonds: list[set[tuple[int, int]]] = []
    for cyc in plaquettes:
        edges: set[tuple[int, int]] = set()
        for k in range(6):
            i = cyc[k]
            j = cyc[(k + 1) % 6]
            key = (i, j) if i < j else (j, i)
            edges.add(key)
        plaq_bonds.append(edges)

    cut_bonds: list[tuple[int, int]] = []
    for a, b in zip(path[:-1], path[1:]):
        shared = plaq_bonds[a].intersection(plaq_bonds[b])
        if not shared:
            raise RuntimeError(f"Adjacent plaquettes {a},{b} share no bond")
        # Pick an arbitrary shared bond
        e = sorted(shared)[0]
        cut_bonds.append(e)

    return cut_bonds


def vison_pair_gauge(Lx: int, Ly: int, pA: int, pB: int) -> dict[tuple[int, int], int]:
    """Gauge field with exactly two visons at plaquettes pA and pB.

    Start from the trivial configuration u_ij = +1 and flip u_ij on all
    bonds along a branch cut between pA and pB.
    """
    u = initial_gauge(Lx, Ly)
    for (i, j) in branch_cut_bonds_for_pair(Lx, Ly, pA, pB):
        key = (i, j) if i < j else (j, i)
        u[key] = -u[key]
    return u


def vison_loop_gauge_path(Lx: int, Ly: int) -> list[dict[tuple[int, int], int]]:
    """Construct a closed gauge path: one vison encircles another.

    We use the hexagon centers to impose a 2D grid structure on the
    plaquettes (3x3 for Lx=Ly=4), choose the central plaquette as pA,
    and let pB follow a loop around it through neighboring plaquettes.
    At each step we rebuild the branch cut between pA and the current pB,
    ensuring that W_p = -1 only at those two plaquettes.
    """
    centers = plaquette_centers(Lx, Ly)
    n_p = len(centers)
    if n_p == 0:
        raise ValueError("No plaquettes found")

    # Order plaquettes by (y, x) of their centers to form a grid
    order = sorted(range(n_p), key=lambda p: (centers[p][1], centers[p][0]))

    # Infer grid size (assume square grid for small systems like 3x3)
    n_side = int(round(np.sqrt(float(n_p))))
    if n_side * n_side != n_p:
        raise ValueError(f"Plaquette grid is not square: n_p={n_p}")

    # Map grid coordinates -> plaquette index
    def grid_to_plaq(i_row: int, j_col: int) -> int:
        return order[i_row * n_side + j_col]

    # Choose central plaquette as pA
    c = n_side // 2
    pA = grid_to_plaq(c, c)

    # Define a simple rectangular loop for pB around the center on the grid
    loop_coords = [
        (c, c - 1),  # left
        (c - 1, c - 1),  # up-left
        (c - 1, c),  # up
        (c - 1, c + 1),  # up-right
        (c, c + 1),  # right
        (c + 1, c + 1),  # down-right
        (c + 1, c),  # down
        (c + 1, c - 1),  # down-left
        (c, c - 1),  # back to left to close
    ]

    configs: list[dict[tuple[int, int], int]] = []
    for (i_row, j_col) in loop_coords:
        if not (0 <= i_row < n_side and 0 <= j_col < n_side):
            raise ValueError("Loop coordinates out of range for plaquette grid")
        pB = grid_to_plaq(i_row, j_col)
        u = vison_pair_gauge(Lx, Ly, pA, pB)
        configs.append(u)

    # Ensure the loop is explicitly closed by repeating the first config
    if configs and configs[-1] is not configs[0]:
        configs.append(configs[0])

    return configs


def compute_gauge_loop_holonomy(Lx: int = 4, Ly: int = 4,
                                Jx: float = 1.0, Jy: float = 1.0, Jz: float = 1.0,
                                dim_sub: int = 2,
                                n_steps: int = 20):
    """Compute Berry holonomy for the low-energy subspace along a gauge loop.

    We:
      1. Build a sequence of gauge configurations u^{(k)} forming a closed
         loop in Z2 gauge-field space (toy_gauge_path).
      2. For each u^{(k)}, build the Majorana Hamiltonian H_k.
      3. Extract the lowest-|E| subspace of dimension dim_sub.
      4. Parallel transport this subspace around the loop using overlap
         matrices projected to the closest unitary via polar decomposition.
    """
    # Use a vison-braiding-inspired gauge loop by default. For comparison
    # with the older toy path, one can swap this to toy_gauge_path.
    configs = vison_loop_gauge_path(Lx, Ly)

    # Initial point
    H0 = build_majorana_matrix(Lx, Ly, Jx, Jy, Jz, configs[0])
    w0, V_prev = low_energy_subspace(H0, dim_sub=dim_sub)

    U_holo = np.eye(dim_sub, dtype=complex)

    for k in range(len(configs) - 1):
        Hk = build_majorana_matrix(Lx, Ly, Jx, Jy, Jz, configs[k + 1])
        wk, V_next = low_energy_subspace(Hk, dim_sub=dim_sub)

        W = V_prev.conj().T @ V_next  # dim_sub x dim_sub overlap
        U_step = project_to_unitary(W)
        U_holo = U_step @ U_holo

        V_prev = V_next

    return U_holo, w0


def scan_Jz_F_Dehn(Lx: int = 4, Ly: int = 4,
                   Jx: float = 1.0, Jy: float = 1.0,
                   Jz_values: list[float] | None = None,
                   dim_sub: int = 2):
    """Scan over Jz and compute F_Dehn(Jz) for the vison-loop Berry holonomy.

    For each Jz, we compute the Berry holonomy U_holo(Jz) along the
    vison_loop_gauge_path, normalize it to SU(2), and compare it with the
    Ising Dehn twist matrix to obtain an SU(2) fidelity F_Dehn(Jz).
    """
    if Jz_values is None:
        Jz_values = [0.5, 0.8, 1.0, 1.2, 1.5]

    results = []
    R_ising_su2, Dehn_ising_su2 = ising_R_and_dehn_su2()

    for Jz in Jz_values:
        U_holo, w0 = compute_gauge_loop_holonomy(Lx=Lx, Ly=Ly,
                                                 Jx=Jx, Jy=Jy, Jz=Jz,
                                                 dim_sub=dim_sub)
        evals, _ = np.linalg.eig(U_holo)
        # Use the smallest-|E| eigenvalue at the initial configuration as a
        # proxy for the gap scale for this Jz.
        gap = float(np.min(np.abs(w0))) if w0 is not None else np.nan
        if dim_sub == 2:
            U_su2 = normalize_to_su2(U_holo)
            F_R, tr_R = su2_fidelity(R_ising_su2, U_su2)
            F_Dehn, tr_Dehn = su2_fidelity(Dehn_ising_su2, U_su2)
        else:
            F_R, tr_R = np.nan, np.nan
            F_Dehn, tr_Dehn = np.nan, np.nan

        results.append({
            "Jz": Jz,
            "evals": evals,
            "gap": gap,
            "F_R": F_R,
            "F_Dehn": F_Dehn,
        })

    return results


def scan_Jy_F_Dehn(Lx: int = 4, Ly: int = 4,
                   Jx: float = 1.0, Jz: float = 1.0,
                   Jy_values: list[float] | None = None,
                   dim_sub: int = 2):
    """Scan over Jy and compute F_Dehn(Jy) for the vison-loop Berry holonomy.

    We fix Jx, Jz and vary Jy along a line in coupling space, reusing the
    same vison-loop gauge path. This gives an anisotropy cut complementary
    to the Jz scan.
    """
    if Jy_values is None:
        Jy_values = [0.5, 0.8, 1.0, 1.2, 1.5]

    results = []
    R_ising_su2, Dehn_ising_su2 = ising_R_and_dehn_su2()

    for Jy in Jy_values:
        U_holo, w0 = compute_gauge_loop_holonomy(Lx=Lx, Ly=Ly,
                                                 Jx=Jx, Jy=Jy, Jz=Jz,
                                                 dim_sub=dim_sub)
        evals, _ = np.linalg.eig(U_holo)
        gap = float(np.min(np.abs(w0))) if w0 is not None else np.nan
        if dim_sub == 2:
            U_su2 = normalize_to_su2(U_holo)
            F_R, tr_R = su2_fidelity(R_ising_su2, U_su2)
            F_Dehn, tr_Dehn = su2_fidelity(Dehn_ising_su2, U_su2)
        else:
            F_R, tr_R = np.nan, np.nan
            F_Dehn, tr_Dehn = np.nan, np.nan

        results.append({
            "Jy": Jy,
            "evals": evals,
            "gap": gap,
            "F_R": F_R,
            "F_Dehn": F_Dehn,
        })

    return results


def scan_triangle_F_Dehn(Lx: int = 4, Ly: int = 4,
                         t_values: list[float] | None = None,
                         dim_sub: int = 2):
    """Scan along a "triangle" path Jx=Jy=(1-t)/2, Jz=t with Jx+Jy+Jz=1.

    This path is motivated by the standard Kitaev honeycomb phase diagram
    where couplings lie on the simplex Jx+Jy+Jz=1 and gapless lines are
    given (schematically) by J_alpha = J_beta + J_gamma. As t increases
    from 0 to 1, the system interpolates between an xy-dominated regime
    and a z-dominated regime, crossing a putative gapless line near
    t ~ 0.5 where Jz ~ Jx+Jy.
    """
    if t_values is None:
        t_values = [0.1 * k for k in range(1, 10)]  # 0.1,0.2,...,0.9

    results = []
    R_ising_su2, Dehn_ising_su2 = ising_R_and_dehn_su2()

    for t in t_values:
        Jz = float(t)
        Jx = Jy = 0.5 * (1.0 - Jz)
        U_holo, w0 = compute_gauge_loop_holonomy(Lx=Lx, Ly=Ly,
                                                 Jx=Jx, Jy=Jy, Jz=Jz,
                                                 dim_sub=dim_sub)
        evals, _ = np.linalg.eig(U_holo)
        gap = float(np.min(np.abs(w0))) if w0 is not None else np.nan
        if dim_sub == 2:
            U_su2 = normalize_to_su2(U_holo)
            F_R, tr_R = su2_fidelity(R_ising_su2, U_su2)
            F_Dehn, tr_Dehn = su2_fidelity(Dehn_ising_su2, U_su2)
        else:
            F_R, tr_R = np.nan, np.nan
            F_Dehn, tr_Dehn = np.nan, np.nan

        results.append({
            "t": t,
            "Jx": Jx,
            "Jy": Jy,
            "Jz": Jz,
            "evals": evals,
            "gap": gap,
            "F_R": F_R,
            "F_Dehn": F_Dehn,
        })

    return results


def plot_F_Dehn_vs_Jz(results: list[dict], out_path: str = "honeycomb_vison_F_Dehn_vs_Jz.png"):
    """Plot F_Dehn(Jz) (and optionally F_R) for the vison loop and save as PNG."""
    Jz_vals = np.array([r["Jz"] for r in results], dtype=float)
    F_R_vals = np.array([r["F_R"] for r in results], dtype=float)
    F_Dehn_vals = np.array([r["F_Dehn"] for r in results], dtype=float)

    plt.figure(figsize=(5, 3.5))
    plt.plot(Jz_vals, F_Dehn_vals, marker="o", label=r"$F_{\mathrm{Dehn}}(J_z)$")
    plt.plot(Jz_vals, F_R_vals, marker="s", linestyle="--", label=r"$F_R(J_z)$")
    plt.xlabel(r"$J_z$")
    plt.ylabel(r"Fidelity")
    plt.ylim(0.0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved honeycomb vison-loop F vs Jz plot to {out_path}")


def plot_gap_vs_Jz(results: list[dict], out_path: str = "honeycomb_vison_gap_vs_Jz.png"):
    """Plot the smallest-|E| scale (gap proxy) vs Jz for the vison configuration."""
    Jz_vals = np.array([r["Jz"] for r in results], dtype=float)
    gaps = np.array([r["gap"] for r in results], dtype=float)

    plt.figure(figsize=(5, 3.5))
    plt.plot(Jz_vals, gaps, marker="o")
    plt.xlabel(r"$J_z$")
    plt.ylabel(r"$\min |E|(J_z)$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved honeycomb vison-loop gap vs Jz plot to {out_path}")


def plot_F_Dehn_vs_Jy(results: list[dict], out_path: str = "honeycomb_vison_F_Dehn_vs_Jy.png"):
    """Plot F_Dehn(Jy) and F_R(Jy) for the vison loop and save as PNG."""
    Jy_vals = np.array([r["Jy"] for r in results], dtype=float)
    F_R_vals = np.array([r["F_R"] for r in results], dtype=float)
    F_Dehn_vals = np.array([r["F_Dehn"] for r in results], dtype=float)

    plt.figure(figsize=(5, 3.5))
    plt.plot(Jy_vals, F_Dehn_vals, marker="o", label=r"$F_{\mathrm{Dehn}}(J_y)$")
    plt.plot(Jy_vals, F_R_vals, marker="s", linestyle="--", label=r"$F_R(J_y)$")
    plt.xlabel(r"$J_y$")
    plt.ylabel(r"Fidelity")
    plt.ylim(0.0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved honeycomb vison-loop F vs Jy plot to {out_path}")


def plot_gap_vs_Jy(results: list[dict], out_path: str = "honeycomb_vison_gap_vs_Jy.png"):
    """Plot the smallest-|E| scale (gap proxy) vs Jy for the vison configuration."""
    Jy_vals = np.array([r["Jy"] for r in results], dtype=float)
    gaps = np.array([r["gap"] for r in results], dtype=float)

    plt.figure(figsize=(5, 3.5))
    plt.plot(Jy_vals, gaps, marker="o")
    plt.xlabel(r"$J_y$")
    plt.ylabel(r"$\min |E|(J_y)$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved honeycomb vison-loop gap vs Jy plot to {out_path}")


def plot_F_Dehn_vs_triangle_t(results: list[dict], out_path: str = "honeycomb_vison_F_Dehn_vs_triangle_t.png"):
    """Plot F_Dehn(t) and F_R(t) along the triangle path and save as PNG."""
    t_vals = np.array([r["t"] for r in results], dtype=float)
    F_R_vals = np.array([r["F_R"] for r in results], dtype=float)
    F_Dehn_vals = np.array([r["F_Dehn"] for r in results], dtype=float)

    plt.figure(figsize=(5, 3.5))
    plt.plot(t_vals, F_Dehn_vals, marker="o", label=r"$F_{\mathrm{Dehn}}(t)$")
    plt.plot(t_vals, F_R_vals, marker="s", linestyle="--", label=r"$F_R(t)$")
    plt.xlabel(r"$t$")
    plt.ylabel(r"Fidelity")
    plt.ylim(0.0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved honeycomb vison-loop F vs triangle t plot to {out_path}")


def plot_gap_vs_triangle_t(results: list[dict], out_path: str = "honeycomb_vison_gap_vs_triangle_t.png"):
    """Plot the smallest-|E| scale (gap proxy) vs t along the triangle path."""
    t_vals = np.array([r["t"] for r in results], dtype=float)
    gaps = np.array([r["gap"] for r in results], dtype=float)

    plt.figure(figsize=(5, 3.5))
    plt.plot(t_vals, gaps, marker="o")
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\min |E|(t)$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved honeycomb vison-loop gap vs triangle t plot to {out_path}")


def main():
    # Lattice and couplings: small system in a generic anisotropic regime.
    Lx, Ly = 4, 4
    Jx, Jy, Jz = 1.0, 1.0, 1.0
    dim_sub = 2
    n_steps = 20

    print("=== Honeycomb Majorana toy model: initial spectrum ===")
    u0 = initial_gauge(Lx, Ly)
    H0 = build_majorana_matrix(Lx, Ly, Jx, Jy, Jz, u0)
    w, _ = np.linalg.eigh(H0)
    idx = np.argsort(np.abs(w))
    sel = idx[:8]
    print("Lowest-|E| eigenvalues (first 8):")
    print(np.round(w[sel], 6))

    # Optional diagnostic: compute Z2 fluxes W_p on hexagonal plaquettes
    fluxes0 = plaquette_fluxes(Lx, Ly, u0)
    fluxes0 = np.array(fluxes0, dtype=int)
    print("Number of plaquettes:", len(fluxes0))
    if len(fluxes0) > 0:
        print("Initial W_p (first up to 8):", fluxes0[:8])

    print("\n=== Berry holonomy along a vison-loop gauge path ===")
    U_holo, w0 = compute_gauge_loop_holonomy(Lx=Lx, Ly=Ly,
                                             Jx=Jx, Jy=Jy, Jz=Jz,
                                             dim_sub=dim_sub,
                                             n_steps=n_steps)
    evals, _ = np.linalg.eig(U_holo)
    print("U_holo ({}x{}) =".format(dim_sub, dim_sub))
    print(np.round(U_holo, 6))
    print("Eigenvalues of U_holo:")
    print(np.round(np.sort_complex(evals), 6))

    # Optional: compare to Ising Dehn twist in SU(2)
    if dim_sub == 2:
        U_holo_su2 = normalize_to_su2(U_holo)
        R_ising_su2, Dehn_ising_su2 = ising_R_and_dehn_su2()
        F_R, tr_R = su2_fidelity(R_ising_su2, U_holo_su2)
        F_Dehn, tr_Dehn = su2_fidelity(Dehn_ising_su2, U_holo_su2)

        print("\n=== Phase-agnostic comparison with Ising SU(2) data ===")
        print("Tr(R_ising_su2^† U_holo_su2) =", tr_R)
        print("F_R    = |Tr|/2 =", F_R)
        print("Tr(Dehn_ising_su2^† U_holo_su2) =", tr_Dehn)
        print("F_Dehn = |Tr|/2 =", F_Dehn)

    # Small Jz scan around the isotropic point to probe F_Dehn(Jz)
    print("\n=== Small Jz scan: F_Dehn(Jz) along vison loop ===")
    Jz_list = [0.5, 0.8, 1.0, 1.2, 1.5]
    results = scan_Jz_F_Dehn(Lx=Lx, Ly=Ly, Jx=Jx, Jy=Jy,
                             Jz_values=Jz_list, dim_sub=dim_sub)
    print("Jz\tF_R\tF_Dehn")
    for r in results:
        print(f"{r['Jz']:.3f}\t{r['F_R']:.6f}\t{r['F_Dehn']:.6f}")

    # Generate a simple PNG plot summarizing the scan
    plot_F_Dehn_vs_Jz(results)
    plot_gap_vs_Jz(results)

    # Anisotropy scan in Jy with Jx=Jz=1
    print("\n=== Small Jy scan: F_Dehn(Jy) along vison loop ===")
    Jy_list = [0.5, 0.8, 1.0, 1.2, 1.5]
    results_Jy = scan_Jy_F_Dehn(Lx=Lx, Ly=Ly, Jx=Jx, Jz=Jz,
                                Jy_values=Jy_list, dim_sub=dim_sub)
    print("Jy\tF_R\tF_Dehn")
    for r in results_Jy:
        print(f"{r['Jy']:.3f}\t{r['F_R']:.6f}\t{r['F_Dehn']:.6f}")

    plot_F_Dehn_vs_Jy(results_Jy)
    plot_gap_vs_Jy(results_Jy)


if __name__ == "__main__":
    main()
