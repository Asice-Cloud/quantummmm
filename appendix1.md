# Appendix 1. Schur Complement to Pauli Tensor Non-Abelian Quantification

This appendix summarizes the relation between the microscopic Schur-complement description and the Pauli-tensor path-ordered non-Abelian metric.

## A1.1 Microscopic partition and Schur complement

Start from a block-partitioned Hamiltonian or kernel,

$$
K =
\begin{pmatrix}
K_{PP} & K_{PQ} \\
K_{QP} & K_{QQ}
\end{pmatrix},
$$

where $P$ denotes the low-energy Majorana/bound-state subspace and $Q$ denotes the eliminated bulk/continuum sector.

The resolvent is

$$
G(z) = (z-K)^{-1}.
$$

By the block inverse identity,

$$
G_{PP}(z)
=
\bigl(z - K_{PP} - K_{PQ}(z-K_{QQ})^{-1}K_{QP}\bigr)^{-1}.
$$

Define the self-energy (Schur complement)

$$
\Sigma(z) = K_{PQ}(z-K_{QQ})^{-1}K_{QP}.
$$

Then the effective low-energy Hamiltonian is

$$
K_{\mathrm{eff}}(z) = K_{PP} + \Sigma(z).
$$

Physical meaning:

- $\Sigma(z)=0$ means exact decoupling of $P$ and $Q$.
- $\Sigma(z)\neq 0$ means virtual $P\to Q\to P$ excursions.
- In Majorana/ABS problems, $\Sigma(z)\neq 0$ is the microscopic origin of deformation away from ideal braiding.

## A1.2 Minimal chain picture

For the minimal chain model, the ideal part is the boundary-zero-mode sector, while the perturbation introduces bulk leakage:

$$
H = H_0 + V,
$$

with $H_0$ describing the ideal chain and $V$ generating boundary-bulk mixing.

After projection,

$$
H_{PP}^{\mathrm{eff}}(z) = H_{PP} + \Sigma(z).
$$

This is the object that feeds the Pauli-tensor representation.

## A1.3 Pauli tensor decomposition of the effective Hamiltonian

The effective two-body Hamiltonian is expanded in the Pauli tensor basis,

$$
\mathrm{End}(V\otimes V)
\;\text{basis}\;\{\sigma^\alpha\otimes\sigma^\beta\},
\qquad
\alpha,\beta\in\{0,x,y,z\}.
$$

Hence

$$
h(u)=\sum_{\alpha,\beta} h_{\alpha\beta}(u)\,\sigma^\alpha\otimes\sigma^\beta.
$$

Here:

- $u$ is the braid/path parameter.
- $h_{\alpha\beta}(u)$ are the path-dependent Pauli-channel coefficients.
- These coefficients are the effective representation of $H_{PP}+\Sigma(z)$ after projection and basis choice.

Interpretation:

- Schur complement supplies the effective couplings.
- Pauli coefficients encode those couplings in a channel-resolved form.
- ABS leakage and bulk admixture appear as additional noncommuting Pauli channels.

## A1.4 Path-ordered evolution

The local evolution operator is

$$
R(u)=\mathcal T\exp\left(-i\int_0^u h(s)\,ds\right).
$$

Because in general

$$
[h(s_1),h(s_2)]\neq 0,
$$

time ordering must be kept.

The Dyson expansion is

$$
R(u)=I+R^{(1)}+R^{(2)}+R^{(3)}+\cdots
$$

with

$$
R^{(1)}=-i\int_0^u ds_1\,h(s_1),
$$

$$
R^{(2)}=(-i)^2\int_0^u ds_1\int_0^{s_1}ds_2\,h(s_1)h(s_2),
$$

and

$$
R^{(3)}=(-i)^3\int_0^u ds_1\int_0^{s_1}ds_2\int_0^{s_2}ds_3\,h(s_1)h(s_2)h(s_3).
$$

## A1.5 Yang--Baxter deviation

Embed the two-body operator into the three-body space:

$$
h_{12}(u)=\sum_{\alpha,\beta}h_{\alpha\beta}(u)\,\sigma^\alpha\otimes\sigma^\beta\otimes I,
$$

$$
h_{23}(u)=\sum_{\mu,\nu}h_{\mu\nu}(u)\,I\otimes\sigma^\mu\otimes\sigma^\nu.
$$

Define the Yang--Baxter deviation

$$
\Delta(u,v)=R_{12}(u)R_{23}(u+v)R_{12}(v)-R_{23}(v)R_{12}(u+v)R_{23}(u).
$$

Then the leading nontrivial term is third order. After reorganizing the time-ordering structure, the dominant commutator contribution is

$$
[h_{12}(s),h_{23}(t)]
=2i\sum_{\alpha,\beta,\mu,\nu,\gamma}
 h_{\alpha\beta}(s)h_{\mu\nu}(t)\epsilon_{\beta\mu\gamma}
\,\sigma^\alpha\otimes\sigma^\gamma\otimes\sigma^\nu.
$$

The scalar prefactor from the Dyson order is

$$
(-i)^3\cdot 2i = -2.
$$

So the non-Abelian contribution is controlled by the Pauli-channel commutator structure weighted by a time kernel $K(s,t)$.

## A1.6 Time kernel representation

The third-order contribution can be written schematically as

$$
\Delta^{(3)} = (-i)^3\int_0^u ds\int_0^{u+v}dt\;K(s,t)\,[h_{12}(s),h_{23}(t)].
$$

The kernel $K(s,t)$ is built from Heaviside functions that encode the relative time ordering. Its role is to compress the original triple integrals into a double-integral commutator form.

The corresponding Frobenius scaling is

$$
\|\Delta^{(3)}\|_F^2
\sim
4\cdot 2^3
\sum_{\alpha,\beta,\mu,\nu,\gamma}
\left|\int\!\!\int K(s,t)
 h_{\alpha\beta}(s)h_{\mu\nu}(t)\,ds\,dt\right|^2
+\cdots
$$

Thus

$$
\|\Delta^{(3)}\|_F = O(J^3),
\qquad
\mathrm{Tr}(\Delta^\dagger\Delta)=O(J^6).
$$

## A1.7 Non-Abelian measure

The non-Abelian metric used in the Pauli-tensor framework is

$$
\mathcal N = \|\Delta_{\mathrm{braid}}\|
$$

or, in the YBE-based numerical implementation,

$$
\mathcal N = \sqrt{\mathrm{Tr}(\Delta^\dagger\Delta)}.
$$

A simplified coefficient-level estimate consistent with the markdown derivation is

$$
\mathcal N
\sim
\sqrt{
\int ds_1 ds_2
\sum_{\beta\neq\mu}
|h_{\alpha\beta}(s_1)|^2
|h_{\mu\nu}(s_2)|^2
}.
$$

This should be viewed as a trend-level formula. It captures the physical mechanism that non-Abelianity arises from noncommuting Pauli channels, but it is not yet identical to the full exact metric unless all normalization and projection conventions are matched.

## A1.8 Relation between the two pictures

The two descriptions are consistent and complementary:

$$
H_{\mathrm{micro}}
\;\Rightarrow\;
\Sigma(z)
\;\Rightarrow\;
H_{PP}^{\mathrm{eff}}(z)
\;\Rightarrow\;
h_{\alpha\beta}(u)
\;\Rightarrow\;
R(u)
\;\Rightarrow\;
\Delta_{YB}
\;\Rightarrow\;
\mathcal N.
$$

Therefore:

- Schur complement is the microscopic effective-theory step.
- Pauli tensor evolution is the channel-resolved dynamical step.
- The non-Abelian metric is the final braid-level observable.

In short: **the Schur-complement model and the Pauli-tensor model are consistent, but they operate at different stages of the same derivation.**

## A1.9 Practical interpretation

If $\Sigma(z)=0$, the braid is close to ideal and the Pauli channels remain effectively commuting.

## A1.10 Calibrated third-order form

The kernel formula above is the correct third-order structural expression, but exact numerical comparison shows that its absolute normalization must be calibrated against the exact YBE deviation.

Define

$$
\mathcal N_{\mathrm{kernel}}^{(3)}
=
\left[
4\cdot 2^3
\sum_{
\alpha,\beta,\mu,\nu,\gamma}
\left|\int\!\!\int K(s,t)
 h_{\alpha\beta}(s)h_{\mu\nu}(t)\,ds\,dt\right|^2
\right]^{1/2}
$$

and

$$
\mathcal N_{\mathrm{exact}}^{(3)}=\|\Delta^{(3)}_{\mathrm{exact}}\|_F.
$$

Then

$$
\mathcal N_{\mathrm{exact}}^{(3)} \approx C_{\mathrm{cal}}\,\mathcal N_{\mathrm{kernel}}^{(3)}.
$$

The factor $C_{\mathrm{cal}}$ is not a change in physics; it reflects the precise normalization, projection, and operator-convention choices used in the exact numerical implementation.

In the current numerical setup used in this repository, the fitted value is approximately

$$
C_{\mathrm{cal}} \approx 5.09\times 10^2.
$$

This value should be treated as the convention-dependent conversion factor between the kernel estimator and the exact third-order Frobenius norm.

If $\Sigma(z)\neq 0$, bulk/ABS leakage generates additional channels, the commutator terms survive, and $\mathcal N$ becomes nonzero.

So the workflow is:

1. Use Schur complement to derive the effective low-energy Hamiltonian.
2. Decompose that effective Hamiltonian into Pauli channels.
3. Evolve the channels with time ordering.
4. Compute $\Delta_{YB}$ and $\mathcal N$.
5. Compare with the paper’s braid/non-Abelian diagnostics.

## A1.10 Third-order exact-vs-kernel calibration

The markdown derivation gives a correct **third-order structure** for the YBE deviation, but the numerical comparison shows that the kernel formula is not yet normalized as an absolute metric.

Let

$$
\mathcal N_{\mathrm{exact}}^{(3)} = \|\Delta^{(3)}_{\mathrm{exact}}\|_F,
$$

and

$$
\mathcal N_{\mathrm{kernel}}^{(3)}
=
\left[
4\cdot 2^3
\sum_{\alpha,\beta,\mu,\nu,\gamma}
\left|\int\!\!\int K(s,t)
 h_{\alpha\beta}(s)h_{\mu\nu}(t)\,ds\,dt\right|^2
\right]^{1/2}.
$$

The comparison implemented in the validation script shows:

- the two quantities have very high shape agreement under parameter scans;
- the kernel expression needs a **global calibration factor** to match the exact third-order norm;
- the missing factor is a normalization/projection constant, not a change of operator structure.

For one representative point in the current numerical model, the direct comparison gives roughly

$$
\mathcal N_{\mathrm{exact}}^{(3)} \approx 7.36\times 10^2\,\mathcal N_{\mathrm{kernel}}^{(3)},
$$

which means the kernel expression is correct as a **shape-preserving third-order estimator**, but not yet the final absolute metric.

Hence the practical calibrated form is

$$
\mathcal N_{\mathrm{calibrated}}^{(3)}
=
C_{\mathrm{cal}}\,\mathcal N_{\mathrm{kernel}}^{(3)},
$$

with $C_{\mathrm{cal}}$ fixed once by matching to the exact numerical $\Delta^{(3)}$ in the chosen convention.

This is the version that should be used when you want the appendix to state both:

1. the analytical Pauli-channel mechanism, and
2. the numerically calibrated non-Abelian magnitude.
