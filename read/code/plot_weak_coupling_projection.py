"""Illustrate weak-coupling projection onto the moving PRB111 low-energy frame.

All curves evaluate explicit instantaneous coefficients, for eta/E1 = 1.
No Schrodinger evolution, leakage probability, or fidelity is calculated.
"""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss_frame_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from plot_parallel_transport_frame import (
    configure_style, geometry, INK, MUTED, GRID, FRAME1, FRAME2,
)


OUT = Path(__file__).resolve().parents[1] / "figures" / "parallel_transport"
OFF = "#C95653"
LOW = "#078777"
STRONG = "#2976AC"


def weak_components(stage, s, ratio=1.):
    """Coefficients in units of E1, in the continuous (f1,f2,n) frame."""
    controls, n, f1, f2 = geometry(stage, s)
    w = np.column_stack((ratio*controls[:, 1], np.ones(len(n)), np.zeros(len(n))))
    omega = np.einsum("ij,ij->i", w, n)
    w1 = np.einsum("ij,ij->i", w, f1)
    w2 = np.einsum("ij,ij->i", w, f2)
    return omega, w1, w2


def arrow3(ax, start, vector, color, lw=2.7):
    ax.quiver(*start, *vector, color=color, linewidth=lw,
              arrow_length_ratio=.14, normalize=False)


def draw_projection(ax):
    ax.set_axis_off()
    ax.set_proj_type("ortho")
    ax.view_init(elev=25, azim=30)
    ax.set_box_aspect((1.15, 1.15, 1.2))
    ax.set(xlim=(-1.05, 1.1), ylim=(-.8, 1.2), zlim=(-.25, 1.3))
    corners = np.array([[-.95, -.65, 0], [1., -.65, 0],
                        [1., 1.05, 0], [-.95, 1.05, 0]])
    ax.add_collection3d(Poly3DCollection([corners], facecolor="#DCEEEB",
                                        edgecolor="#A4CCC6", alpha=.43,
                                        linewidth=.8))
    for vector, color, label in (([1.05, 0, 0], FRAME1, r"$\mathbf{f}_1$"),
                                  ([0, 1.15, 0], FRAME2, r"$\mathbf{f}_2$"),
                                  ([0, 0, 1.23], MUTED, r"$\mathbf{n}$")):
        arrow3(ax, np.zeros(3), vector, color, lw=1.3)
        ax.text(*(np.asarray(vector)*1.08), label, color=color,
                fontsize=15, ha="center")
    omega, w1, w2 = weak_components(1, [.5])
    w = np.array([w1[0], w2[0], omega[0]])
    perp = np.array([w1[0], w2[0], 0.])
    parallel = np.array([0., 0., omega[0]])
    arrow3(ax, np.zeros(3), w, INK, lw=3.)
    arrow3(ax, np.zeros(3), perp, LOW, lw=3.1)
    arrow3(ax, np.zeros(3), parallel, OFF, lw=3.1)
    ax.plot(*np.array([perp, w]).T, ls="--", color=OFF, lw=1.5)
    ax.plot(*np.array([parallel, w]).T, ls="--", color=LOW, lw=1.5)
    ax.scatter(0, 0, 0, s=24, color=INK, depthshade=False)
    ax.text(*(w + np.array([-.10, .04, .14])), r"$\mathbf{w}$", fontsize=19, color=INK)
    ax.text(*(perp + np.array([-.20, -.08, -.14])), r"$\mathbf{w}_{\perp}$",
            fontsize=17, color=LOW)
    ax.text(.10, -.17, .48, r"$\omega\mathbf{n}$", fontsize=16, color=OFF)
    ax.text2D(.045, .085, r"$\mathbf{w}=\omega\mathbf{n}+w_1\mathbf{f}_1+w_2\mathbf{f}_2$",
              transform=ax.transAxes, fontsize=15)
    ax.text2D(.045, .015, "绿色平面：剩余低能方向；法向：辅助模式方向",
              transform=ax.transAxes, fontsize=10.5, color=MUTED)


def node(ax, x, y, label, color):
    ax.add_patch(Circle((x, y), .061, fc="white", ec=color, lw=2.1, zorder=4))
    ax.text(x, y, label, ha="center", va="center", color=color, fontsize=18, zorder=5)


def draw_couplings(ax):
    ax.set_axis_off()
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.add_patch(FancyBboxPatch((.05, .70), .90, .24,
                                boxstyle="round,pad=.015,rounding_size=.025",
                                fc="#EFF5FA", ec="#B9D1E2", lw=1.))
    ax.add_patch(FancyBboxPatch((.05, .115), .90, .42,
                                boxstyle="round,pad=.015,rounding_size=.025",
                                fc="#F0F7F5", ec="#BAD8D1", lw=1.))
    ax.text(.5, .90, "辅助费米子模式：激发能隙 2R", ha="center",
            fontsize=11.5, color=STRONG)
    ax.plot([.23, .77], [.79, .79], color=STRONG, lw=5.)
    ax.text(.5, .805, r"$R$", ha="center", va="bottom", fontsize=17, color=STRONG)
    ax.plot([.50, .77], [.46, .79], color=OFF, lw=2.7, ls="--")
    ax.text(.78, .62, r"$\omega$", color=OFF, ha="center", fontsize=18)
    ax.text(.12, .61, "可连接高、低能空间", color=OFF, fontsize=10.8)
    ax.plot([.50, .22], [.46, .26], color=FRAME1, lw=3.)
    ax.plot([.50, .78], [.46, .26], color=FRAME2, lw=3.)
    ax.text(.28, .40, r"$w_1$", color=FRAME1, fontsize=17)
    ax.text(.69, .40, r"$w_2$", color=FRAME2, fontsize=17)
    node(ax, .23, .79, r"$\gamma_a$", STRONG)
    node(ax, .77, .79, r"$\gamma_c$", STRONG)
    node(ax, .50, .46, r"$\gamma_1$", INK)
    node(ax, .22, .26, r"$\widetilde{\gamma}_2$", FRAME1)
    node(ax, .78, .26, r"$\widetilde{\gamma}_3$", FRAME2)
    ax.text(.5, .16, "这两条耦合直接作用于低能内部", color=LOW,
            ha="center", fontsize=11.)
    ax.text(.5, .045, r"$h_K=-w_1 Z+w_2 Y$", ha="center", fontsize=20)


def draw_coefficients(ax):
    s = np.linspace(0, 1, 301)
    for stage, fill in enumerate(("#F2F6FB", "#EDF7F4", "#FCF3EF")):
        ax.axvspan(stage, stage+1, color=fill, zorder=0)
        omega, w1, w2 = weak_components(stage, s)
        t = stage+s
        ax.plot(t, -w1, color=FRAME1, lw=2.8, zorder=4)
        ax.plot(t, w2, color=FRAME2, lw=2.8, zorder=4)
        ax.plot(t, omega, color=OFF, ls="--", lw=2.1, zorder=3)
        ax.text(stage+.5, 1.23, ("I", "II", "III")[stage],
                ha="center", fontsize=14, color=MUTED)
    ax.axhline(0, color=GRID, lw=.9, zorder=1)
    ax.set(xlim=(0, 3), ylim=(-1.12, 1.38),
           xlabel=r"一圈中的时间 $t/\tau$", ylabel=r"耦合系数 / $E_1$")
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([-1, 0, 1])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=11)
    handles = [
        Line2D([0], [0], color=FRAME1, lw=2.8, label=r"$Z$ 系数 $-w_1$"),
        Line2D([0], [0], color=FRAME2, lw=2.8, label=r"$Y$ 系数 $w_2$"),
        Line2D([0], [0], color=OFF, ls="--", lw=2.1, label=r"跨能区耦合 $\omega$"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, 1.01),
               frameon=False, fontsize=10.5, ncol=1, handlelength=2.8,
               borderaxespad=0, labelspacing=.3)


def main():
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(18, 9.3), facecolor="white")
    fig.text(.04, .951, "弱耦合留下什么？从模式分解到低能动态演化",
             fontsize=23, fontweight="bold")
    fig.text(.04, .908, "同一个弱耦合既可连接辅助模式，也可直接改变逻辑态；低能近似保留后者的领先作用",
             fontsize=12.5, color=MUTED)
    for x, title, sub in (
        (.035, "A  将弱耦合向量分解", "示例：第二段中点；η = E₁；系数以 E₁ 为单位"),
        (.365, "B  这些分量分别耦合谁？", "连线表示耦合，不表示已经发生了跃迁"),
        (.704, "C  低能动态作用随时间变化", "原三段波形，η / E₁ = 1；直接绘制解析系数"),
    ):
        fig.text(x, .845, title, fontsize=15, fontweight="bold")
        fig.text(x, .807, sub, fontsize=10.5, color=MUTED)
    projection = fig.add_axes([.016, .325, .325, .46], projection="3d")
    couplings = fig.add_axes([.365, .32, .29, .46])
    coefficients = fig.add_axes([.745, .385, .225, .30])
    draw_projection(projection)
    draw_couplings(couplings)
    draw_coefficients(coefficients)
    fig.text(.735, .325, "紫色：相对相位；金色：逻辑态混合\n红色虚线是耦合强度，不是泄漏概率",
             fontsize=11, color=MUTED, linespacing=1.6, va="top")
    fig.add_artist(plt.Line2D([.04, .96], [.274, .274], transform=fig.transFigure,
                              color="#DCE4EA", lw=1))
    fig.text(.045, .228, r"$\mathbf{w}=(u_1,E_1,0),\quad \omega=\mathbf{w}\cdot\mathbf{n},\quad w_i=\mathbf{w}\cdot\mathbf{f}_i$",
             fontsize=16)
    fig.text(.54, .228, r"$V=i\gamma_1(\omega\gamma_c+w_1\widetilde{\gamma}_2+w_2\widetilde{\gamma}_3)$",
             fontsize=16)
    fig.text(.045, .172, "领先低能投影保留：", fontsize=12, color=LOW)
    fig.text(.205, .172, r"$h_{\mathrm{I}}=\frac{\eta a^2-E_1b}{q}Z\qquad h_{\mathrm{II}}=\eta bZ+E_1gY\qquad h_{\mathrm{III}}=E_1Y$",
             fontsize=19)
    fig.text(.045, .105, "Z、Y 是随动逻辑基底中的 Pauli 轴。第二段的旋转轴随时间变化，必须保留时间有序性。",
             fontsize=12, color=INK)
    fig.text(.045, .059, "跨能区耦合并未从完整模型中消失：它可带来高低能混合、虚跃迁修正和适当条件下的真实激发。",
             fontsize=11, color=MUTED)
    fig.text(.045, .026, "本图不计算泄漏概率或保真度；由强控制基底运动引起的有限速度激发还需另行考虑。",
             fontsize=11, color=MUTED)
    for suffix in ("png", "pdf", "svg"):
        path = OUT / f"weak_coupling_projection.{suffix}"
        fig.savefig(path, dpi=190, facecolor="white")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
