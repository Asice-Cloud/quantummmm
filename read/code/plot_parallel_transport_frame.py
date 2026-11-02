"""Draw the PRB111 reference path and its parallel-transported Majorana frame.

This evaluates the explicit frame formulas; it does not solve dynamics or
validate fidelity. The arrows describe Majorana coefficient vectors, not
many-body state vectors. Their spin representation transports the two columns
of the logical basis matrix Q. Export PNG, PDF and SVG for the report.
"""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss_frame_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


OUT = Path(__file__).resolve().parents[1] / "figures" / "parallel_transport"
INK = "#19334A"
MUTED = "#607487"
GRID = "#BCCBD4"
STAGES = ("#2374BC", "#079280", "#E16B45")
FRAME1 = "#7954AB"
FRAME2 = "#C08A18"


def configure_style():
    font_path = Path("/usr/share/fonts/adobe-source-han-sans/SourceHanSans.ttc")
    if font_path.exists():
        # Matplotlib otherwise uses the collection's first, ExtraLight face.
        from fontTools.ttLib import TTCollection
        font_cache = Path(os.environ["MPLCONFIGDIR"]) / "frame_fonts"
        font_cache.mkdir(parents=True, exist_ok=True)
        font_files = [font_cache / f"SourceHanSansSC-{w}.otf"
                      for w in ("Regular", "Bold")]
        if not all(p.exists() for p in font_files):
            collection = TTCollection(font_path, lazy=True)
            for font in collection.fonts:
                names = font["name"]
                family = names.getBestFamilyName()
                style = names.getBestSubFamilyName()
                if family == "Source Han Sans SC" and style in ("Regular", "Bold"):
                    font.save(font_cache / f"SourceHanSansSC-{style}.otf")
            collection.close()
        for path in font_files:
            font_manager.fontManager.addfont(str(path))
        chinese_font = "Source Han Sans SC"
    else:
        chinese_font = "sans-serif"
    plt.rcParams.update({
        "font.family": [chinese_font, "DejaVu Sans"],
        "font.size": 12,
        "text.color": INK,
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",
    })


def geometry(stage, s):
    """Return physical controls (d,u2,u3)/C and n,f1,f2 in (b,2,3)."""
    s = np.atleast_1d(s)
    a = (1 - np.cos(np.pi * s)) / 2
    b = 1 - a
    q = np.hypot(a, b)
    z = np.zeros_like(a)
    o = np.ones_like(a)
    if stage == 0:
        controls = np.column_stack((b, a, z))
        n = np.column_stack((b / q, a / q, z))
        f1 = np.column_stack((-a / q, b / q, z))
        f2 = np.column_stack((z, z, o))
    elif stage == 1:
        controls = np.column_stack((z, b, a))
        n = np.column_stack((z, b / q, -a / q))
        f1 = np.column_stack((-o, z, z))
        f2 = np.column_stack((z, a / q, b / q))
    elif stage == 2:
        controls = np.column_stack((a, z, b))
        n = np.column_stack((a / q, z, -b / q))
        f1 = np.column_stack((-b / q, z, -a / q))
        f2 = np.column_stack((z, o, z))
    else:
        raise ValueError("stage must be 0, 1 or 2")
    return controls, n, f1, f2


def line3(ax, points, **kwargs):
    ax.plot(points[:, 0], points[:, 1], points[:, 2], **kwargs)


def arrow3(ax, origin, vector, color, length=1., width=2., alpha=1.):
    vector = np.asarray(vector) * length
    ax.quiver(*origin, *vector, color=color, linewidth=width,
              arrow_length_ratio=.24, normalize=False, alpha=alpha)


def path_arrow(ax, stage, field, color):
    points = geometry(stage, np.array([.46, .65]))[field]
    ax.quiver(*points[0], *(points[1] - points[0]), color=color,
              linewidth=2.5, arrow_length_ratio=.52, normalize=False)


def panel_heading(ax, number, title, subtitle):
    ax.text2D(.01, 1.01, f"{number}  {title}", transform=ax.transAxes,
              fontsize=15, fontweight="bold")
    ax.text2D(.01, .955, subtitle, transform=ax.transAxes,
              fontsize=11, color=MUTED)


def draw_controls(ax):
    panel_heading(ax, "A", "门控参数走过闭环", "坐标是耦合强度，不是器件中的空间位置")
    ax.set_axis_off()
    ax.set_proj_type("ortho")
    ax.view_init(elev=23, azim=37)
    ax.set_box_aspect((1, 1, 1))
    ax.set(xlim=(-.13, 1.2), ylim=(-.13, 1.2), zlim=(-.13, 1.2))
    origin = np.zeros(3)
    for axis in np.eye(3):
        arrow3(ax, origin, axis, GRID, length=1.21, width=1.2)
    for xyz, label in (([1.25, 0, 0], r"$d/C$"),
                       ([0, 1.25, 0], r"$u_2/C$"),
                       ([.15, 0, 1.16], r"$u_3/C$")):
        ax.text(*xyz, label, fontsize=15, ha="center", color=MUTED)
    vertices = np.eye(3)
    ax.add_collection3d(Poly3DCollection([vertices], facecolor="#E8F0F5",
                                        alpha=.45, edgecolor="none"))
    for stage, color in enumerate(STAGES):
        controls = geometry(stage, np.linspace(0, 1, 160))[0]
        line3(ax, controls, color=color, lw=3.8)
        path_arrow(ax, stage, 0, color)
    for xyz, label, offset in (
        (vertices[0], r"$A_0$", [-.02, .10, .08]),
        (vertices[1], r"$A_1$", [.12, -.02, .08]),
        (vertices[2], r"$A_2$", [.10, .10, -.05]),
    ):
        ax.scatter(*xyz, color=INK, s=30, depthshade=False)
        ax.text(*(xyz + np.array(offset)), label, fontsize=15, ha="center")
    ax.text2D(.04, .11, r"$A_0\;(1,0,0)\ \to\ A_1\;(0,1,0)$",
              transform=ax.transAxes, fontsize=13)
    ax.text2D(.04, .055, r"$\to\ A_2\;(0,0,1)\ \to\ A_0$",
              transform=ax.transAxes, fontsize=13)
    ax.text2D(.04, -.015, r"$u_1=(\eta/C)u_2$；超导相位差保持 $\pi$",
              transform=ax.transAxes, fontsize=11, color=MUTED)


def tangent_frame(ax, n, f1, f2, scale=.26, plane=True):
    if plane:
        corners = [n + .22 * (i*f1 + j*f2)
                   for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
        ax.add_collection3d(Poly3DCollection([corners], facecolor="#E8EDF3",
                                            edgecolor="#C9D4DF", linewidth=.55,
                                            alpha=.5))
    arrow3(ax, n, f1, FRAME1, length=scale, width=2.6)
    arrow3(ax, n, f2, FRAME2, length=scale, width=2.6)
    ax.scatter(*n, s=13, color=INK, depthshade=False)


def draw_sphere(ax):
    panel_heading(ax, "B", "低能方向随路径平行输运", "球面法向是强耦合方向；两个箭头始终与它正交")
    ax.set_axis_off()
    ax.set_proj_type("ortho")
    ax.view_init(elev=-24, azim=40)
    ax.set_box_aspect((1, 1, 1))
    ax.set(xlim=(-.23, 1.3), ylim=(-.23, 1.3), zlim=(-1.3, .35))
    theta, phi = np.meshgrid(np.linspace(0, np.pi/2, 35),
                             np.linspace(0, np.pi/2, 35))
    ax.plot_surface(np.cos(theta), np.sin(theta)*np.cos(phi),
                    -np.sin(theta)*np.sin(phi), color="#BEDBE6",
                    alpha=.19, linewidth=0, shade=False)
    for angle in np.linspace(0, np.pi/2, 5):
        u = np.linspace(0, np.pi/2, 90)
        meridian = np.column_stack((np.cos(u), np.sin(u)*np.cos(angle),
                                    -np.sin(u)*np.sin(angle)))
        parallel = np.column_stack((np.full_like(u, np.cos(angle)),
                                    np.sin(angle)*np.cos(u),
                                    -np.sin(angle)*np.sin(u)))
        line3(ax, meridian, color=GRID, alpha=.55, lw=.65)
        line3(ax, parallel, color=GRID, alpha=.55, lw=.65)
    for vector, label in (([1, 0, 0], r"$+\gamma_b$"),
                          ([0, 1, 0], r"$+\gamma_2$"),
                          ([0, 0, -1], r"$-\gamma_3$")):
        arrow3(ax, np.zeros(3), vector, GRID, length=1.2, width=1.1)
        ax.text(*(1.28*np.asarray(vector)), label, fontsize=14, ha="center",
                color=MUTED)
    for stage, color in enumerate(STAGES):
        n = geometry(stage, np.linspace(0, 1, 180))[1]
        line3(ax, n, color=color, lw=3.5)
        path_arrow(ax, stage, 1, color)
        _, point, f1, f2 = geometry(stage, [.43])
        tangent_frame(ax, point[0], f1[0], f2[0], scale=.32)
    for stage, end in ((0, 0.), (0, 1.), (1, 1.)):
        point = geometry(stage, [end])[1][0]
        ax.scatter(*point, color=INK, s=23, depthshade=False)
    ax.text2D(.08, .09, r"$\mathbf{n}=(d,u_2,-u_3)/R$",
              transform=ax.transAxes, fontsize=15)
    ax.text2D(.08, .03, r"$\widetilde{\gamma}_2=\mathbf{f}_1\!\cdot\!(\gamma_b,\gamma_2,\gamma_3)$",
              transform=ax.transAxes, color=FRAME1, fontsize=13)
    ax.text2D(.08, -.02, r"$\widetilde{\gamma}_3=\mathbf{f}_2\!\cdot\!(\gamma_b,\gamma_2,\gamma_3)$",
              transform=ax.transAxes, color=FRAME2, fontsize=13)


def arrow2(ax, start, end, color, style="-", width=2.8):
    patch = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=17,
                            linewidth=width, linestyle=style, color=color)
    ax.add_patch(patch)


def draw_endpoint(ax):
    ax.set_axis_off()
    ax.set_aspect("equal")
    ax.set(xlim=(-1.55, 1.55), ylim=(-1.7, 1.7))
    ax.text(.01, 1.01, "C  回到同一点，标架旋转了", transform=ax.transAxes,
            fontsize=15, fontweight="bold")
    ax.text(.01, .955, "共同切平面：横向为原来的 γ₂，纵向为原来的 γ₃",
            transform=ax.transAxes, fontsize=10.5, color=MUTED)
    ax.plot([-1.13, 1.15], [0, 0], color=GRID, lw=1)
    ax.plot([0, 0], [-1.13, 1.15], color=GRID, lw=1)
    arrow2(ax, (0, 0), (.93, 0), FRAME1, style="--", width=2.3)
    arrow2(ax, (0, 0), (0, .93), FRAME2, style="--", width=2.3)
    arrow2(ax, (0, 0), (0, -.93), FRAME1)
    arrow2(ax, (0, -.035), (1.11, -.035), FRAME2)
    ax.scatter(0, 0, s=30, color=INK, zorder=5)
    ax.text(.57, .14, r"$\widetilde{\gamma}_2(0)=\gamma_2$", color=FRAME1,
            fontsize=12, ha="center")
    ax.text(-.1, 1.06, r"$\widetilde{\gamma}_3(0)=\gamma_3$", color=FRAME2,
            fontsize=12, ha="center")
    ax.text(.67, -.26, r"$\widetilde{\gamma}_3(T)=\gamma_2$", color=FRAME2,
            fontsize=12, ha="center")
    ax.text(-.1, -1.11, r"$\widetilde{\gamma}_2(T)=-\gamma_3$", color=FRAME1,
            fontsize=12, ha="center")
    t = np.linspace(np.pi/2, 0, 80)
    x, y = -.67 + .43*np.cos(t), -.25 + .43*np.sin(t)
    ax.plot(x, y, color=INK, lw=1.5)
    arrow2(ax, (x[-5], y[-5]), (x[-1], y[-1]), INK, width=1.5)
    ax.text(-.87, .36, r"$-\pi/2$", fontsize=17)
    ax.text(.5, .075, "虚线：起点    实线：终点", transform=ax.transAxes,
            ha="center", color=MUTED, fontsize=11)
    ax.text(.5, .01, r"$T=3\tau,\qquad \mathbf{n}(T)=\mathbf{n}(0)$",
            transform=ax.transAxes, ha="center", fontsize=13)


def main():
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(17.8, 8.5), facecolor="white")
    fig.text(.04, .952, "从耦合闭环到编织门：我们的平行输运标架",
             fontsize=23, fontweight="bold")
    fig.text(.04, .905, "原三段协议的理想强耦合参考几何；弱耦合的动态演化另由 D₁、D₂、D₃ 描述",
             fontsize=12, color=MUTED)
    axes = (
        fig.add_axes([.018, .255, .29, .59], projection="3d"),
        fig.add_axes([.315, .255, .35, .59], projection="3d"),
        fig.add_axes([.692, .255, .285, .59]),
    )
    draw_controls(axes[0])
    draw_sphere(axes[1])
    draw_endpoint(axes[2])
    handles = [Line2D([0], [0], color=c, lw=3, label=f"第 {n} 段")
               for n, c in zip(("I", "II", "III"), STAGES)]
    handles += [Line2D([0], [0], color=FRAME1, lw=3, label=r"低能方向 $\mathbf{f}_1$"),
                Line2D([0], [0], color=FRAME2, lw=3, label=r"低能方向 $\mathbf{f}_2$")]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.495, .213),
               ncol=5, frameon=False, fontsize=12, handlelength=2.2,
               columnspacing=2.)
    fig.add_artist(plt.Line2D([.04, .96], [.174, .174], transform=fig.transFigure,
                              color="#DCE4EA", lw=1))
    fig.text(.045, .126, r"局部平行：$\mathbf{f}_i\cdot\dot{\mathbf{f}}_j=0$",
             fontsize=15)
    fig.text(.39, .126, r"闭环结果：$\mathcal{Q}(T)=\mathcal{Q}(0)B$",
             fontsize=15)
    fig.text(.72, .126, r"$B=e^{i\pi X/4},\quad X=-i\gamma_2\gamma_3$",
             fontsize=16)
    fig.text(.045, .065, "箭头是 Majorana 线性组合的方向，不是量子态振幅；γ₁ 与旁观 γ₄ 未画出。",
             fontsize=11, color=MUTED)
    fig.text(.045, .031, "这是推导所用参考空间的几何图；并非完整五模在非零 E₁、t₁ 下无误差演化的示意。",
             fontsize=11, color=MUTED)
    for suffix in ("png", "pdf", "svg"):
        path = OUT / f"parallel_transport_frame.{suffix}"
        fig.savefig(path, dpi=190, facecolor="white")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
