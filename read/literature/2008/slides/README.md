# Majorana CNOT 文献汇报幻灯片

根据 `read/literature/2008.md` 整理。16:9，中文 Beamer，13 页主报告（含封面）及 2 页备用内容，适合约 10–15 分钟汇报。主报告含 Fig. 3(b) 的专页解释，以及附录 A、B 的简要介绍。

## Overleaf 使用

1. 新建项目，选择 **Upload Project**，上传完整的 `2008_cnot_overleaf.zip`。不能只复制 `main.tex`：五张配图也需要上传。
2. 在项目设置中，将 **Compiler** 设为 **XeLaTeX**，主文件设为 **main.tex**。源文件开头的 `% !TeX program = xelatex` 注释不会自动改变 Overleaf 的编译器设置。
3. 确认项目根目录中有 `main.tex`、`latexmkrc` 和五张 PNG，所有文件位于同一级目录。
4. 点击 **Recompile**。如果之前使用了 pdfLaTeX，可选择 **Recompile from scratch** 后重新编译。

中文使用 TeX Live 自带的 Fandol 字体，不需要上传字体或开启 shell escape。图片均为项目内的相对路径；TikZ 电路和流程图直接由 TeX 绘制，无需外部转换。

源文件会检查编译引擎，避免选错引擎后继续产生大量字体错误。上传包使用平铺结构，图片与 `main.tex` 同级；本地仓库仍把图片放在 `figures/` 中，源文件兼容两种布局。

如果要修复已有 Overleaf 项目，先在电脑上解压 ZIP，再把其中的 `main.tex`、`latexmkrc` 和五张 PNG 一起上传到项目根目录，覆盖旧的同名文件。只替换 TeX 文件无法解决缺图问题。

上传后的文件结构应为：

```text
main.tex
latexmkrc
README.md
paper_fig1.png
paper_fig3b.png
paper_fig4.png
paper_fig6b.png
double_tunneling_channels.png
```

## 编译日志中的常见问题

- `This is pdfTeX ... preloaded format=pdflatex`、`CTeX fontset 'fandol' is unavailable in current mode`、`Font uniyou... not found`：当前选用了 pdfLaTeX，应在项目设置中改为 XeLaTeX，不必删除中文或更换字体。
- `File 'paper_fig1.png' not found`、`Unable to load picture or PDF file` 等：编译器没有找到图片。将五张 PNG 上传到 `main.tex` 同级目录。随后的 `Package graphics Error: Division by 0` 是图片尺寸无法读取引发的连带错误。
- `Overfull \vbox`：内容超出可用页高的排版警告。已根据 XeLaTeX 日志减少标题、图示和段落留白，并为较密的页面加入自动缩放；修改内容后仍应检查最终 PDF。

## 本地编译

也可以在具有完整 TeX Live 的本地环境中，在本目录连续运行两次：

```bash
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

主报告最后一页后为备用内容：“有效二能级基态”与参考文献。源文件按页编号注释，便于修改；封面未填写汇报人，可自行添加。

已根据用户提供的 XeLaTeX 日志修订，并检查 TeX 结构、图片路径和压缩包完整性。当前工作环境缺少 Beamer、中文宏包及 XeLaTeX 格式文件，修订版未完成本地 PDF 编译；最终版式请以 Overleaf 编译结果为准。

## 本地仓库中的文件

- `main.tex`：完整幻灯片源文件。
- `latexmkrc`：供 latexmk 使用的 XeLaTeX 编译配置；Overleaf 项目设置仍应选择 XeLaTeX。
- `figures/paper_fig1.png`：原论文 Fig. 1。
- `figures/paper_fig3b.png`：原论文 Fig. 3(b)，双隧穿协议的瞬时能谱。
- `figures/paper_fig4.png`：原论文 Fig. 4。
- `figures/paper_fig6b.png`：原论文 Fig. 6(b)，不对称通道下的保真度。
- `figures/double_tunneling_channels.png`：已有双隧穿示意图的通道结构部分，根据论文模型绘制。

打包时，以上五张图片均放入 ZIP 根目录，与 `main.tex` 同级。

论文：A. Calzona, N. P. Bauer, B. Trauzettel, *Holonomic implementation of CNOT gate on topological Majorana qubits*, SciPost Phys. Core **3**, 014 (2020)，[DOI](https://doi.org/10.21468/SciPostPhysCore.3.2.014)。论文属于理论方案与数值模拟研究；幻灯片中的保真度图为原文结果。
