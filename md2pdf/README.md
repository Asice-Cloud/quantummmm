# md2pdf — 用 Slidev 把 Markdown 转成逐页 PPT / PDF

基于 [Slidev](https://sli.dev)（slide + dev）。把 Markdown 按 `---` 分页写成幻灯片，
开发时浏览器实时预览，最后导出成**一页一张幻灯片**的 PDF / PPTX / PNG。

## 目录结构

```
md2pdf/
├── slides.md        # 幻灯片源文件（唯一需要改的文件）
├── package.json     # 脚本与依赖
├── components/      # 可选：自定义 Vue 组件
├── public/          # 可选：图片等静态资源
└── style.css        # 可选：自定义样式
```

## 安装

```bash
npm install
# 首次导出需要下载 Chromium（约 150MB）
npx playwright install chromium
```

## 开发预览

```bash
npm run dev          # 打开 http://localhost:3030，热更新
```

## 导出

```bash
npm run export:pdf              # 生成 ./slides-export.pdf（逐页）
npm run export:pptx             # 生成真正的 PPTX（可再编辑）
npm run export:png              # 每页一张 PNG

# 直接传参数
npx slidev export --format pdf -o out.pdf
npx slidev export --with-clicks # 动画的每一步也各导一页
npx slidev export --dark        # 深色主题
```

也可以不用命令行：`npm run dev` 后打开 `http://localhost:3030/export`，
在网页 UI 里点击导出。

## Markdown 写法

```markdown
---
theme: default
title: 我的汇报
---

# 第一页标题

正文内容，支持 **粗体**、`代码`、$公式$、表格、图片

---

# 第二页

- 要点一
- 要点二
```

要点：

- 全局配置写在文件**最开头**的 frontmatter（`---` 包裹）
- 每张幻灯片之间用单独一行的 `---` 分隔
- 单页内可用 `<style>`、`<div>`、Vue 组件做精细排版
- 中文若显示异常，在 frontmatter 里配置 `fonts` 或在 `style.css` 指定中文字体

## 把已有 md 改造成 slides.md

1. 把原文档标题/章节用 `---` 分隔成独立页面
2. 每页只保留一个主题，内容过多就拆成两页
3. 需要图表就插入 `![alt](./public/xxx.png)`，把图片放进 `public/`
4. 代码块、公式、Mermaid 直接沿用原 Markdown

## 常见问题

| 问题 | 解决 |
|---|---|
| 导出报 Playwright 错误 | 执行 `npx playwright install chromium` |
| 中文乱码/方块 | frontmatter 配置 `fonts:` 或 `style.css` 指定中文字体 |
| 一页内容溢出 | 拆成两页，或调小字号 |
| 想要可编辑 PPT | `npm run export:pptx`（复杂布局可能偏差） |
