const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType,
        Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType, PageBreak } = require("docx");

const REPO = "C:/Users/drric/TIC_TOC";
const MD = `${REPO}/manuscript/manuscript_tictoc.md`;

// ---- inline markdown -> TextRun[] ----
function inline(s) {
  s = s.replace(/\[([^\]]+)\]\([^)]+\)/g, "$1");   // [text](url) -> text
  const runs = [], re = /(\*\*([^*]+)\*\*|\*([^*]+)\*|`([^`]+)`)/g;
  let last = 0, m;
  while ((m = re.exec(s))) {
    if (m.index > last) runs.push(new TextRun(s.slice(last, m.index)));
    if (m[2] !== undefined) runs.push(new TextRun({ text: m[2], bold: true }));
    else if (m[3] !== undefined) runs.push(new TextRun({ text: m[3], italics: true }));
    else if (m[4] !== undefined) runs.push(new TextRun({ text: m[4], font: "Consolas", size: 20 }));
    last = re.lastIndex;
  }
  if (last < s.length) runs.push(new TextRun(s.slice(last)));
  return runs.length ? runs : [new TextRun(s)];
}
const para = (t) => new Paragraph({ children: inline(t), spacing: { after: 120, line: 276 } });
const bullet = (t) => new Paragraph({ children: inline(t), bullet: { level: 0 }, spacing: { after: 60 } });
const quote = (t) => new Paragraph({ children: inline(t).map(r => r), indent: { left: 360 },
  spacing: { after: 120 }, border: { left: { style: BorderStyle.SINGLE, size: 12, color: "999999", space: 8 } } });

function heading(t) {
  const h = t.match(/^#+/)[0].length, txt = t.replace(/^#+\s*/, "");
  if (h === 1) return new Paragraph({ children: [new TextRun({ text: txt, bold: true, size: 34 })],
    spacing: { after: 200 } });
  return new Paragraph({ heading: h === 2 ? HeadingLevel.HEADING_1 : HeadingLevel.HEADING_2,
    children: inline(txt) });
}

// ---- parse markdown body ----
const lines = fs.readFileSync(MD, "utf8").split(/\r?\n/);
const body = []; let buf = [], i = 0;
const flush = () => { if (buf.length) { body.push(para(buf.join(" "))); buf = []; } };
while (i < lines.length) {
  const raw = lines[i], t = raw.trim();
  if (t.startsWith("<!--")) { while (i < lines.length && !lines[i].includes("-->")) i++; i++; continue; }
  if (t === "" || t === "---") { flush(); i++; continue; }
  if (t.startsWith("#")) { flush(); body.push(heading(t)); i++; continue; }
  if (/^[-*]\s+/.test(t)) { flush(); body.push(bullet(t.replace(/^[-*]\s+/, ""))); i++; continue; }
  if (t.startsWith(">")) { flush(); body.push(quote(t.replace(/^>\s?/, ""))); i++; continue; }
  buf.push(t); i++;
}
flush();

// ---- figures ----
const pngSize = (f) => { const b = fs.readFileSync(f); return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) }; };
const FIGS = [
  ["Figure 1", `${REPO}/manuscript/figures/Fig1_design.png`, "Experimental design and analysis workflow."],
  ["Figure 2", `${REPO}/manuscript/figures/Fig2_root_architecture.png`, "Root architecture over days 3–6 (Flight vs Ground by genotype)."],
  ["Figure 3", `${REPO}/manuscript/figures/Fig3b_DEG_counts.png`, "Spaceflight DEG counts per genotype × tissue (up/down)."],
  ["Figure 4", `${REPO}/manuscript/Fig4_module_trait_named.png`, "Root co-expression module–trait correlations (group-level)."],
  ["Figure 5", `${REPO}/manuscript/figures/Fig5_GO_enrichment.png`, "GO biological-process enrichment across key contrasts."],
  ["Figure 6", `${REPO}/manuscript/figures/Fig6_physioscores.png`, "PhysioSpace stress-pattern decoding (stress axis × genotype × tissue)."],
  ["Figure 7", `${REPO}/manuscript/figures/Fig7_integrative_model.png`, "Integrative model — AVP-OX as a stress attenuator."],
  ["Figure 8", `${REPO}/manuscript/figures/Fig8_program_atlas.png`, "Autoencoder expression-program atlas (bulk RNA-seq; programs, not cells)."],
  ["Figure 9", `${REPO}/manuscript/figures/Fig9_program_stress_WTvsAVPOX.png`, "Program stress response, WT vs AVP-OX."],
];
const figEls = [new Paragraph({ children: [new PageBreak()] }),
  new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun("Figures")] })];
for (const [name, file, cap] of FIGS) {
  const { w, h } = pngSize(file), W = Math.min(600, w), H = Math.round(h * W / w);
  figEls.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160 },
    children: [new ImageRun({ type: "png", data: fs.readFileSync(file), transformation: { width: W, height: H },
      altText: { title: name, description: cap, name } })] }));
  figEls.push(new Paragraph({ spacing: { after: 160 }, children: [
    new TextRun({ text: name + ". ", bold: true }), new TextRun(cap)] }));
}

// ---- Table 2: DEG counts ----
const csv = fs.readFileSync(`${REPO}/deseq2/contrasts/DEG_counts_summary.csv`, "utf8").trim().split(/\r?\n/)
  .map(r => r.split(",").map(c => c.replace(/^"|"$/g, "")));
const bd = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: bd, bottom: bd, left: bd, right: bd };
const cell = (txt, hdr, wdt) => new TableCell({ borders, width: { size: wdt, type: WidthType.DXA },
  shading: hdr ? { fill: "D5E8F0", type: ShadingType.CLEAR } : undefined,
  margins: { top: 60, bottom: 60, left: 100, right: 100 },
  children: [new Paragraph({ children: [new TextRun({ text: txt, bold: !!hdr })] })] });
const cw = [4200, 2000, 1580, 1580];
const t2 = new Table({ width: { size: 9360, type: WidthType.DXA }, columnWidths: cw,
  rows: csv.map((r, ri) => new TableRow({ children: [
    cell(r[0], ri === 0, cw[0]), cell(r[1], ri === 0, cw[1]), cell(r[2], ri === 0, cw[2]), cell(r[3], ri === 0, cw[3]) ] })) });
const tableEls = [new Paragraph({ children: [new PageBreak()] }),
  new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun("Table 2. DEG counts per contrast")] }),
  new Paragraph({ spacing: { after: 120 }, children: [new TextRun({ italics: true,
    text: "DESeq2, P.adj < 0.05, |log2FC| >= 1 (ashr-shrunk). Root responds far more than shoot; AVP-OX interaction contrasts down-biased." })] }),
  t2];

// ---- assemble ----
const doc = new Document({
  styles: { default: { document: { run: { font: "Arial", size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: "Arial", color: "12324f" },
        paragraph: { spacing: { before: 260, after: 140 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: "Arial", color: "1b4a63" },
        paragraph: { spacing: { before: 180, after: 100 }, outlineLevel: 1 } } ] },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    children: [...body, ...tableEls, ...figEls]
  }]
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(`${REPO}/manuscript/TICTOC_manuscript_draft.docx`, b);
  console.log("wrote TICTOC_manuscript_draft.docx", b.length, "bytes;", body.length, "body elements,", FIGS.length, "figures"); });
