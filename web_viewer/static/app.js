const canvas = document.getElementById("plot");
const ctx = canvas.getContext("2d");
const dotplotCanvas = document.getElementById("dotplot");
const dotCtx = dotplotCanvas.getContext("2d");

const MAX_SCALE = 0.3;
const MIN_SCALE = 0.0005;
const PANEL_PADDING = 160;
const PANEL_GUTTER = 420;
const PANEL_ROW_GUTTER = 520;
const PANEL_LABEL_HEIGHT = 220;
const OUTLINE_MIN_SCREEN_PX = 0.15;
const OUTLINE_MAX_SCREEN_PX = 1.25;
const OUTLINE_BASE_SCREEN_PX = 0.22;

const state = {
  samples: {},
  expressions: {},
  currentSample: "S1",
  currentGene: "",
  expressionRange: { vmin: 0, vmax: 0 },
  view: { scale: 1, x: 0, y: 0 },
  dragging: false,
  lastPointer: null,
  devicePixelRatio: window.devicePixelRatio || 1,
  drawPending: false,
  dotplot: { payload: null, error: "", loading: false },
  dotplotDrawPending: false,
};

const els = {
  form: document.getElementById("geneForm"),
  input: document.getElementById("geneInput"),
  geneList: document.getElementById("geneList"),
  status: document.getElementById("status"),
  currentGene: document.getElementById("currentGene"),
  currentSample: document.getElementById("currentSample"),
  cellCount: document.getElementById("cellCount"),
  nonzeroCount: document.getElementById("nonzeroCount"),
  rangeText: document.getElementById("rangeText"),
  legendMax: document.getElementById("legendMax"),
  legendMin: document.getElementById("legendMin"),
  scaleText: document.getElementById("scaleText"),
};

const REDS = [
  [255, 245, 240],
  [254, 224, 210],
  [252, 187, 161],
  [252, 146, 114],
  [251, 106, 74],
  [239, 59, 44],
  [203, 24, 29],
  [165, 15, 21],
  [103, 0, 13],
];

function setStatus(message) {
  els.status.textContent = message;
}

function formatNumber(value, digits = 3) {
  if (!Number.isFinite(value)) return "-";
  if (Math.abs(value) >= 1000) return value.toLocaleString();
  return Number(value.toFixed(digits)).toString();
}

function tileKey(x, y) {
  return `${x},${y}`;
}

function bboxWidth(bbox) {
  return bbox[2] - bbox[0] + 1;
}

function bboxHeight(bbox) {
  return bbox[3] - bbox[1] + 1;
}

function bboxIntersects(a, b) {
  return !(a[2] < b.left || a[0] > b.right || a[3] < b.top || a[1] > b.bottom);
}

function rectsIntersect(a, b) {
  return !(a.right < b.left || a.left > b.right || a.bottom < b.top || a.top > b.bottom);
}

function panelColumns(sample, replicateCount) {
  if (sample === "S1") return 4;
  return Math.max(1, replicateCount);
}

function prepareSpatial(manifest, replicatePayload) {
  const tileMap = new Map();
  for (const tile of manifest.tiles || []) {
    tileMap.set(tileKey(tile.x, tile.y), tile);
  }

  const replicates = replicatePayload.replicates || [];
  const panels = [];
  const cellToPanel = new Map();
  const tileToPanels = new Map();
  const referenceRep = replicates.find((rep) => rep.id === `${manifest.sample.toLowerCase()}_rep1`) || replicates[0];
  const displayWidth = referenceRep ? bboxWidth(referenceRep.bbox) : manifest.width;
  const displayHeight = referenceRep ? bboxHeight(referenceRep.bbox) : manifest.height;
  const columns = panelColumns(manifest.sample, replicates.length);
  const rows = Math.max(1, Math.ceil(replicates.length / columns));
  let assignedCellCount = 0;

  for (const [index, rep] of replicates.entries()) {
    const col = index % columns;
    const row = Math.floor(index / columns);
    const sourceWidth = bboxWidth(rep.bbox);
    const sourceHeight = bboxHeight(rep.bbox);
    const panel = {
      ...rep,
      x: PANEL_PADDING + col * (displayWidth + PANEL_GUTTER),
      y: PANEL_LABEL_HEIGHT + PANEL_PADDING + row * (displayHeight + PANEL_ROW_GUTTER),
      width: displayWidth,
      height: displayHeight,
      sourceWidth,
      sourceHeight,
      scaleX: displayWidth / sourceWidth,
      scaleY: displayHeight / sourceHeight,
    };
    panel.bounds = {
      left: panel.x,
      top: panel.y,
      right: panel.x + displayWidth,
      bottom: panel.y + displayHeight,
    };

    panels.push(panel);
    assignedCellCount += rep.assignedCellCount || rep.cellIds.length;

    for (const cellId of rep.cellIds) {
      cellToPanel.set(cellId, panel);
    }
    for (const key of rep.tileKeys || []) {
      if (!tileToPanels.has(key)) tileToPanels.set(key, []);
      tileToPanels.get(key).push(panel);
    }

  }

  const layoutWidth = Math.max(
    1,
    panels.length ? PANEL_PADDING * 2 + columns * displayWidth + (columns - 1) * PANEL_GUTTER : manifest.width
  );
  const layoutHeight = panels.length
    ? PANEL_LABEL_HEIGHT + PANEL_PADDING * 2 + rows * displayHeight + (rows - 1) * PANEL_ROW_GUTTER
    : manifest.height;

  return {
    ...manifest,
    panels,
    panelColumns: columns,
    panelRows: rows,
    cellToPanel,
    tileToPanels,
    layoutWidth,
    layoutHeight,
    assignedCellCount,
    tileMap,
    loadedTiles: new Map(),
    loadingTiles: new Set(),
    failedTiles: new Set(),
    pathCache: new Map(),
  };
}

function requestDraw() {
  if (state.drawPending) return;
  state.drawPending = true;
  window.requestAnimationFrame(() => {
    state.drawPending = false;
    draw();
  });
}

function requestDotplotDraw() {
  if (state.dotplotDrawPending) return;
  state.dotplotDrawPending = true;
  window.requestAnimationFrame(() => {
    state.dotplotDrawPending = false;
    drawDotplot();
  });
}

function resizeCanvasElement(targetCanvas) {
  const rect = targetCanvas.getBoundingClientRect();
  targetCanvas.width = Math.max(1, Math.floor(rect.width * state.devicePixelRatio));
  targetCanvas.height = Math.max(1, Math.floor(rect.height * state.devicePixelRatio));
}

function resizeCanvas() {
  state.devicePixelRatio = window.devicePixelRatio || 1;
  resizeCanvasElement(canvas);
  resizeCanvasElement(dotplotCanvas);
  requestDraw();
  requestDotplotDraw();
}

function currentSpatial() {
  return state.samples[state.currentSample];
}

function currentExpression() {
  return state.expressions[state.currentSample] || { map: new Map(), max: 0, min: 0, nonzero: 0 };
}

function fitView() {
  const spatial = currentSpatial();
  if (!spatial) return;

  const rect = canvas.getBoundingClientRect();
  const margin = 36;
  const sx = Math.max(1, rect.width - margin * 2) / spatial.layoutWidth;
  const sy = Math.max(1, rect.height - margin * 2) / spatial.layoutHeight;
  const scale = Math.max(MIN_SCALE, Math.min(sx, sy));

  spatial.fitScale = scale;
  state.view.scale = scale;
  state.view.x = (rect.width - spatial.layoutWidth * scale) / 2;
  state.view.y = (rect.height - spatial.layoutHeight * scale) / 2;
  updateScaleText();
}

function updateScaleText() {
  const spatial = currentSpatial();
  if (!spatial) {
    els.scaleText.textContent = "-";
    return;
  }
  const fitted = Math.min(canvas.clientWidth / spatial.layoutWidth, canvas.clientHeight / spatial.layoutHeight);
  els.scaleText.textContent = `${Math.round((state.view.scale / fitted) * 100)}%`;
}

function colorForValue(value, range) {
  const min = Number.isFinite(range.vmin) ? range.vmin : 0;
  const max = Number.isFinite(range.vmax) ? range.vmax : 0;
  const t = max > min ? Math.max(0, Math.min(1, (value - min) / (max - min))) : 0;
  const scaled = t * (REDS.length - 1);
  const idx = Math.min(REDS.length - 2, Math.floor(scaled));
  const local = scaled - idx;
  const a = REDS[idx];
  const b = REDS[idx + 1];
  return a.map((component, i) => Math.round(component + (b[i] - component) * local));
}

function rgbCss(rgb) {
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

function invertRgb(rgb) {
  return [255 - rgb[0], 255 - rgb[1], 255 - rgb[2]];
}

function outlineScreenWidth(spatial, scale) {
  const fitScale = spatial.fitScale || scale || 1;
  const relativeZoom = scale / fitScale;
  return Math.max(
    OUTLINE_MIN_SCREEN_PX,
    Math.min(OUTLINE_MAX_SCREEN_PX, OUTLINE_BASE_SCREEN_PX * relativeZoom)
  );
}

function createCellPath(cell) {
  const path = new Path2D();
  const x0 = cell.bbox[0];
  const y0 = cell.bbox[1];

  for (const contour of cell.contours || []) {
    if (!contour.length) continue;
    path.moveTo(x0 + contour[0][0], y0 + contour[0][1]);
    for (let i = 1; i < contour.length; i += 1) {
      path.lineTo(x0 + contour[i][0], y0 + contour[i][1]);
    }
    path.closePath();
  }

  return path;
}

function getCellPath(spatial, cell) {
  let path = spatial.pathCache.get(cell.id);
  if (!path) {
    path = createCellPath(cell);
    spatial.pathCache.set(cell.id, path);
  }
  return path;
}

function visibleBounds(rect) {
  const { scale, x, y } = state.view;
  return {
    left: -x / scale,
    top: -y / scale,
    right: (-x + rect.width) / scale,
    bottom: (-y + rect.height) / scale,
  };
}

function originalBoundsForPanel(panel, bounds) {
  const left = Math.max(panel.bounds.left, bounds.left);
  const top = Math.max(panel.bounds.top, bounds.top);
  const right = Math.min(panel.bounds.right, bounds.right);
  const bottom = Math.min(panel.bounds.bottom, bounds.bottom);
  return {
    left: Math.max(panel.bbox[0], panel.bbox[0] + (left - panel.x) / panel.scaleX),
    top: Math.max(panel.bbox[1], panel.bbox[1] + (top - panel.y) / panel.scaleY),
    right: Math.min(panel.bbox[2], panel.bbox[0] + (right - panel.x) / panel.scaleX),
    bottom: Math.min(panel.bbox[3], panel.bbox[1] + (bottom - panel.y) / panel.scaleY),
  };
}

function tileKeysForOriginalBounds(spatial, bounds) {
  const keys = [];
  const maxTileX = Math.ceil(spatial.width / spatial.tileSize) - 1;
  const maxTileY = Math.ceil(spatial.height / spatial.tileSize) - 1;
  const tx0 = Math.max(0, Math.floor(Math.max(0, bounds.left) / spatial.tileSize));
  const ty0 = Math.max(0, Math.floor(Math.max(0, bounds.top) / spatial.tileSize));
  const tx1 = Math.min(maxTileX, Math.floor(Math.min(spatial.width - 1, bounds.right) / spatial.tileSize));
  const ty1 = Math.min(maxTileY, Math.floor(Math.min(spatial.height - 1, bounds.bottom) / spatial.tileSize));

  for (let ty = ty0; ty <= ty1; ty += 1) {
    for (let tx = tx0; tx <= tx1; tx += 1) {
      const key = tileKey(tx, ty);
      if (spatial.tileMap.has(key)) keys.push(key);
    }
  }
  return keys;
}

function visibleTileKeys(spatial, bounds) {
  const keys = new Set();
  for (const panel of spatial.panels) {
    if (!rectsIntersect(panel.bounds, bounds)) continue;
    const originalBounds = originalBoundsForPanel(panel, bounds);
    for (const key of tileKeysForOriginalBounds(spatial, originalBounds)) {
      keys.add(key);
    }
  }
  return Array.from(keys);
}

async function loadTile(spatial, key) {
  const tile = spatial.tileMap.get(key);
  if (!tile || spatial.loadedTiles.has(key) || spatial.loadingTiles.has(key) || spatial.failedTiles.has(key)) {
    return;
  }

  spatial.loadingTiles.add(key);
  try {
    const response = await fetch(tile.url);
    if (!response.ok) throw new Error(`无法加载轮廓块 ${tile.url}`);
    const payload = await response.json();
    spatial.loadedTiles.set(key, payload);
  } catch (error) {
    spatial.failedTiles.add(key);
    console.error(error);
    setStatus(error.message);
  } finally {
    spatial.loadingTiles.delete(key);
    requestDraw();
  }
}

function ensureTiles(spatial, keys) {
  let missing = 0;
  for (const key of keys) {
    if (!spatial.loadedTiles.has(key) && !spatial.failedTiles.has(key)) {
      missing += 1;
      loadTile(spatial, key);
    }
  }
  return missing;
}

function translatedBBox(cell, panel) {
  return [
    panel.x + (cell.bbox[0] - panel.bbox[0]) * panel.scaleX,
    panel.y + (cell.bbox[1] - panel.bbox[1]) * panel.scaleY,
    panel.x + (cell.bbox[2] - panel.bbox[0] + 1) * panel.scaleX,
    panel.y + (cell.bbox[3] - panel.bbox[1] + 1) * panel.scaleY,
  ];
}

function drawPanels(spatial, bounds, scale) {
  ctx.save();
  ctx.lineWidth = 1 / scale;

  for (const panel of spatial.panels) {
    if (!rectsIntersect(panel.bounds, bounds)) continue;

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(panel.x, panel.y, panel.width, panel.height);
    ctx.strokeStyle = "#b8c1cc";
    ctx.strokeRect(panel.x, panel.y, panel.width, panel.height);
  }

  ctx.restore();
}

function drawPanelLabels(spatial, bounds, scale, viewX, viewY, rect) {
  ctx.save();
  ctx.font = "13px Arial, \"Noto Sans SC\", sans-serif";
  ctx.textBaseline = "bottom";

  for (const panel of spatial.panels) {
    if (!rectsIntersect(panel.bounds, bounds)) continue;

    const screenLeft = viewX + panel.x * scale;
    const screenTop = viewY + panel.y * scale;
    const screenRight = viewX + (panel.x + panel.width) * scale;
    if (screenRight < 0 || screenLeft > rect.width || screenTop < 0 || screenTop > rect.height) {
      continue;
    }

    ctx.fillStyle = "rgba(255,255,255,0.88)";
    ctx.fillRect(screenLeft, Math.max(0, screenTop - 23), 170, 20);
    ctx.fillStyle = "#253044";
    ctx.fillText(
      `${panel.label} · ${panel.assignedCellCount.toLocaleString()} cells`,
      screenLeft + 6,
      Math.max(15, screenTop - 7)
    );
  }

  ctx.restore();
}

function drawCells(spatial, expression, keys, bounds, scale) {
  const drawn = new Set();
  let drawnCells = 0;
  const outlinePx = outlineScreenWidth(spatial, scale);

  ctx.lineJoin = "round";
  ctx.lineCap = "round";

  for (const key of keys) {
    const tile = spatial.loadedTiles.get(key);
    if (!tile) continue;

    for (const cell of tile.cells || []) {
      if (drawn.has(cell.id)) continue;
      const panel = spatial.cellToPanel.get(cell.id);
      if (!panel) continue;

      const layoutBBox = translatedBBox(cell, panel);
      if (!bboxIntersects(layoutBBox, bounds)) continue;
      drawn.add(cell.id);

      const value = expression.map.get(cell.id) || 0;
      const rgb = colorForValue(value, state.expressionRange);
      const path = getCellPath(spatial, cell);

      ctx.save();
      ctx.translate(panel.x, panel.y);
      ctx.scale(panel.scaleX, panel.scaleY);
      ctx.translate(-panel.bbox[0], -panel.bbox[1]);
      ctx.lineWidth = outlinePx / (scale * Math.max(panel.scaleX, panel.scaleY));
      ctx.fillStyle = rgbCss(rgb);
      ctx.strokeStyle = rgbCss(invertRgb(rgb));
      ctx.fill(path);
      ctx.stroke(path);
      ctx.restore();

      drawnCells += 1;
    }
  }

  return drawnCells;
}

function draw() {
  const spatial = currentSpatial();
  const rect = canvas.getBoundingClientRect();
  ctx.setTransform(state.devicePixelRatio, 0, 0, state.devicePixelRatio, 0, 0);
  ctx.clearRect(0, 0, rect.width, rect.height);
  ctx.fillStyle = "#eef1f4";
  ctx.fillRect(0, 0, rect.width, rect.height);

  if (!spatial) return;

  const expression = currentExpression();
  const bounds = visibleBounds(rect);
  const keys = visibleTileKeys(spatial, bounds);
  const missingTiles = ensureTiles(spatial, keys);
  const { scale, x, y } = state.view;

  ctx.save();
  ctx.translate(x, y);
  ctx.scale(scale, scale);

  drawPanels(spatial, bounds, scale);
  const drawnCells = drawCells(spatial, expression, keys, bounds, scale);

  ctx.restore();
  drawPanelLabels(spatial, bounds, scale, x, y, rect);

  updateScaleText();
  if (missingTiles > 0) {
    setStatus(`加载 ${spatial.sample} 轮廓块 ${spatial.loadedTiles.size}/${spatial.tileCount}...`);
  } else if (state.currentGene) {
    setStatus(`已加载 ${state.currentGene}，当前视图绘制 ${drawnCells.toLocaleString()} 个细胞`);
  }
}

function ellipsizeText(context, text, maxWidth) {
  if (context.measureText(text).width <= maxWidth) return text;
  const ellipsis = "...";
  let low = 0;
  let high = text.length;
  while (low < high) {
    const mid = Math.ceil((low + high) / 2);
    const candidate = `${text.slice(0, mid)}${ellipsis}`;
    if (context.measureText(candidate).width <= maxWidth) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }
  return `${text.slice(0, low)}${ellipsis}`;
}

function dotRadius(pctExpr, pctRange) {
  const minRadius = 3;
  const maxRadius = 15;
  if (!Number.isFinite(pctExpr)) return minRadius;
  if (!pctRange || pctRange.max <= pctRange.min) return (minRadius + maxRadius) / 2;
  const t = Math.max(0, Math.min(1, (pctExpr - pctRange.min) / (pctRange.max - pctRange.min)));
  return minRadius + t * (maxRadius - minRadius);
}

function drawDotplotMessage(rect, message) {
  dotCtx.fillStyle = "#ffffff";
  dotCtx.fillRect(0, 0, rect.width, rect.height);
  dotCtx.fillStyle = "#667085";
  dotCtx.font = "13px Arial, \"Noto Sans SC\", sans-serif";
  dotCtx.textAlign = "center";
  dotCtx.textBaseline = "middle";
  dotCtx.fillText(message, rect.width / 2, rect.height / 2);
}

function scaledColorForValue(value, range) {
  const min = Number.isFinite(range.vmin) ? range.vmin : -2.5;
  const max = Number.isFinite(range.vmax) ? range.vmax : 2.5;
  const t = max > min ? Math.max(0, Math.min(1, (value - min) / (max - min))) : 0.5;
  const scaled = t * (REDS.length - 1);
  const idx = Math.min(REDS.length - 2, Math.floor(scaled));
  const local = scaled - idx;
  const a = REDS[idx];
  const b = REDS[idx + 1];
  return a.map((component, i) => Math.round(component + (b[i] - component) * local));
}

function drawDotplotLegend(colorRange, pctRange, right, centerY) {
  const legendWidth = 116;
  const legendHeight = 10;
  const x = Math.max(170, right - legendWidth);
  const y = Math.max(28, centerY - 36);
  const gradient = dotCtx.createLinearGradient(x, 0, x + legendWidth, 0);

  for (let i = 0; i < REDS.length; i += 1) {
    gradient.addColorStop(i / (REDS.length - 1), rgbCss(REDS[i]));
  }

  dotCtx.fillStyle = "#253044";
  dotCtx.font = "11px Arial, \"Noto Sans SC\", sans-serif";
  dotCtx.textAlign = "left";
  dotCtx.textBaseline = "bottom";
  dotCtx.fillText("Scaled avg expr", x, y - 4);
  dotCtx.fillStyle = gradient;
  dotCtx.fillRect(x, y, legendWidth, legendHeight);
  dotCtx.strokeStyle = "#d0b4aa";
  dotCtx.strokeRect(x, y, legendWidth, legendHeight);
  dotCtx.fillStyle = "#667085";
  dotCtx.textBaseline = "top";
  dotCtx.fillText(formatNumber(colorRange.vmin, 1), x, y + legendHeight + 3);
  dotCtx.textAlign = "right";
  dotCtx.fillText(formatNumber(colorRange.vmax, 1), x + legendWidth, y + legendHeight + 3);

  const sizeY = y + 64;
  const sizes = pctRange.max > pctRange.min
    ? [
        pctRange.min,
        pctRange.min + (pctRange.max - pctRange.min) / 2,
        pctRange.max,
      ]
    : [pctRange.min];
  dotCtx.textAlign = "left";
  dotCtx.textBaseline = "middle";
  dotCtx.fillStyle = "#253044";
  dotCtx.fillText("% cells", x, sizeY - 22);
  for (const [index, pct] of sizes.entries()) {
    const cx = sizes.length === 1 ? x + 56 : x + 12 + index * 40;
    const radius = dotRadius(pct, pctRange);
    dotCtx.beginPath();
    dotCtx.arc(cx, sizeY, radius, 0, Math.PI * 2);
    dotCtx.fillStyle = "#fcbba1";
    dotCtx.fill();
    dotCtx.strokeStyle = "#8a1d18";
    dotCtx.stroke();
    dotCtx.fillStyle = "#667085";
    dotCtx.textAlign = "center";
    dotCtx.fillText(formatNumber(pct, 1), cx, sizeY + 21);
  }
}

function drawDotplot() {
  const rect = dotplotCanvas.getBoundingClientRect();
  dotCtx.setTransform(state.devicePixelRatio, 0, 0, state.devicePixelRatio, 0, 0);
  dotCtx.clearRect(0, 0, rect.width, rect.height);

  if (state.dotplot.loading) {
    drawDotplotMessage(rect, "加载 cluster dotplot...");
    return;
  }
  if (state.dotplot.error) {
    drawDotplotMessage(rect, state.dotplot.error);
    return;
  }

  const payload = state.dotplot.payload;
  const clusters = payload ? payload.clusters || [] : [];
  if (!clusters.length) {
    drawDotplotMessage(rect, "暂无 dotplot 数据");
    return;
  }

  dotCtx.fillStyle = "#ffffff";
  dotCtx.fillRect(0, 0, rect.width, rect.height);

  const compact = rect.width < 680;
  const left = compact ? 118 : 158;
  const right = compact ? rect.width - 24 : rect.width - 176;
  const titleY = 20;
  const plotWidth = Math.max(1, right - left);
  const centerY = Math.min(rect.height - 42, titleY + 44);
  const band = plotWidth / clusters.length;
  const scaledValues = clusters.map((cluster) => Number(cluster.avgExprScaled)).filter(Number.isFinite);
  const pctValues = clusters.map((cluster) => Number(cluster.pctExpr)).filter(Number.isFinite);
  const colorRange = scaledValues.length
    ? { vmin: Math.min(...scaledValues), vmax: Math.max(...scaledValues) }
    : { vmin: 0, vmax: Math.max(0, ...clusters.map((cluster) => Number(cluster.avgExpr) || 0)) };
  const pctRange = pctValues.length
    ? { min: Math.min(...pctValues), max: Math.max(...pctValues) }
    : { min: 0, max: 100 };

  dotCtx.fillStyle = "#253044";
  dotCtx.font = "700 18px Arial, \"Noto Sans SC\", sans-serif";
  dotCtx.textAlign = "left";
  dotCtx.textBaseline = "middle";
  dotCtx.fillText("Seurat Clusters", left, titleY);

  dotCtx.strokeStyle = "#d7dde5";
  dotCtx.lineWidth = 1;
  dotCtx.beginPath();
  dotCtx.moveTo(left, centerY);
  dotCtx.lineTo(right, centerY);
  dotCtx.stroke();

  dotCtx.font = "700 14px Arial, \"Noto Sans SC\", sans-serif";
  dotCtx.textAlign = "right";
  dotCtx.textBaseline = "middle";
  dotCtx.fillStyle = "#253044";
  dotCtx.fillText(ellipsizeText(dotCtx, payload.gene || state.currentGene || "-", left - 28), left - 14, centerY);

  for (const [index, cluster] of clusters.entries()) {
    const x = left + band * (index + 0.5);
    const avgExpr = Number(cluster.avgExpr) || 0;
    const avgExprScaled = Number(cluster.avgExprScaled);
    const pctExpr = Number(cluster.pctExpr) || 0;
    const radius = dotRadius(pctExpr, pctRange);
    const rgb = Number.isFinite(avgExprScaled)
      ? scaledColorForValue(avgExprScaled, colorRange)
      : colorForValue(avgExpr, colorRange);

    dotCtx.beginPath();
    dotCtx.arc(x, centerY, radius, 0, Math.PI * 2);
    dotCtx.fillStyle = rgbCss(rgb);
    dotCtx.fill();
    dotCtx.strokeStyle = "#7a1c16";
    dotCtx.lineWidth = 0.8;
    dotCtx.stroke();

    if (!compact) {
      dotCtx.fillStyle = "#415064";
      dotCtx.font = "13px Arial, \"Noto Sans SC\", sans-serif";
      dotCtx.textAlign = "center";
      dotCtx.textBaseline = "top";
      dotCtx.fillText(String(cluster.label || cluster.id), x, centerY + 22);
    }
  }

  if (compact) {
    for (const [index, cluster] of clusters.entries()) {
      const x = left + band * (index + 0.5);
      dotCtx.save();
      dotCtx.translate(x, centerY + 30);
      dotCtx.rotate(-Math.PI / 4);
      dotCtx.fillStyle = "#415064";
      dotCtx.font = "12px Arial, \"Noto Sans SC\", sans-serif";
      dotCtx.textAlign = "right";
      dotCtx.textBaseline = "middle";
      dotCtx.fillText(String(cluster.label || cluster.id), 0, 0);
      dotCtx.restore();
    }
  }

  if (!compact) {
    drawDotplotLegend(colorRange, pctRange, rect.width - 26, centerY);
  }
}

async function loadSpatial() {
  if (typeof Path2D === "undefined") {
    throw new Error("当前浏览器不支持 Path2D，无法绘制细胞轮廓");
  }

  setStatus("加载 S1/S2 轮廓和重复索引...");
  const [s1, s2, replicates] = await Promise.all([
    fetch("/data/contours/S1/manifest.json").then((response) => {
      if (!response.ok) throw new Error("缺少 S1 轮廓数据，请先运行 web_viewer/export_contours.py");
      return response.json();
    }),
    fetch("/data/contours/S2/manifest.json").then((response) => {
      if (!response.ok) throw new Error("缺少 S2 轮廓数据，请先运行 web_viewer/export_contours.py");
      return response.json();
    }),
    fetch("/api/replicates").then((response) => {
      if (!response.ok) throw new Error("缺少重复信息，请先运行 web_viewer/export_replicates.py");
      return response.json();
    }),
  ]);

  state.samples.S1 = prepareSpatial(s1, replicates.samples.S1);
  state.samples.S2 = prepareSpatial(s2, replicates.samples.S2);
  fitView();
  updateStats();
  requestDraw();
}

async function loadGenes() {
  const response = await fetch("/api/genes");
  if (!response.ok) return;
  const payload = await response.json();
  els.geneList.innerHTML = "";
  for (const gene of payload.genes.slice(0, 2000)) {
    const option = document.createElement("option");
    option.value = gene;
    els.geneList.appendChild(option);
  }
}

function unpackExpression(samplePayload) {
  const values = samplePayload.values || [];
  const map = new Map();
  let min = Infinity;
  let max = 0;
  for (const [cellId, value] of values) {
    const oldValue = map.get(cellId);
    map.set(cellId, oldValue === undefined ? value : Math.max(oldValue, value));
    if (value > 0) {
      min = Math.min(min, value);
      max = Math.max(max, value);
    }
  }
  return {
    map,
    min: Number.isFinite(min) ? min : 0,
    max,
    nonzero: samplePayload.nonzero || values.length,
  };
}

async function queryGene(gene) {
  gene = gene.trim();
  if (!gene) return;

  setStatus(`查询 ${gene} 表达量和 cluster dotplot...`);
  state.dotplot = { payload: null, error: "", loading: true };
  requestDotplotDraw();

  const encodedGene = encodeURIComponent(gene);
  const [geneResponse, dotplotResponse] = await Promise.all([
    fetch(`/api/gene?gene=${encodedGene}`),
    fetch(`/api/dotplot?gene=${encodedGene}`),
  ]);
  const payload = await geneResponse.json();
  const dotplotPayload = await dotplotResponse.json();

  state.dotplot.loading = false;
  if (!geneResponse.ok) {
    state.dotplot = { payload: null, error: "等待有效基因", loading: false };
    requestDotplotDraw();
    setStatus(payload.error || "查询失败");
    return;
  }

  state.currentGene = gene;
  state.expressionRange = payload.range || { vmin: 0, vmax: 0 };
  state.expressions.S1 = unpackExpression(payload.samples.S1);
  state.expressions.S2 = unpackExpression(payload.samples.S2);
  if (dotplotResponse.ok) {
    state.dotplot = { payload: dotplotPayload, error: "", loading: false };
  } else {
    state.dotplot = {
      payload: null,
      error: dotplotPayload.error || "dotplot 数据不可用",
      loading: false,
    };
  }
  updateStats();
  requestDraw();
  requestDotplotDraw();
}

function updateStats() {
  const spatial = currentSpatial();
  const expression = currentExpression();
  const { vmin, vmax } = state.expressionRange;
  els.currentGene.textContent = state.currentGene || "-";
  els.currentSample.textContent = spatial ? `${state.currentSample} (${spatial.panels.length} reps)` : state.currentSample;
  els.cellCount.textContent = spatial ? spatial.assignedCellCount.toLocaleString() : "-";
  els.nonzeroCount.textContent = expression.nonzero ? expression.nonzero.toLocaleString() : "0";
  els.rangeText.textContent = vmax ? `${formatNumber(vmin)} - ${formatNumber(vmax)}` : "-";
  els.legendMax.textContent = vmax ? formatNumber(vmax) : "max";
  els.legendMin.textContent = Number.isFinite(vmin) ? formatNumber(vmin) : "0";
  updateScaleText();
}

function setSample(sample) {
  state.currentSample = sample;
  document.querySelectorAll(".sample-toggle button").forEach((button) => {
    button.classList.toggle("active", button.dataset.sample === sample);
  });
  fitView();
  updateStats();
  requestDraw();
}

function zoomAt(factor, centerX = canvas.clientWidth / 2, centerY = canvas.clientHeight / 2) {
  const beforeX = (centerX - state.view.x) / state.view.scale;
  const beforeY = (centerY - state.view.y) / state.view.scale;
  state.view.scale = Math.max(MIN_SCALE, Math.min(MAX_SCALE, state.view.scale * factor));
  state.view.x = centerX - beforeX * state.view.scale;
  state.view.y = centerY - beforeY * state.view.scale;
  updateScaleText();
  requestDraw();
}

els.form.addEventListener("submit", (event) => {
  event.preventDefault();
  queryGene(els.input.value);
});

document.querySelectorAll(".sample-toggle button").forEach((button) => {
  button.addEventListener("click", () => setSample(button.dataset.sample));
});

document.getElementById("zoomIn").addEventListener("click", () => zoomAt(1.25));
document.getElementById("zoomOut").addEventListener("click", () => zoomAt(0.8));
document.getElementById("resetView").addEventListener("click", () => {
  fitView();
  requestDraw();
});

canvas.addEventListener("wheel", (event) => {
  event.preventDefault();
  const rect = canvas.getBoundingClientRect();
  const factor = event.deltaY < 0 ? 1.15 : 0.87;
  zoomAt(factor, event.clientX - rect.left, event.clientY - rect.top);
}, { passive: false });

canvas.addEventListener("pointerdown", (event) => {
  state.dragging = true;
  state.lastPointer = { x: event.clientX, y: event.clientY };
  canvas.setPointerCapture(event.pointerId);
});

canvas.addEventListener("pointermove", (event) => {
  if (!state.dragging || !state.lastPointer) return;
  state.view.x += event.clientX - state.lastPointer.x;
  state.view.y += event.clientY - state.lastPointer.y;
  state.lastPointer = { x: event.clientX, y: event.clientY };
  requestDraw();
});

canvas.addEventListener("pointerup", (event) => {
  state.dragging = false;
  state.lastPointer = null;
  canvas.releasePointerCapture(event.pointerId);
});

window.addEventListener("resize", () => {
  resizeCanvas();
  fitView();
  requestDraw();
  requestDotplotDraw();
});

(async function init() {
  try {
    resizeCanvas();
    await loadSpatial();
    await loadGenes();
    await queryGene(els.input.value);
  } catch (error) {
    setStatus(error.message);
    console.error(error);
  }
})();
