/* ============================================================
   NumberGraph — Grapher v2.0
   Canvas rendering: Cartesian, Polar, Parametric, Complex, Implicit
   ============================================================ */

const Grapher = (() => {
    'use strict';

    let canvas, ctx;
    let dpr = 1;

    // View state
    let viewX = -10, viewY = -8, viewW = 20, viewH = 16;

    // Pan state
    let isPanning = false;
    let panStart = { x: 0, y: 0 }, panViewStart = { x: 0, y: 0 };

    // Settings
    let settings = {
        showGrid: true,
        showLabels: true,
        showMinorGrid: false,
        lineWidth: 2.5,
        pointSize: 4,
        density: 2,
    };

    // Graph mode: 'cartesian', 'polar', 'complex', 'parametric'
    let graphMode = 'cartesian';

    let expressions = [];
    let mousePos = null;
    let onCoordsUpdate = null;
    let onPointHover = null;

    // Complex plane cache
    let complexImageData = null;
    let complexCacheBounds = null;

    // ========================================================
    // INIT
    // ========================================================
    function init(canvasElement) {
        canvas = canvasElement;
        ctx = canvas.getContext('2d', { willReadFrequently: true });
        dpr = window.devicePixelRatio || 1;
        resize();
        setupEventListeners();
        render();
    }

    function resize() {
        const rect = canvas.parentElement.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        canvas.style.width = rect.width + 'px';
        canvas.style.height = rect.height + 'px';
        const aspect = rect.width / rect.height;
        viewH = viewW / aspect;
        complexImageData = null; // invalidate cache
    }

    // ========================================================
    // COORDINATE TRANSFORMS
    // ========================================================
    function mathToPixelX(mx) { return ((mx - viewX) / viewW) * canvas.width; }
    function mathToPixelY(my) { return canvas.height - ((my - viewY) / viewH) * canvas.height; }
    function pixelToMathX(px) { return viewX + (px / canvas.width) * viewW; }
    function pixelToMathY(py) { return viewY + ((canvas.height - py) / canvas.height) * viewH; }

    // ========================================================
    // GRID
    // ========================================================
    function getTickSpacing(range) {
        const rough = range / 10;
        const mag = Math.pow(10, Math.floor(Math.log10(rough)));
        const norm = rough / mag;
        let spacing;
        if (norm < 1.5) spacing = 1;
        else if (norm < 3.5) spacing = 2;
        else if (norm < 7.5) spacing = 5;
        else spacing = 10;
        return spacing * mag;
    }

    function drawGrid() {
        if (!settings.showGrid) return;
        const theme = getTheme();
        const tickX = getTickSpacing(viewW);
        const tickY = getTickSpacing(viewH);

        if (settings.showMinorGrid) {
            ctx.strokeStyle = theme.gridMinor;
            ctx.lineWidth = 0.5 * dpr;
            const minorX = tickX / 5, minorY = tickY / 5;
            for (let x = Math.floor(viewX / minorX) * minorX; x <= viewX + viewW; x += minorX) {
                const px = mathToPixelX(x);
                ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, canvas.height); ctx.stroke();
            }
            for (let y = Math.floor(viewY / minorY) * minorY; y <= viewY + viewH; y += minorY) {
                const py = mathToPixelY(y);
                ctx.beginPath(); ctx.moveTo(0, py); ctx.lineTo(canvas.width, py); ctx.stroke();
            }
        }

        // Major grid
        ctx.strokeStyle = theme.gridMajor;
        ctx.lineWidth = 0.8 * dpr;
        const startX = Math.floor(viewX / tickX) * tickX;
        const endX = viewX + viewW;
        for (let x = startX; x <= endX; x += tickX) {
            if (Math.abs(x) < tickX * 0.01) continue;
            const px = mathToPixelX(x);
            ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, canvas.height); ctx.stroke();
        }
        const startY = Math.floor(viewY / tickY) * tickY;
        const endY = viewY + viewH;
        for (let y = startY; y <= endY; y += tickY) {
            if (Math.abs(y) < tickY * 0.01) continue;
            const py = mathToPixelY(y);
            ctx.beginPath(); ctx.moveTo(0, py); ctx.lineTo(canvas.width, py); ctx.stroke();
        }

        // Axes
        ctx.strokeStyle = theme.axis;
        ctx.lineWidth = 1.5 * dpr;
        if (viewY <= 0 && viewY + viewH >= 0) {
            const py = mathToPixelY(0);
            ctx.beginPath(); ctx.moveTo(0, py); ctx.lineTo(canvas.width, py); ctx.stroke();
        }
        if (viewX <= 0 && viewX + viewW >= 0) {
            const px = mathToPixelX(0);
            ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, canvas.height); ctx.stroke();
        }

        // Labels
        if (settings.showLabels) {
            ctx.fillStyle = theme.axisLabel;
            ctx.font = `${11 * dpr}px 'Inter', sans-serif`;
            const axisYPx = mathToPixelY(0);
            const axisXPx = mathToPixelX(0);
            const oY = 6 * dpr, oX = 6 * dpr;

            ctx.textAlign = 'center';
            ctx.textBaseline = 'top';
            for (let x = startX; x <= endX; x += tickX) {
                if (Math.abs(x) < tickX * 0.01) continue;
                const px = mathToPixelX(x);
                const ly = (viewY <= 0 && viewY + viewH >= 0)
                    ? Math.min(Math.max(axisYPx + oY, 0), canvas.height - 20 * dpr)
                    : canvas.height - 20 * dpr;
                ctx.fillText(formatLabel(x, tickX), px, ly);
            }
            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            for (let y = startY; y <= endY; y += tickY) {
                if (Math.abs(y) < tickY * 0.01) continue;
                const py = mathToPixelY(y);
                const lx = (viewX <= 0 && viewX + viewW >= 0) ? Math.max(axisXPx - oX, 35 * dpr) : 35 * dpr;
                ctx.fillText(formatLabel(y, tickY), lx, py);
            }
            if (viewX <= 0 && viewX + viewW >= 0 && viewY <= 0 && viewY + viewH >= 0) {
                ctx.textAlign = 'right'; ctx.textBaseline = 'top';
                ctx.fillText('0', axisXPx - oX, axisYPx + oY);
            }
        }
    }

    function drawPolarGrid() {
        if (!settings.showGrid) return;
        const theme = getTheme();
        const cx = mathToPixelX(0), cy = mathToPixelY(0);
        const maxR = Math.max(Math.abs(viewX), Math.abs(viewX + viewW), Math.abs(viewY), Math.abs(viewY + viewH));
        const tickR = getTickSpacing(maxR);

        // Concentric circles
        ctx.strokeStyle = theme.gridMajor;
        ctx.lineWidth = 0.8 * dpr;
        for (let r = tickR; r <= maxR * 1.5; r += tickR) {
            const rpx = (r / viewW) * canvas.width;
            ctx.beginPath();
            ctx.arc(cx, cy, rpx, 0, Math.PI * 2);
            ctx.stroke();
            if (settings.showLabels) {
                ctx.fillStyle = theme.axisLabel;
                ctx.font = `${11 * dpr}px 'Inter', sans-serif`;
                ctx.textAlign = 'left';
                ctx.textBaseline = 'bottom';
                ctx.fillText(formatLabel(r, tickR), cx + rpx + 4 * dpr, cy - 4 * dpr);
            }
        }

        // Radial lines every 30°
        ctx.strokeStyle = theme.gridMajor;
        ctx.lineWidth = 0.5 * dpr;
        const maxPx = Math.max(canvas.width, canvas.height) * 1.5;
        for (let a = 0; a < 360; a += 30) {
            const rad = a * Math.PI / 180;
            ctx.beginPath();
            ctx.moveTo(cx, cy);
            ctx.lineTo(cx + Math.cos(rad) * maxPx, cy - Math.sin(rad) * maxPx);
            ctx.stroke();
        }

        // Axes (thicker)
        ctx.strokeStyle = theme.axis;
        ctx.lineWidth = 1.5 * dpr;
        ctx.beginPath(); ctx.moveTo(0, cy); ctx.lineTo(canvas.width, cy); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(cx, 0); ctx.lineTo(cx, canvas.height); ctx.stroke();
    }

    function formatLabel(value, spacing) {
        const absVal = Math.abs(value);
        if (spacing >= 1 && absVal === Math.round(absVal)) return Math.round(value).toString();
        const decimals = Math.max(0, Math.ceil(-Math.log10(spacing)) + 1);
        return value.toFixed(Math.min(decimals, 6));
    }

    // ========================================================
    // CARTESIAN GRAPH RENDERING
    // ========================================================
    function drawExpressions() {
        for (const expr of expressions) {
            if (!expr.visible || !expr.parsed) continue;
            const p = expr.parsed;
            if (p.type === 'error' || p.type === 'slider') continue;

            if (p.type === 'polar') {
                drawPolarGraph(p.ast, expr.color, expr.sliderVars || {});
            } else if (p.type === 'parametric') {
                drawParametricGraph(p.astX, p.astY, expr.color, expr.sliderVars || {});
            } else if (p.type === 'implicit') {
                drawImplicitGraph(p.ast, expr.color, expr.sliderVars || {});
            } else if (p.type === 'function') {
                const isDiscrete = expr.plotMode === 'discrete' || (expr.plotMode === 'auto' && p.isDiscrete);
                const isBar = expr.plotMode === 'bar';
                if (isDiscrete || isBar) {
                    drawDiscreteGraph(p.ast, expr.color, expr.sliderVars || {}, isBar);
                } else {
                    drawContinuousGraph(p.ast, expr.color, expr.sliderVars || {});
                }
            }
        }
    }

    function drawContinuousGraph(ast, color, sliderVars) {
        const width = canvas.width;
        const samplesPerPixel = settings.density;
        const totalSamples = Math.floor(width / dpr) * samplesPerPixel;
        ctx.strokeStyle = color;
        ctx.lineWidth = settings.lineWidth * dpr;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';

        const dx = viewW / totalSamples;
        let prevY = null;
        let drawing = false;
        ctx.beginPath();

        for (let i = 0; i <= totalSamples; i++) {
            const mx = viewX + i * dx;
            const vars = { x: mx, n: mx, t: mx, ...sliderVars };
            let my;
            try { my = MathEngine.evaluate(ast, vars); } catch { my = NaN; }

            if (!isFinite(my) || Math.abs(my) > 1e15) {
                if (drawing) { ctx.stroke(); ctx.beginPath(); drawing = false; }
                prevY = null;
                continue;
            }

            const px = mathToPixelX(mx);
            const py = mathToPixelY(my);

            if (prevY !== null && Math.abs(py - prevY) > canvas.height * 0.8) {
                ctx.stroke(); ctx.beginPath(); ctx.moveTo(px, py); drawing = true;
            } else if (!drawing) {
                ctx.moveTo(px, py); drawing = true;
            } else {
                ctx.lineTo(px, py);
            }
            prevY = py;
        }
        if (drawing) ctx.stroke();
    }

    function drawDiscreteGraph(ast, color, sliderVars, isBar) {
        const xStart = Math.max(1, Math.floor(viewX));
        const xEnd = Math.min(1100000, Math.ceil(viewX + viewW));

        if (isBar) {
            const barPixelWidth = Math.max(1, (mathToPixelX(1) - mathToPixelX(0)) * 0.7);
            ctx.fillStyle = color + '88';
            ctx.strokeStyle = color;
            ctx.lineWidth = 1 * dpr;
            for (let n = xStart; n <= xEnd; n++) {
                const vars = { x: n, n, ...sliderVars };
                let y;
                try { y = MathEngine.evaluate(ast, vars); } catch { continue; }
                if (!isFinite(y)) continue;
                const px = mathToPixelX(n), py = mathToPixelY(y), py0 = mathToPixelY(0);
                ctx.fillRect(px - barPixelWidth / 2, Math.min(py, py0), barPixelWidth, Math.abs(py0 - py));
                ctx.strokeRect(px - barPixelWidth / 2, Math.min(py, py0), barPixelWidth, Math.abs(py0 - py));
            }
        } else {
            const ptRadius = settings.pointSize * dpr;
            ctx.fillStyle = color;
            ctx.strokeStyle = color + '40';
            ctx.lineWidth = 1 * dpr;
            ctx.beginPath();
            let first = true;
            for (let n = xStart; n <= xEnd; n++) {
                const vars = { x: n, n, ...sliderVars };
                let y;
                try { y = MathEngine.evaluate(ast, vars); } catch { first = true; continue; }
                if (!isFinite(y)) { first = true; continue; }
                const px = mathToPixelX(n), py = mathToPixelY(y);
                if (first) { ctx.moveTo(px, py); first = false; } else { ctx.lineTo(px, py); }
            }
            ctx.stroke();
            for (let n = xStart; n <= xEnd; n++) {
                const vars = { x: n, n, ...sliderVars };
                let y;
                try { y = MathEngine.evaluate(ast, vars); } catch { continue; }
                if (!isFinite(y)) continue;
                const px = mathToPixelX(n), py = mathToPixelY(y);
                ctx.beginPath(); ctx.arc(px, py, ptRadius, 0, Math.PI * 2); ctx.fill();
            }
        }
    }

    // ========================================================
    // POLAR GRAPH RENDERING
    // ========================================================
    function drawPolarGraph(ast, color, sliderVars) {
        ctx.strokeStyle = color;
        ctx.lineWidth = settings.lineWidth * dpr;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';

        const totalSteps = 2000 * settings.density;
        const thetaMax = 4 * Math.PI;
        const dTheta = thetaMax / totalSteps;
        let drawing = false;
        ctx.beginPath();

        for (let i = 0; i <= totalSteps; i++) {
            const theta = i * dTheta;
            const vars = { theta, t: theta, x: theta, n: theta, ...sliderVars };
            let r;
            try { r = MathEngine.evaluate(ast, vars); } catch { r = NaN; }
            if (!isFinite(r)) {
                if (drawing) { ctx.stroke(); ctx.beginPath(); drawing = false; }
                continue;
            }
            const mx = r * Math.cos(theta);
            const my = r * Math.sin(theta);
            const px = mathToPixelX(mx);
            const py = mathToPixelY(my);
            if (!drawing) { ctx.moveTo(px, py); drawing = true; } else { ctx.lineTo(px, py); }
        }
        if (drawing) ctx.stroke();
    }

    // ========================================================
    // PARAMETRIC GRAPH RENDERING
    // ========================================================
    function drawParametricGraph(astX, astY, color, sliderVars) {
        ctx.strokeStyle = color;
        ctx.lineWidth = settings.lineWidth * dpr;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';

        const totalSteps = 2000 * settings.density;
        const tMin = -10, tMax = 10;
        const dt = (tMax - tMin) / totalSteps;
        let drawing = false;
        ctx.beginPath();

        for (let i = 0; i <= totalSteps; i++) {
            const t = tMin + i * dt;
            const vars = { t, x: t, n: t, theta: t, ...sliderVars };
            let mx, my;
            try { mx = MathEngine.evaluate(astX, vars); my = MathEngine.evaluate(astY, vars); } catch { mx = NaN; my = NaN; }
            if (!isFinite(mx) || !isFinite(my)) {
                if (drawing) { ctx.stroke(); ctx.beginPath(); drawing = false; }
                continue;
            }
            const px = mathToPixelX(mx);
            const py = mathToPixelY(my);
            if (!drawing) { ctx.moveTo(px, py); drawing = true; } else { ctx.lineTo(px, py); }
        }
        if (drawing) ctx.stroke();
    }

    // ========================================================
    // IMPLICIT CURVE RENDERING (Marching Squares)
    // ========================================================
    function drawImplicitGraph(ast, color, sliderVars) {
        ctx.strokeStyle = color;
        ctx.lineWidth = settings.lineWidth * dpr;

        const resolution = Math.max(1, Math.floor(3 / settings.density));
        const w = Math.floor(canvas.width / dpr / resolution);
        const h = Math.floor(canvas.height / dpr / resolution);

        // Evaluate on grid
        const grid = new Float64Array((w + 1) * (h + 1));
        for (let j = 0; j <= h; j++) {
            for (let ii = 0; ii <= w; ii++) {
                const mx = pixelToMathX(ii * resolution * dpr);
                const my = pixelToMathY(j * resolution * dpr);
                const vars = { x: mx, y: my, ...sliderVars };
                try { grid[j * (w + 1) + ii] = MathEngine.evaluate(ast, vars); } catch { grid[j * (w + 1) + ii] = NaN; }
            }
        }

        // March
        ctx.beginPath();
        for (let j = 0; j < h; j++) {
            for (let ii = 0; ii < w; ii++) {
                const v00 = grid[j * (w + 1) + ii];
                const v10 = grid[j * (w + 1) + ii + 1];
                const v01 = grid[(j + 1) * (w + 1) + ii];
                const v11 = grid[(j + 1) * (w + 1) + ii + 1];

                if (!isFinite(v00) || !isFinite(v10) || !isFinite(v01) || !isFinite(v11)) continue;

                const s00 = v00 > 0 ? 1 : 0;
                const s10 = v10 > 0 ? 1 : 0;
                const s01 = v01 > 0 ? 1 : 0;
                const s11 = v11 > 0 ? 1 : 0;
                const cellCase = s00 | (s10 << 1) | (s01 << 2) | (s11 << 3);

                if (cellCase === 0 || cellCase === 15) continue;

                const x0 = ii * resolution * dpr;
                const y0 = j * resolution * dpr;
                const sz = resolution * dpr;

                // Interpolation helpers
                const lerpX = (va, vb) => {
                    if (Math.abs(va - vb) < 1e-20) return 0.5;
                    return Math.max(0, Math.min(1, -va / (vb - va)));
                };

                const top = lerpX(v00, v10);
                const bottom = lerpX(v01, v11);
                const left = lerpX(v00, v01);
                const right = lerpX(v10, v11);

                const segments = marchingSquaresSegments(cellCase, x0, y0, sz, top, bottom, left, right);
                for (const seg of segments) {
                    ctx.moveTo(seg[0], seg[1]);
                    ctx.lineTo(seg[2], seg[3]);
                }
            }
        }
        ctx.stroke();
    }

    function marchingSquaresSegments(c, x, y, s, top, bot, left, right) {
        const segs = [];
        const t = [x + top * s, y];
        const b = [x + bot * s, y + s];
        const l = [x, y + left * s];
        const r = [x + s, y + right * s];

        const lineSegs = {
            1: [[l, t]], 2: [[t, r]], 3: [[l, r]], 4: [[l, b]],
            5: [[l, t], [b, r]], 6: [[t, b]], 7: [[l, b]], 8: [[b, r]],
            9: [[t, b]], 10: [[t, l], [b, r]], 11: [[b, r]], 12: [[l, r]],
            13: [[t, r]], 14: [[l, t]],
        };
        if (lineSegs[c]) {
            for (const [p1, p2] of lineSegs[c]) {
                segs.push([p1[0], p1[1], p2[0], p2[1]]);
            }
        }
        return segs;
    }

    // ========================================================
    // COMPLEX PLANE RENDERING (Domain Coloring)
    // ========================================================
    function drawComplexPlane() {
        const visibleComplex = expressions.filter(e =>
            e.visible && e.parsed && e.parsed.type === 'function' && !e.parsed.isDiscrete
        );
        if (visibleComplex.length === 0) return;

        // Use first expression for domain coloring
        const expr = visibleComplex[0];
        const ast = expr.parsed.ast;

        // Check cache validity
        const cacheKey = `${viewX},${viewY},${viewW},${viewH},${canvas.width},${canvas.height}`;
        if (complexImageData && complexCacheBounds === cacheKey) {
            ctx.putImageData(complexImageData, 0, 0);
            return;
        }

        const w = canvas.width;
        const h = canvas.height;
        const step = Math.max(1, Math.floor(4 / settings.density));
        const imgData = ctx.createImageData(w, h);
        const data = imgData.data;

        for (let py = 0; py < h; py += step) {
            for (let px = 0; px < w; px += step) {
                const mx = pixelToMathX(px);
                const my = pixelToMathY(py);
                const z = new MathEngine.Complex(mx, my);

                let result;
                try {
                    result = MathEngine.evaluateComplex(ast, z);
                } catch {
                    result = new MathEngine.Complex(NaN, NaN);
                }

                let r = 0, g = 0, b = 0;
                if (result.isFinite()) {
                    const hue = (result.arg() / (2 * Math.PI) + 1) % 1;
                    const mag = result.abs();
                    const lightness = 1 - 1 / (1 + mag * 0.3);
                    const sat = 0.8;
                    [r, g, b] = hslToRgb(hue, sat, lightness);

                    // Add contour lines for magnitude
                    const logMag = Math.log2(mag + 1e-10);
                    const frac = logMag - Math.floor(logMag);
                    if (frac < 0.05 || frac > 0.95) {
                        r = Math.floor(r * 0.7);
                        g = Math.floor(g * 0.7);
                        b = Math.floor(b * 0.7);
                    }
                    // Phase contour lines
                    const phaseFrac = ((result.arg() / (Math.PI / 6)) % 1 + 1) % 1;
                    if (phaseFrac < 0.03 || phaseFrac > 0.97) {
                        r = Math.floor(r * 0.75);
                        g = Math.floor(g * 0.75);
                        b = Math.floor(b * 0.75);
                    }
                }

                // Fill step x step block
                for (let dy = 0; dy < step && py + dy < h; dy++) {
                    for (let dx = 0; dx < step && px + dx < w; dx++) {
                        const idx = ((py + dy) * w + (px + dx)) * 4;
                        data[idx] = r;
                        data[idx + 1] = g;
                        data[idx + 2] = b;
                        data[idx + 3] = 255;
                    }
                }
            }
        }

        ctx.putImageData(imgData, 0, 0);
        complexImageData = imgData;
        complexCacheBounds = cacheKey;
    }

    function hslToRgb(h, s, l) {
        let r, g, b;
        if (s === 0) {
            r = g = b = l;
        } else {
            const hue2rgb = (p, q, t) => {
                if (t < 0) t += 1;
                if (t > 1) t -= 1;
                if (t < 1/6) return p + (q - p) * 6 * t;
                if (t < 1/2) return q;
                if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
                return p;
            };
            const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
            const p = 2 * l - q;
            r = hue2rgb(p, q, h + 1/3);
            g = hue2rgb(p, q, h);
            b = hue2rgb(p, q, h - 1/3);
        }
        return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
    }

    // ========================================================
    // CROSSHAIR / HOVER
    // ========================================================
    function drawCrosshair() {
        if (!mousePos) return;
        const theme = getTheme();
        const mx = pixelToMathX(mousePos.x * dpr);
        const my = pixelToMathY(mousePos.y * dpr);

        ctx.strokeStyle = theme.crosshair;
        ctx.lineWidth = 0.8 * dpr;
        ctx.setLineDash([4 * dpr, 4 * dpr]);
        ctx.beginPath(); ctx.moveTo(mousePos.x * dpr, 0); ctx.lineTo(mousePos.x * dpr, canvas.height); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(0, mousePos.y * dpr); ctx.lineTo(canvas.width, mousePos.y * dpr); ctx.stroke();
        ctx.setLineDash([]);

        // Snap to nearest point
        let nearestInfo = null, nearestDist = Infinity;
        for (const expr of expressions) {
            if (!expr.visible || !expr.parsed || expr.parsed.type !== 'function') continue;
            const isDiscrete = expr.plotMode === 'discrete' || (expr.plotMode === 'auto' && expr.parsed.isDiscrete);
            if (isDiscrete) {
                const n = Math.round(mx);
                if (n >= 1) {
                    const vars = { x: n, n, ...(expr.sliderVars || {}) };
                    try {
                        const y = MathEngine.evaluate(expr.parsed.ast, vars);
                        if (isFinite(y)) {
                            const px = mathToPixelX(n), py = mathToPixelY(y);
                            const dist = Math.hypot(mousePos.x * dpr - px, mousePos.y * dpr - py);
                            if (dist < nearestDist && dist < 30 * dpr) {
                                nearestDist = dist;
                                nearestInfo = { x: n, y, px, py, color: expr.color, expr: expr.raw };
                            }
                        }
                    } catch {}
                }
            } else {
                const vars = { x: mx, n: mx, t: mx, ...(expr.sliderVars || {}) };
                try {
                    const y = MathEngine.evaluate(expr.parsed.ast, vars);
                    if (isFinite(y)) {
                        const py = mathToPixelY(y);
                        const dist = Math.abs(mousePos.y * dpr - py);
                        if (dist < nearestDist && dist < 20 * dpr) {
                            nearestDist = dist;
                            nearestInfo = { x: mx, y, px: mousePos.x * dpr, py, color: expr.color, expr: expr.raw };
                        }
                    }
                } catch {}
            }
        }

        if (nearestInfo) {
            ctx.fillStyle = nearestInfo.color;
            ctx.beginPath(); ctx.arc(nearestInfo.px, nearestInfo.py, 5 * dpr, 0, Math.PI * 2); ctx.fill();
            ctx.strokeStyle = getTheme().bg;
            ctx.lineWidth = 2 * dpr;
            ctx.beginPath(); ctx.arc(nearestInfo.px, nearestInfo.py, 5 * dpr, 0, Math.PI * 2); ctx.stroke();
            if (onPointHover) onPointHover({
                x: nearestInfo.x, y: nearestInfo.y,
                screenX: nearestInfo.px / dpr, screenY: nearestInfo.py / dpr,
                expr: nearestInfo.expr,
            });
        } else {
            if (onPointHover) onPointHover(null);
        }

        if (onCoordsUpdate) {
            if (graphMode === 'complex') {
                onCoordsUpdate(mx, my, true); // complex mode flag
            } else {
                onCoordsUpdate(mx, my, false);
            }
        }
    }

    // ========================================================
    // MAIN RENDER
    // ========================================================
    function render() {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        const theme = getTheme();
        ctx.fillStyle = theme.bg;
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        if (graphMode === 'complex') {
            drawComplexPlane();
            // Draw grid overlay lightly
            const oldAlpha = ctx.globalAlpha;
            ctx.globalAlpha = 0.3;
            drawGrid();
            ctx.globalAlpha = oldAlpha;
        } else if (graphMode === 'polar') {
            drawPolarGrid();
            drawExpressions();
        } else {
            drawGrid();
            drawExpressions();
        }

        drawCrosshair();
    }

    // ========================================================
    // THEME
    // ========================================================
    function getTheme() {
        const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
        return isDark ? {
            bg: '#0d1117',
            gridMajor: 'rgba(255,255,255,0.07)',
            gridMinor: 'rgba(255,255,255,0.025)',
            axis: 'rgba(255,255,255,0.25)',
            axisLabel: 'rgba(255,255,255,0.45)',
            crosshair: 'rgba(255,255,255,0.15)',
        } : {
            bg: '#ffffff',
            gridMajor: 'rgba(0,0,0,0.06)',
            gridMinor: 'rgba(0,0,0,0.025)',
            axis: 'rgba(0,0,0,0.2)',
            axisLabel: 'rgba(0,0,0,0.5)',
            crosshair: 'rgba(0,0,0,0.12)',
        };
    }

    // ========================================================
    // EVENTS
    // ========================================================
    function setupEventListeners() {
        canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const rect = canvas.getBoundingClientRect();
            const mx = pixelToMathX((e.clientX - rect.left) * dpr);
            const my = pixelToMathY((e.clientY - rect.top) * dpr);
            const factor = e.deltaY > 0 ? 1.15 : 1 / 1.15;
            const newW = viewW * factor, newH = viewH * factor;
            viewX = mx - (mx - viewX) * (newW / viewW);
            viewY = my - (my - viewY) * (newH / viewH);
            viewW = newW; viewH = newH;
            complexImageData = null;
            render();
        }, { passive: false });

        canvas.addEventListener('mousedown', (e) => {
            if (e.button === 0) {
                isPanning = true;
                panStart = { x: e.clientX, y: e.clientY };
                panViewStart = { x: viewX, y: viewY };
                canvas.classList.add('panning');
            }
        });

        window.addEventListener('mousemove', (e) => {
            if (isPanning) {
                const rect = canvas.getBoundingClientRect();
                viewX = panViewStart.x - (e.clientX - panStart.x) / rect.width * viewW;
                viewY = panViewStart.y + (e.clientY - panStart.y) / rect.height * viewH;
                complexImageData = null;
                render();
            }
            const rect = canvas.getBoundingClientRect();
            if (e.clientX >= rect.left && e.clientX <= rect.right &&
                e.clientY >= rect.top && e.clientY <= rect.bottom) {
                mousePos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
                render();
            }
        });

        window.addEventListener('mouseup', () => { isPanning = false; canvas.classList.remove('panning'); });

        canvas.addEventListener('mouseleave', () => {
            mousePos = null;
            if (onCoordsUpdate) onCoordsUpdate(null, null);
            if (onPointHover) onPointHover(null);
            render();
        });

        // Touch
        let touchDist = 0, touchCenter = { x: 0, y: 0 };
        canvas.addEventListener('touchstart', (e) => {
            e.preventDefault();
            if (e.touches.length === 1) {
                isPanning = true;
                panStart = { x: e.touches[0].clientX, y: e.touches[0].clientY };
                panViewStart = { x: viewX, y: viewY };
            } else if (e.touches.length === 2) {
                isPanning = false;
                touchDist = Math.hypot(e.touches[1].clientX - e.touches[0].clientX, e.touches[1].clientY - e.touches[0].clientY);
                touchCenter = { x: (e.touches[0].clientX + e.touches[1].clientX) / 2, y: (e.touches[0].clientY + e.touches[1].clientY) / 2 };
            }
        }, { passive: false });

        canvas.addEventListener('touchmove', (e) => {
            e.preventDefault();
            if (e.touches.length === 1 && isPanning) {
                const rect = canvas.getBoundingClientRect();
                viewX = panViewStart.x - (e.touches[0].clientX - panStart.x) / rect.width * viewW;
                viewY = panViewStart.y + (e.touches[0].clientY - panStart.y) / rect.height * viewH;
                complexImageData = null;
                render();
            } else if (e.touches.length === 2) {
                const newDist = Math.hypot(e.touches[1].clientX - e.touches[0].clientX, e.touches[1].clientY - e.touches[0].clientY);
                const factor = touchDist / newDist;
                const rect = canvas.getBoundingClientRect();
                const cx = pixelToMathX((touchCenter.x - rect.left) * dpr);
                const cy = pixelToMathY((touchCenter.y - rect.top) * dpr);
                const newW = viewW * factor, newH = viewH * factor;
                viewX = cx - (cx - viewX) * (newW / viewW);
                viewY = cy - (cy - viewY) * (newH / viewH);
                viewW = newW; viewH = newH;
                touchDist = newDist;
                complexImageData = null;
                render();
            }
        }, { passive: false });

        canvas.addEventListener('touchend', () => { isPanning = false; });
        window.addEventListener('resize', () => { resize(); render(); });
    }

    // ========================================================
    // PUBLIC API
    // ========================================================
    return {
        init, render, resize,
        setExpressions(exprs) { expressions = exprs; complexImageData = null; render(); },
        setGraphMode(mode) { graphMode = mode; complexImageData = null; render(); },
        getGraphMode() { return graphMode; },
        getView() {
            return { viewX, viewY, viewW, viewH };
        },
        setView(view) {
            if (!view || typeof view !== 'object') return;
            const nx = Number(view.viewX);
            const ny = Number(view.viewY);
            const nw = Number(view.viewW);
            const nh = Number(view.viewH);
            if (!isFinite(nx) || !isFinite(ny) || !isFinite(nw) || !isFinite(nh)) return;
            if (nw <= 0 || nh <= 0) return;
            viewX = nx;
            viewY = ny;
            viewW = nw;
            viewH = nh;
            complexImageData = null;
            render();
        },
        zoomIn() {
            const cx = viewX + viewW / 2, cy = viewY + viewH / 2;
            viewW /= 1.3; viewH /= 1.3;
            viewX = cx - viewW / 2; viewY = cy - viewH / 2;
            complexImageData = null; render();
        },
        zoomOut() {
            const cx = viewX + viewW / 2, cy = viewY + viewH / 2;
            viewW *= 1.3; viewH *= 1.3;
            viewX = cx - viewW / 2; viewY = cy - viewH / 2;
            complexImageData = null; render();
        },
        resetView() {
            viewX = -10; viewW = 20;
            viewH = viewW / (canvas.width / canvas.height);
            viewY = -viewH / 2;
            complexImageData = null; render();
        },
        updateSettings(s) { Object.assign(settings, s); complexImageData = null; render(); },
        exportPNG() { return canvas.toDataURL('image/png'); },
        setOnCoordsUpdate(fn) { onCoordsUpdate = fn; },
        setOnPointHover(fn) { onPointHover = fn; },
    };
})();
