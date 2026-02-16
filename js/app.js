/* ============================================================
   MathGraph — App Controller v2.0
   Expression management, UI bindings, graph mode switching
   ============================================================ */

(() => {
    'use strict';

    // ========================================================
    // STATE
    // ========================================================
    const COLORS = [
        '#58a6ff', '#f47067', '#a371f7', '#3fb950',
        '#d29922', '#f778ba', '#79c0ff', '#ff9f43',
        '#22d3ee', '#a3e635', '#e879f9', '#fb923c'
    ];
    let colorIdx = 0;

    let expressions = [];
    let sliders = {};
    let nextId = 1;
    const STORAGE_KEY = 'mathgraph.state.v3';
    const POWER_SLOT = '□';
    let saveTimer = null;

    // ========================================================
    // DOM REFS
    // ========================================================
    const $exprList = document.getElementById('expr-list');
    const $btnAdd = document.getElementById('btn-add');
    const $presetsToggle = document.getElementById('presets-toggle');
    const $presetsList = document.getElementById('presets-list');
    const $coordsDisplay = document.getElementById('coords-display');
    const $coordsText = document.getElementById('coords-text');
    const $tooltip = document.getElementById('point-tooltip');
    const $modeIndicator = document.getElementById('mode-indicator');
    const $computeInput = document.getElementById('compute-input');
    const $computeResult = document.getElementById('compute-result');
    const $autocomplete = document.getElementById('autocomplete');
    const $presetSearch = document.getElementById('preset-search');
    const $btnClearAll = document.getElementById('btn-clear-all');
    const $btnExport = document.getElementById('btn-export');
    const $btnSidebarToggle = document.getElementById('btn-sidebar-toggle');

    // ========================================================
    // GRAPH INIT
    // ========================================================
    Grapher.init(document.getElementById('graph-canvas'));

    Grapher.setOnCoordsUpdate((x, y, isComplex) => {
        if (x == null || y == null) {
            $coordsDisplay.classList.add('hidden');
            return;
        }
        $coordsDisplay.classList.remove('hidden');
        if (isComplex) {
            const sign = y >= 0 ? '+' : '';
            $coordsText.textContent = `z = ${x.toFixed(3)} ${sign} ${y.toFixed(3)}i`;
        } else {
            $coordsText.textContent = `(${x.toFixed(4)}, ${y.toFixed(4)})`;
        }
    });

    Grapher.setOnPointHover(info => {
        if (!info) { $tooltip.classList.add('hidden'); return; }
        $tooltip.classList.remove('hidden');
        $tooltip.innerHTML = `<span class="tt-label">${info.expr}</span><br><span class="tt-value">f(${formatNum(info.x)}) = ${formatNum(info.y)}</span>`;
        const rect = document.getElementById('graph-area').getBoundingClientRect();
        let tx = info.screenX + rect.left + 16;
        let ty = info.screenY + rect.top - 10;
        if (tx + 200 > window.innerWidth) tx -= 232;
        if (ty < 10) ty = 10;
        $tooltip.style.left = tx + 'px';
        $tooltip.style.top = ty + 'px';
    });

    function formatNum(v) {
        if (!isFinite(v)) return String(v);
        if (Number.isInteger(v)) return String(v);
        return Math.abs(v) > 1e6 || (Math.abs(v) < 1e-4 && v !== 0) ? v.toExponential(4) : v.toFixed(4);
    }

    function isMathFieldElement(el) {
        return !!el && el.tagName === 'MATH-FIELD';
    }

    function getFieldValue(el) {
        if (!el) return '';
        const value = typeof el.value === 'string' ? el.value : (el.getAttribute('value') || '');
        return value;
    }

    function setFieldValue(el, value) {
        if (!el) return;
        if (typeof el.value === 'string' || isMathFieldElement(el)) {
            el.value = value;
            return;
        }
        el.setAttribute('value', value);
    }

    function configureMathField(el) {
        if (!isMathFieldElement(el)) return;
        try {
            if (typeof el.setOptions === 'function') {
                el.setOptions({
                    smartMode: false,
                    letterShapeStyle: 'tex',
                    virtualKeyboardMode: 'onfocus',
                });
            }
        } catch {}
    }

    function insertTextAtCursor($input, text, selectionOffset = 0, selectionLength = 0) {
        const start = $input.selectionStart ?? $input.value.length;
        const end = $input.selectionEnd ?? start;
        const before = $input.value.slice(0, start);
        const after = $input.value.slice(end);
        $input.value = before + text + after;
        const selStart = start + selectionOffset;
        const selEnd = selStart + selectionLength;
        $input.selectionStart = selStart;
        $input.selectionEnd = selEnd;
        $input.dispatchEvent(new Event('input'));
    }

    function insertExponentSlot($input) {
        const start = $input.selectionStart ?? $input.value.length;
        const end = $input.selectionEnd ?? start;
        const selected = $input.value.slice(start, end);

        if (selected) {
            insertTextAtCursor($input, `^{${selected}}`, selected.length + 3, 0);
            return;
        }

        insertTextAtCursor($input, `^{${POWER_SLOT}}`, 2, 1);
    }

    function insertBracePair($input) {
        insertTextAtCursor($input, '{}', 1, 0);
    }

    function handleLatexInputShortcuts(e, $input) {
        if (isMathFieldElement($input)) return false;
        if (e.ctrlKey || e.altKey || e.metaKey) return false;

        if (e.key === '^') {
            e.preventDefault();
            insertExponentSlot($input);
            return true;
        }

        if (e.key === '{') {
            e.preventDefault();
            insertBracePair($input);
            return true;
        }

        return false;
    }

    // ========================================================
    // EXPRESSION MANAGEMENT
    // ========================================================
    function addExpression(text = '', focusInput = true) {
        const id = nextId++;
        const color = COLORS[colorIdx % COLORS.length];
        colorIdx++;

        const expr = {
            id, text, color,
            visible: true,
            parsed: null,
            plotMode: 'auto',
            sliderVars: {},
        };
        expressions.push(expr);
        renderExprItem(expr);
        if (text) updateExpression(expr);
        if (focusInput) {
            const input = $exprList.querySelector(`[data-id="${id}"] .expr-input`);
            if (input) setTimeout(() => input.focus(), 50);
        }
        syncGrapher();
        return expr;
    }

    function clearExpressions(startWithEmpty = true) {
        expressions = [];
        sliders = {};
        $exprList.innerHTML = '';
        if (startWithEmpty) {
            addExpression('', false);
        } else {
            syncGrapher();
        }
    }

    function renderExprItem(expr) {
        const div = document.createElement('div');
        div.className = 'expr-item';
        div.dataset.id = expr.id;

        div.innerHTML = `
            <div class="expr-color-bar${expr.visible ? '' : ' invisible'}" style="background:${expr.color}">
                <input type="color" class="expr-color-picker" value="${expr.color}" title="Change color">
            </div>
            <div class="expr-body">
                <div class="expr-input-row">
                    <span class="expr-label">f<sub>${expr.id}</sub></span>
                    <math-field class="expr-input expr-math-field" virtual-keyboard-mode="onfocus" smart-mode="false" placeholder="Type LaTeX expression..."></math-field>
                </div>
            </div>
            <div class="expr-actions">
                <button class="expr-action-btn toggle-visibility" title="Toggle Visibility">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"/><circle cx="12" cy="12" r="3"/></svg>
                </button>
                <button class="expr-action-btn duplicate" title="Duplicate">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>
                </button>
                <button class="expr-action-btn delete" title="Delete">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="18" y1="6" x2="6" y2="18"/><line x1="6" y1="6" x2="18" y2="18"/></svg>
                </button>
            </div>
        `;

        const $input = div.querySelector('.expr-input');
        configureMathField($input);
        setFieldValue($input, expr.text);

        $input.addEventListener('input', () => {
            expr.text = getFieldValue($input);
            updateExpression(expr);
            syncGrapher();
            if (!isMathFieldElement($input)) showAutocomplete($input);
        });
        $input.addEventListener('focus', () => {
            div.classList.add('active');
            if (!isMathFieldElement($input)) showAutocomplete($input);
        });
        $input.addEventListener('blur', () => {
            div.classList.remove('active');
            setTimeout(() => hideAutocomplete(), 200);
        });
        $input.addEventListener('keydown', (e) => {
            if (handleLatexInputShortcuts(e, $input)) return;
            if (e.key === 'Enter') { hideAutocomplete(); addExpression(); }
            else if (e.key === 'Escape') { hideAutocomplete(); $input.blur(); }
            else if (e.key === 'Tab') handleAutocompleteTab(e);
            else if (e.key === 'ArrowDown' || e.key === 'ArrowUp') handleAutocompleteNav(e);
        });

        div.querySelector('.toggle-visibility').addEventListener('click', () => {
            expr.visible = !expr.visible;
            div.querySelector('.expr-color-bar').classList.toggle('invisible', !expr.visible);
            syncGrapher();
        });

        div.querySelector('.duplicate').addEventListener('click', () => {
            addExpression(expr.text, false);
        });

        div.querySelector('.delete').addEventListener('click', () => {
            deleteExpression(expr.id);
        });

        div.querySelector('.expr-color-picker').addEventListener('input', (e) => {
            expr.color = e.target.value;
            div.querySelector('.expr-color-bar').style.background = expr.color;
            syncGrapher();
        });

        $exprList.appendChild(div);
    }

    function updateExpression(expr) {
        const parsed = MathEngine.parseExpression(expr.text);
        expr.parsed = parsed;

        const $item = $exprList.querySelector(`[data-id="${expr.id}"]`);
        if (!$item) return;
        const $input = $item.querySelector('.expr-input');
        const $bodyDiv = $item.querySelector('.expr-body');

        // Remove old type badge
        const oldBadge = $bodyDiv.querySelector('.expr-type-badge');
        if (oldBadge) oldBadge.remove();

        // Remove old slider
        const oldSlider = $bodyDiv.querySelector('.expr-slider-row');
        if (oldSlider) oldSlider.remove();

        if (!parsed) {
            $item.classList.remove('has-error');
            return;
        }

        if (parsed.type === 'error') {
            $item.classList.add('has-error');
            return;
        }
        $item.classList.remove('has-error');

        // Type label
        let badgeType = null;
        if (parsed.type === 'polar') badgeType = 'polar';
        else if (parsed.type === 'parametric') badgeType = 'parametric';
        else if (parsed.type === 'implicit') badgeType = 'implicit';
        else if (parsed.type === 'function' && parsed.isComplexHint) badgeType = 'complex';

        if (badgeType) {
            const badge = document.createElement('span');
            badge.className = `expr-type-badge ${badgeType}`;
            badge.textContent = badgeType;
            const row = $bodyDiv.querySelector('.expr-input-row');
            row.insertBefore(badge, $input);
        }

        // Slider handling
        if (parsed.type === 'slider') {
            sliders[parsed.name] = parsed.value;
            updateAllSliderVars();

            const sliderRow = document.createElement('div');
            sliderRow.className = 'expr-slider-row';
            sliderRow.innerHTML = `
                <span style="font-family:var(--font-mono);font-size:12px;color:var(--text-muted)">${parsed.name} =</span>
                <input type="range" class="expr-slider" min="-10" max="10" step="0.01" value="${parsed.value}">
                <span class="expr-slider-value">${Number(parsed.value).toFixed(2)}</span>
            `;
            $bodyDiv.appendChild(sliderRow);

            const slider = sliderRow.querySelector('.expr-slider');
            const valDisplay = sliderRow.querySelector('.expr-slider-value');
            slider.addEventListener('input', () => {
                const v = parseFloat(slider.value);
                sliders[parsed.name] = v;
                valDisplay.textContent = v.toFixed(2);
                updateAllSliderVars();
                syncGrapher();
            });
        }
    }

    function deleteExpression(id) {
        const idx = expressions.findIndex(e => e.id === id);
        if (idx === -1) return;
        const expr = expressions[idx];

        // Clean up slider
        if (expr.parsed && expr.parsed.type === 'slider') {
            delete sliders[expr.parsed.name];
            updateAllSliderVars();
        }

        expressions.splice(idx, 1);
        const $item = $exprList.querySelector(`[data-id="${id}"]`);
        if ($item) $item.remove();
        syncGrapher();
    }

    function updateAllSliderVars() {
        for (const expr of expressions) {
            expr.sliderVars = { ...sliders };
        }
    }

    function syncGrapher() {
        updateAllSliderVars();
        const grapherExprs = expressions.map(e => ({
            visible: e.visible,
            parsed: e.parsed,
            color: e.color,
            plotMode: e.plotMode,
            sliderVars: e.sliderVars,
            raw: e.text,
        }));
        Grapher.setExpressions(grapherExprs);
        queueSave();
    }

    function escHtml(str) {
        return str.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    }

    // ========================================================
    // AUTOCOMPLETE
    // ========================================================
    let acItems = [];
    let acSelected = -1;
    let acTarget = null;

    function showAutocomplete($input) {
        if (isMathFieldElement($input)) {
            hideAutocomplete();
            return;
        }
        acTarget = $input;
        const val = $input.value;
        const cursor = $input.selectionStart || 0;

        // Find the current word being typed
        const before = val.substring(0, cursor);
        const match = before.match(/\\?([a-zA-Z_]\w*)$/);
        if (!match || match[1].length < 1) { hideAutocomplete(); return; }

        const prefix = match[1];
        acItems = MathEngine.getAutocompleteSuggestions(prefix);

        if (acItems.length === 0) { hideAutocomplete(); return; }

        acSelected = -1;
        $autocomplete.classList.remove('hidden');
        $autocomplete.innerHTML = acItems.map((item, i) =>
            `<div class="ac-item" data-index="${i}">
                <span class="ac-name">${item.name}${item.args > 0 ? '(' : ''}</span>
                <span class="ac-desc">${item.desc || ''}</span>
            </div>`
        ).join('');

        // Position
        const rect = $input.getBoundingClientRect();
        $autocomplete.style.left = rect.left + 'px';
        $autocomplete.style.top = (rect.bottom + 4) + 'px';
        $autocomplete.style.minWidth = rect.width + 'px';

        // Click handlers
        $autocomplete.querySelectorAll('.ac-item').forEach(el => {
            el.addEventListener('mousedown', (e) => {
                e.preventDefault();
                const idx = parseInt(el.dataset.index);
                applyAutocomplete(idx);
            });
        });
    }

    function hideAutocomplete() {
        $autocomplete.classList.add('hidden');
        acItems = [];
        acSelected = -1;
    }

    function applyAutocomplete(idx) {
        if (!acTarget || idx < 0 || idx >= acItems.length) return;
        if (isMathFieldElement(acTarget)) return;
        const item = acItems[idx];
        const val = acTarget.value;
        const cursor = acTarget.selectionStart || 0;
        const before = val.substring(0, cursor);
        const match = before.match(/\\?([a-zA-Z_]\w*)$/);
        if (!match) return;

        const start = cursor - match[1].length;
        let insert = item.name;
        if (item.args > 0) insert += '(';
        const after = val.substring(cursor);
        acTarget.value = val.substring(0, start) + insert + after;
        acTarget.selectionStart = acTarget.selectionEnd = start + insert.length;

        // Trigger update
        acTarget.dispatchEvent(new Event('input'));
        hideAutocomplete();
    }

    function handleAutocompleteTab(e) {
        if (acItems.length === 0) return;
        e.preventDefault();
        if (acSelected < 0) acSelected = 0;
        applyAutocomplete(acSelected);
    }

    function handleAutocompleteNav(e) {
        if (acItems.length === 0) return;
        e.preventDefault();
        if (e.key === 'ArrowDown') {
            acSelected = (acSelected + 1) % acItems.length;
        } else {
            acSelected = acSelected <= 0 ? acItems.length - 1 : acSelected - 1;
        }
        $autocomplete.querySelectorAll('.ac-item').forEach((el, i) => {
            el.classList.toggle('selected', i === acSelected);
        });
    }

    // ========================================================
    // COMPUTATION PANEL
    // ========================================================
    configureMathField($computeInput);

    $computeInput.addEventListener('keydown', (e) => {
        if (handleLatexInputShortcuts(e, $computeInput)) return;
        if (e.key === 'Enter') {
            e.preventDefault();
            computeExpression(getFieldValue($computeInput));
        }
    });

    function computeExpression(input) {
        input = input.trim();
        if (!input) {
            $computeResult.classList.add('hidden');
            return;
        }

        $computeResult.classList.remove('hidden');
        $computeResult.classList.remove('error');

        try {
            const parsed = MathEngine.parseExpression(input);
            if (!parsed || parsed.type === 'error') {
                $computeResult.classList.add('error');
                $computeResult.innerHTML = `<span class="result-label">Error</span>${parsed ? parsed.error : 'Empty expression'}`;
                return;
            }

            let result;
            if (parsed.type === 'function' || parsed.type === 'slider') {
                const ast = parsed.ast || parsed;
                // Try to evaluate without variables
                try {
                    result = MathEngine.evaluate(ast, { ...sliders });
                } catch {
                    // If it needs x, try a few values
                    const vars = parsed.variables;
                    if (vars && vars.has('x')) {
                        $computeResult.innerHTML = `<span class="result-label">Function</span>
                            Needs variable x. Add to graph panel to plot.<br>
                            x=0 → ${tryEval(ast, 0)}<br>
                            x=1 → ${tryEval(ast, 1)}<br>
                            x=π → ${tryEval(ast, Math.PI)}`;
                        return;
                    }
                    throw new Error('Cannot evaluate');
                }
            } else {
                $computeResult.innerHTML = `<span class="result-label">Expression</span>Type: ${parsed.type}`;
                return;
            }

            if (typeof result === 'number' && isFinite(result)) {
                const isInt = Number.isInteger(result);
                let display = isInt ? result.toString() : result.toPrecision(12).replace(/0+$/, '').replace(/\.$/, '');
                $computeResult.innerHTML = `<span class="result-label">Result</span>= ${display}`;
            } else if (result instanceof MathEngine.Complex) {
                $computeResult.innerHTML = `<span class="result-label">Complex Result</span>= ${result.re.toFixed(6)} + ${result.im.toFixed(6)}i`;
            } else {
                $computeResult.innerHTML = `<span class="result-label">Result</span>= ${result}`;
            }
        } catch (e) {
            $computeResult.classList.add('error');
            $computeResult.innerHTML = `<span class="result-label">Error</span>${e.message}`;
        }
    }

    function tryEval(ast, x) {
        try {
            const v = MathEngine.evaluate(ast, { x, n: x, t: x, ...sliders });
            return formatNum(v);
        } catch { return '?'; }
    }

    // ========================================================
    // PERSISTENCE
    // ========================================================
    function queueSave() {
        clearTimeout(saveTimer);
        saveTimer = setTimeout(saveState, 140);
    }

    function saveState() {
        try {
            const theme = document.documentElement.getAttribute('data-theme') || 'dark';
            const graphMode = Grapher.getGraphMode();
            const view = Grapher.getView ? Grapher.getView() : null;
            const settings = {
                showGrid: document.getElementById('setting-grid').checked,
                showLabels: document.getElementById('setting-labels').checked,
                showMinorGrid: document.getElementById('setting-minor-grid').checked,
                lineWidth: parseFloat(document.getElementById('setting-line-width').value),
                pointSize: parseFloat(document.getElementById('setting-point-size').value),
                density: parseInt(document.getElementById('setting-density').value, 10),
            };

            const activeModeBtn = document.querySelector('.mode-btn.active');
            let plotMode = 'continuous';
            if (activeModeBtn?.id === 'btn-mode-discrete') plotMode = 'discrete';
            else if (activeModeBtn?.id === 'btn-mode-bar') plotMode = 'bar';

            const payload = {
                theme,
                graphMode,
                plotMode,
                settings,
                view,
                sidebarCollapsed: document.body.classList.contains('sidebar-collapsed'),
                presetsOpen: !$presetsList.classList.contains('hidden'),
                expressions: expressions.map(expr => ({
                    text: expr.text,
                    color: expr.color,
                    visible: expr.visible,
                    plotMode: expr.plotMode,
                })),
            };

            localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
        } catch {}
    }

    function loadState() {
        try {
            const raw = localStorage.getItem(STORAGE_KEY);
            if (!raw) return null;
            const parsed = JSON.parse(raw);
            return parsed && typeof parsed === 'object' ? parsed : null;
        } catch {
            return null;
        }
    }

    function setGlobalPlotMode(mode) {
        $modeBtns.forEach(b => b.classList.remove('active'));
        if (mode === 'discrete') document.getElementById('btn-mode-discrete').classList.add('active');
        else if (mode === 'bar') document.getElementById('btn-mode-bar').classList.add('active');
        else document.getElementById('btn-mode-continuous').classList.add('active');

        expressions.forEach(e => e.plotMode = mode);
        syncGrapher();
    }

    function applySavedState(saved) {
        if (!saved) return false;

        if (saved.theme === 'dark' || saved.theme === 'light') {
            document.documentElement.setAttribute('data-theme', saved.theme);
        }

        if (saved.settings) {
            if (typeof saved.settings.showGrid === 'boolean') {
                document.getElementById('setting-grid').checked = saved.settings.showGrid;
            }
            if (typeof saved.settings.showLabels === 'boolean') {
                document.getElementById('setting-labels').checked = saved.settings.showLabels;
            }
            if (typeof saved.settings.showMinorGrid === 'boolean') {
                document.getElementById('setting-minor-grid').checked = saved.settings.showMinorGrid;
            }
            if (typeof saved.settings.lineWidth === 'number') {
                document.getElementById('setting-line-width').value = String(saved.settings.lineWidth);
            }
            if (typeof saved.settings.pointSize === 'number') {
                document.getElementById('setting-point-size').value = String(saved.settings.pointSize);
            }
            if (typeof saved.settings.density === 'number') {
                document.getElementById('setting-density').value = String(saved.settings.density);
            }
            Grapher.updateSettings({
                showGrid: document.getElementById('setting-grid').checked,
                showLabels: document.getElementById('setting-labels').checked,
                showMinorGrid: document.getElementById('setting-minor-grid').checked,
                lineWidth: parseFloat(document.getElementById('setting-line-width').value),
                pointSize: parseFloat(document.getElementById('setting-point-size').value),
                density: parseInt(document.getElementById('setting-density').value, 10),
            });
        }

        clearExpressions(false);
        const savedExprs = Array.isArray(saved.expressions) ? saved.expressions : [];
        if (savedExprs.length > 0) {
            for (const item of savedExprs) {
                const expr = addExpression(typeof item.text === 'string' ? item.text : '', false);
                if (typeof item.color === 'string') expr.color = item.color;
                if (typeof item.visible === 'boolean') expr.visible = item.visible;
                if (typeof item.plotMode === 'string') expr.plotMode = item.plotMode;

                const $item = $exprList.querySelector(`[data-id="${expr.id}"]`);
                if ($item) {
                    const $bar = $item.querySelector('.expr-color-bar');
                    const $picker = $item.querySelector('.expr-color-picker');
                    $bar.style.background = expr.color;
                    $bar.classList.toggle('invisible', !expr.visible);
                    $picker.value = expr.color;
                }
            }
        } else {
            addExpression('sin(x)', false);
        }

        if (typeof saved.plotMode === 'string') {
            setGlobalPlotMode(saved.plotMode);
        }
        if (typeof saved.graphMode === 'string') {
            setGraphMode(saved.graphMode);
        }
        if (saved.view && Grapher.setView) {
            Grapher.setView(saved.view);
        }

        if (saved.sidebarCollapsed) {
            document.body.classList.add('sidebar-collapsed');
        }

        if (saved.presetsOpen) {
            $presetsList.classList.remove('hidden');
            $presetsToggle.classList.add('open');
        }

        syncGrapher();
        return true;
    }

    function toggleSidebar() {
        document.body.classList.toggle('sidebar-collapsed');
        queueSave();
    }

    function exportGraphPng() {
        if (!Grapher.exportPNG) return;
        const dataUrl = Grapher.exportPNG();
        const link = document.createElement('a');
        const stamp = new Date().toISOString().replace(/[:.]/g, '-');
        link.href = dataUrl;
        link.download = `mathgraph-${stamp}.png`;
        document.body.appendChild(link);
        link.click();
        link.remove();
    }

    function filterPresets(query) {
        const q = query.trim().toLowerCase();
        const children = Array.from($presetsList.children);
        let currentButtons = [];
        let currentLabel = null;

        const finalizeCategory = () => {
            if (!currentLabel) return;
            const anyVisible = currentButtons.some(btn => !btn.classList.contains('hidden-by-search'));
            currentLabel.classList.toggle('hidden-by-search', !anyVisible && q.length > 0);
        };

        for (const child of children) {
            if (child.classList.contains('preset-category-label')) {
                finalizeCategory();
                currentLabel = child;
                currentButtons = [];
                continue;
            }

            if (child.classList.contains('preset-btn')) {
                currentButtons.push(child);
                if (!q) {
                    child.classList.remove('hidden-by-search');
                } else {
                    const hay = `${child.textContent} ${child.dataset.expr || ''} ${child.dataset.expr2 || ''}`.toLowerCase();
                    child.classList.toggle('hidden-by-search', !hay.includes(q));
                }
            }
        }
        finalizeCategory();
    }

    // ========================================================
    // GRAPH MODE SWITCHING
    // ========================================================
    const $modeTabs = document.querySelectorAll('.graph-mode-tab');
    $modeTabs.forEach(tab => {
        tab.addEventListener('click', () => {
            $modeTabs.forEach(t => t.classList.remove('active'));
            tab.classList.add('active');
            const mode = tab.dataset.mode;
            Grapher.setGraphMode(mode);
            $modeIndicator.textContent = mode.charAt(0).toUpperCase() + mode.slice(1);
            queueSave();
        });
    });

    function setGraphMode(mode) {
        $modeTabs.forEach(t => {
            t.classList.toggle('active', t.dataset.mode === mode);
        });
        Grapher.setGraphMode(mode);
        $modeIndicator.textContent = mode.charAt(0).toUpperCase() + mode.slice(1);
        queueSave();
    }

    // ========================================================
    // PLOT MODE TOGGLE
    // ========================================================
    const $modeBtns = document.querySelectorAll('.mode-btn');
    $modeBtns.forEach(btn => {
        btn.addEventListener('click', () => {
            let mode = 'auto';
            if (btn.id === 'btn-mode-continuous') mode = 'continuous';
            else if (btn.id === 'btn-mode-discrete') mode = 'discrete';
            else if (btn.id === 'btn-mode-bar') mode = 'bar';

            setGlobalPlotMode(mode);
        });
    });

    // ========================================================
    // PRESETS
    // ========================================================
    $presetsToggle.addEventListener('click', () => {
        $presetsList.classList.toggle('hidden');
        $presetsToggle.classList.toggle('open');
        queueSave();
    });

    $presetSearch.addEventListener('input', () => {
        filterPresets($presetSearch.value);
        if ($presetSearch.value.trim() && $presetsList.classList.contains('hidden')) {
            $presetsList.classList.remove('hidden');
            $presetsToggle.classList.add('open');
        }
    });

    document.querySelectorAll('.preset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            const exprText = btn.dataset.expr;
            const mode = btn.dataset.mode;

            // Switch graph mode if specified
            if (mode) setGraphMode(mode);

            // Add expressions
            addExpression(exprText, false);

            // Second expression if present (for comparisons)
            const expr2 = btn.dataset.expr2;
            if (expr2) addExpression(expr2, false);

            queueSave();
        });
    });

    // ========================================================
    // ZOOM CONTROLS
    // ========================================================
    document.getElementById('btn-zoom-in').addEventListener('click', () => { Grapher.zoomIn(); queueSave(); });
    document.getElementById('btn-zoom-out').addEventListener('click', () => { Grapher.zoomOut(); queueSave(); });
    document.getElementById('btn-zoom-reset').addEventListener('click', () => { Grapher.resetView(); queueSave(); });

    // ========================================================
    // THEME TOGGLE
    // ========================================================
    document.getElementById('btn-theme').addEventListener('click', () => {
        const html = document.documentElement;
        const isDark = html.getAttribute('data-theme') === 'dark';
        html.setAttribute('data-theme', isDark ? 'light' : 'dark');
        Grapher.render();
        queueSave();
    });

    $btnSidebarToggle.addEventListener('click', () => {
        toggleSidebar();
    });

    $btnClearAll.addEventListener('click', () => {
        const shouldClear = confirm('Clear all expressions?');
        if (!shouldClear) return;
        clearExpressions(true);
        queueSave();
    });

    $btnExport.addEventListener('click', () => {
        exportGraphPng();
    });

    // ========================================================
    // REFERENCE & SETTINGS PANELS
    // ========================================================
    document.getElementById('btn-reference').addEventListener('click', () => {
        document.getElementById('reference-panel').classList.toggle('hidden');
        document.getElementById('settings-panel').classList.add('hidden');
    });
    document.getElementById('btn-close-ref').addEventListener('click', () => {
        document.getElementById('reference-panel').classList.add('hidden');
    });
    document.getElementById('btn-settings').addEventListener('click', () => {
        document.getElementById('settings-panel').classList.toggle('hidden');
        document.getElementById('reference-panel').classList.add('hidden');
    });
    document.getElementById('btn-close-settings').addEventListener('click', () => {
        document.getElementById('settings-panel').classList.add('hidden');
    });

    // ========================================================
    // SETTINGS BINDINGS
    // ========================================================
    document.getElementById('setting-grid').addEventListener('change', (e) => {
        Grapher.updateSettings({ showGrid: e.target.checked });
        queueSave();
    });
    document.getElementById('setting-labels').addEventListener('change', (e) => {
        Grapher.updateSettings({ showLabels: e.target.checked });
        queueSave();
    });
    document.getElementById('setting-minor-grid').addEventListener('change', (e) => {
        Grapher.updateSettings({ showMinorGrid: e.target.checked });
        queueSave();
    });
    document.getElementById('setting-line-width').addEventListener('input', (e) => {
        Grapher.updateSettings({ lineWidth: parseFloat(e.target.value) });
        queueSave();
    });
    document.getElementById('setting-point-size').addEventListener('input', (e) => {
        Grapher.updateSettings({ pointSize: parseFloat(e.target.value) });
        queueSave();
    });
    document.getElementById('setting-density').addEventListener('change', (e) => {
        Grapher.updateSettings({ density: parseInt(e.target.value) });
        queueSave();
    });

    // ========================================================
    // KEYBOARD SHORTCUTS
    // ========================================================
    document.addEventListener('keydown', (e) => {
        // Ctrl+Enter to add expression
        if (e.ctrlKey && e.key === 'Enter') {
            e.preventDefault();
            addExpression();
        }
        // Ctrl+/ to focus compute input
        if (e.ctrlKey && e.key === '/') {
            e.preventDefault();
            $computeInput.focus();
        }
        if (e.ctrlKey && e.key.toLowerCase() === 'b') {
            e.preventDefault();
            toggleSidebar();
        }
    });

    // ========================================================
    // ADD BUTTON
    // ========================================================
    $btnAdd.addEventListener('click', () => addExpression());

    // ========================================================
    // INIT — Restore previous state if available
    // ========================================================
    const restored = applySavedState(loadState());
    if (!restored) {
        addExpression('sin(x)', false);
    }

    // Close panels on outside click  
    document.addEventListener('click', (e) => {
        const refPanel = document.getElementById('reference-panel');
        const setPanel = document.getElementById('settings-panel');
        if (!refPanel.classList.contains('hidden') &&
            !refPanel.contains(e.target) &&
            !document.getElementById('btn-reference').contains(e.target)) {
            refPanel.classList.add('hidden');
        }
        if (!setPanel.classList.contains('hidden') &&
            !setPanel.contains(e.target) &&
            !document.getElementById('btn-settings').contains(e.target)) {
            setPanel.classList.add('hidden');
        }
    });

    window.addEventListener('beforeunload', saveState);

})();
