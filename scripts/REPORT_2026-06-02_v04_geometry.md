# Отчёт сессии — v0.4 «Consolidation» + геометрия многообразия

**Дата:** 2026-06-02 · **Инструмент:** perceptome v0.2.3 → v0.4 (в работе) · **Автор работы:** сессия с Claude
**Память проекта:** `project_tool_modernization_v04.md` (обновлена по ходу)

---

## 0. Отправная точка и задача

Сессия началась с обзора проекта **perceptome** (фреймворк: ~44 сигнальных модуля как «перцептивная система» клетки → 9-PC eigenspace из 154 нормальных типов HPA; папки papers 1-8; инструмент `tool/perceptome2`, v0.2.3).

Пользователь сформулировал цель: серии **paper 4.0-4.8** — это корпус примеров использования тула, установивших границы и идеи; нужно **свернуть их обратно В инструмент**, сделав «solid tool», который лучше описывает накопленное — **без динамики** (следующая эпоха) — плюс **геометрический анализ** для чёткости результатов.

**План модернизации (v0.4 «Consolidation»):**
- **Входит:** консолидация данных 4.0-4.8 → Factor 2 (предсказание величины ответа) → кодирование границ → геометрия как координатная инвариантность → hardening.
- **Откладывается:** динамика (Langevin, v1.0); тканеспецифичные модули (research-gated).
- **Жёсткие ограничения (locked):** PCA остаётся canonical (беспараметричность = щит от «curve-fit»); геометрические альтернативы только как robustness, unsupervised на 154 норме, рак — проекцией; никаких supervised-методов; diffusion как воспроизводимость, не замена; Factor 2 — порядковый (данных мало, strict 3a FAIL).

---

## 1. Phase 0 — консолидация substrate-серии ✅

**Скрипт:** `scripts/13_consolidate_substrate_series.py` → `substrate_series_v1.{json,csv}`

- **594 наблюдения** module-эффектов из 7 работ (paper4.3 immune, 4.4 muscle, 4.5 epithelial, 4.6/4.7/4.8 organoid, paper4 memory).
- Схема: `study · operation · context_class · condition · cell_type · module · score_type · effect · effect_metric · direction · verdict`.
- **Метрики не смешиваются:** cohens_d (paper4, 4.5) vs log2fc (4.3/4.4/4.6-4.8) — тегируются раздельно.
- Композиты `Cyto/ER_chaperones` помечены как не-catalog, исключены из калибровки.

---

## 2. Phase 1 — Factor 2 (предсказание величины) + границы ✅

**Скрипты:** `scripts/14_build_operation_intensity.py` → `operation_intensity_v1.json` (забандлен в `perceptome/perceptivity/data/`)

- **Factor-2 reference:** 228 записей (operation × module × metric), 7 операций. Порядковые полосы интенсивности (none/low/medium/high/extreme) + направление.
- **`predict_engagement(cell, modules, operation=...)`** расширен (заглушка `predicted_magnitude='unknown'` на floor.py:127 заменена). Обратно совместим: без `operation=` — старое поведение Factor-1.
- **Комбинация:** `capacity_floor × operation_intensity` → порядковая `predicted_magnitude` + направление.
- **Границы зашиты:** `none/flat` для низко-интенсивных операций (граница 4.5/4.6); upward-asymmetric (down при робастной супрессии |eff|≥0.30, напр. atRA→UPR-ATF6; up_blocked для насыщенных; down_weak для шума у нуля).
- **Честность в API:** `n_anchors`, `confidence` (low если n<2 или калибр-клетка насыщена), QC-флаг насыщения; disclaimer «calibration-grade, not validated».
- **Cell→HPA маппинг:** myofiber→myonuclei, organoid→enteric stem cells, ILC3→innate lymphoid cells, DG→brain excitatory neurons, Goblet/Enterocyte/Paneth — точно.
- **Тесты:** 79/79 (+6 регрессионных на substrate-биологию). Валидация воспроизводит вердикты: hypertrophy→extreme up, terminal-diff→boundary, retinoid→down, regeneration→Cell Cycle.

---

## 3. Phase 2 — геометрия как координатная инвариантность ✅

**Модули (аддитивно, PCA не тронут):** `eigenspace/rotate.py` (varimax), `eigenspace/manifold.py` (DiffusionMap: Coifman-Lafon + Nyström OOS, ε=median, не подбирается), `eigenspace/invariance.py` (coordinate_invariance + beacon_compactness).

**Результат:** в родной калибровке beacon'а (**single_cell_scaled**) конвергенция **координатно-инвариантна**:
- PCA percentile 0.000 (z=−3.48), varimax-rotated идентично, diffusion 0.0005 (z=−1.78) → **геометрия, не артефакт PCA**.
- Диагностика: «PCA adequate — diffusion no tighter» (подтверждает стратегию: координаты = подтверждение, не спасение).
- **Caveat (всплыл сам):** тесность beacon'а калибровка-специфична — в pseudobulk диффузна (0.46). Согласуется с обоснованием sc-scaled из CHANGELOG v0.2.3.
- Метрика валидирована на positive controls (линеажи retina/immune/germ тесные).
- **Тесты:** 84/84 (+5 геометрии). Sanity: rotation distance-identical ✓, Nyström round-trip 1e-16 ✓.

---

## 4. Phase 2b — глубокая геометрия многообразия

### 4.1 Форма (`scripts/16_manifold_geometry.py`)
- **Размерность:** внутренняя (TwoNN) ~4.8-6.5 vs линейная (participation ratio) ~8.1 → умеренно низкоразмерно, есть нелинейное сжатие.
- **Кривизна:** geodesic/euclidean ratio 1.6-1.85, но feature-shuffled null уже 1.6-1.7 → **слабая кривизна; PCA адекватен** (не Swiss-roll). Большая часть «искривления» = тяжёлые хвосты модулей + разреженность N=154.
- PC1 = 29% дисперсии.

### 4.2 Интерпретируемый базис (`scripts/18_rotated_axes.py`)
Varimax(PC1-6) → чистые оси: **R1** стресс/протеостаз (UPR/ERK/p53/NF-κB), **R2** развитие (Wnt/NPAS4/SREBP), **R3** сигналинг (cAMP/PI3K/NFAT), **R4** печёночная (PXR/PPARα/FXR — чисто!), **R5** репродуктивно-стероидная (PR/ER), **R6** механо-метаболическая (AMPK/Hippo). Varimax изолирует тканеспецифичные оси, размазанные в сырых PC.

### 4.3 Главная последовательность (`scripts/17_trajectory_mainsequence.py`)
Использован готовый корпус **`transition_vectors/`** (81 вектор, 16 классов, пре-регистрировано — НЕ циркулярно). 14 process-траекторий упорядочены на одной доминирующей оси ≈ R1 (cos −0.76) по полярности **ОСТРОЕ↔ХРОНИЧЕСКОЕ**:
- Острый полюс: ACTIVATION_ACUTE (−0.79), autoimmune, DKD.
- Хронический полюс: cancer (TUMOR_TME +0.85, TRANSFORMATION_EPITH +0.79), inflammaging (+0.92), neurodegen, fibrosis, late-consolidation.
- Подпись оси: Autophagy/UPR-PERK/VDR/NF-κB/HSF1/GR/HIF.
- Доминирующая ось = 38% дисперсии набора, топ-2 = 57%, эфф. dim ~6.75 → **«главная последовательность с разбросом»**, не идеальная 1D-кривая.

**⚠️ Коррекция:** организующая ось — острое↔хроническое, **НЕ дифференцировка**. Чистый тест: cos(рак, differentiation-completion) = −0.14…+0.28 (≈ортогонально). Раннее «−0.62 рак⟂дифференцировке» было артефактом циркулярного `attractor_direction` (его клетки пересекаются с эмбриональным полюсом) + z-score кадра.

### 4.4 HR-вердикт
Аналогия Hertzsprung-Russell **держится в точной форме**: одна доминирующая интерпретируемая ось (R1), вдоль которой гетерогенные процессы упорядочены по острое↔хроническое — видимая в **линейных** координатах (varimax). Универсальность видна в простейших координатах. Организующий принцип — **хроническая реорганизация**, а не дифференцировка.

---

## 5. Раковые траектории — где рак на главной последовательности

### 5.1 Не-циркулярный HCC (`scripts/19_cancer_vector_noncircular.py`)
Пересчёт Sun2021 HCC из сырья (paired, без attractor-прокси): **малигнантный сдвиг +0.47** к хроническому полюсу; TME тоже дрейфует (CAF +0.35, Myeloid +0.50, T/NK +0.61), endothelial/B ≈0.

### 5.2 Картинка 19 раков (`scripts/20_cancer_mainsequence_figure.py` → `fig_cancer_mainsequence.png`)
Попытка спроецировать абсолютные позиции выявила и устранила **две ловушки** (методологическая дисциплина):
1. **Баг имён генов:** Census-кэши хранят символы в `var['feature_name']`, var_names = '0','1',… → 0/15 матчей → коллапс. Исправлено ремапом.
2. **Конфляция:** абсолютная ms-позиция НЕ отражает конвергенцию — Spearman(позиция, conv_8beacon)=−0.30, p=0.28 (ns); нейробластома на крайнем полюсе из-за *нейральной* программы, не дуэта.

→ Итоговая картинка из **робастных величин**: (A) 14 process-классов на главной оси (рак на хроническом полюсе); (B) **все 19 раков** конвергируют на дуэт/beacon (валидированный conv_8beacon) — **19/19 выше базы 5.2%**, медиана 38.9%, макс 92% (MM).

### 5.3 Чистые tumor−normal сдвиги (`scripts/21` + `scripts/22` → `fig_shift_directions.png`)
Калибровка-независимые сдвиги (offset сокращается в разности):
- **TIER 1 (within-dataset, gold):** HCC +0.47, BCC +0.85, BRCA **−0.62** (все 3 подтипа отрицательны).
- **TIER 2 (cross-Census cell-of-origin, приближённо):** RCC +0.55, PDAC +0.90, Glioblastoma +0.63, DLBCL −0.64, Neuroblastoma +0.02.
- **6/8 хронический, 2/8 острый (BRCA, DLBCL) → направление РАСЩЕПЛЕНО, тканеспецифично.**

**⚠️ Главный вывод по ракам:** направление tumor−normal на главной оси **тканеспецифично, НЕ универсально**. Робастное универсальное утверждение — **позиционная конвергенция** (conv_8beacon, 19/19): дуэт — это аттрактор-МИШЕНЬ, к которой ткани приходят с разных сторон (знак сдвига = функция происхождения). Это ровно заключение v7 («direction-тест непостоянен → используем convergence-fraction»), подтверждённое здесь независимо.

### 5.4 Граница охвата (честно)
Чистые сдвиги достижимы только для **3 раков** (within-dataset normal: HCC/BRCA/BCC) + ~5 приближённых (cross-Census). Остальные ~11 — только `*_malignant.h5ad` без нормы; census-атласы иммунно-смещены (печень/кишка/кожа без паренхимы); атласов breast/ovary/prostate/head-neck нет. **Чистые 19 локально невозможны** — нужен полноценный per-cell нормальный атлас с паренхимой (HPA-single-cell / Tabula Sapiens по органам).

---

## 6. Сводка ключевых находок

1. **Factor 2 операционализирован** — двухфакторный capacity-предиктор (порядковый, честный про лимиты данных).
2. **Beacon-конвергенция координатно-инвариантна** (PCA + rotated + diffusion согласны) → геометрия, не артефакт.
3. **Многообразие** умеренно-низкоразмерно (~5), **слабо искривлено** → PCA адекватен.
4. **Клеточная главная последовательность** = ось острое↔хроническое (R1); 14 процессов упорядочены на ней — HR-аналог в линейных координатах.
5. **Направление рака тканеспецифично** (6 хрон./2 остр. из 8); **позиционная конвергенция** — робастный универсальный закон.
6. **Две коррекции:** рак ⟂ дифференцировке (не −0.62, артефакт); beacon → дуэт megakaryocyte+cytotrophoblast (v7, не старый 8-cell).

---

## 7. Артефакты

**Код инструмента (изменён/добавлен):** `perceptivity/floor.py` (+Factor 2), `perceptivity/reference.py` (+load_operation_intensity), `eigenspace/{rotate,manifold,invariance}.py` (новые), оба `__init__.py` (экспорты), `tests/test_perceptivity.py` + `tests/test_geometry.py` (+11 тестов).
**Данные:** `operation_intensity_v1.json` (забандлен). `substrate_series_v1.{json,csv}`, `coordinate_invariance_v1.json` (staging в scripts/).
**Скрипты анализа:** `scripts/13`–`scripts/22`.
**Картинки:** `fig_cancer_mainsequence.png`, `fig_shift_directions.png`.
**Тесты:** 84/84 проходят.
**Память:** `project_tool_modernization_v04.md` — полная хронология Phases 0/1/2/2b + коррекции + caveats.

**Dev-нюанс:** чистый `import perceptome` из стороннего cwd резолвится в СТАРЫЙ v0.1 (`tool/perceptome`, pip-установлен) — perceptome2 берётся при cwd=perceptome2 или sys.path.insert (скрипты 13-22 это учитывают).

---

## 8. Что дальше (на выбор)

1. **Оформить геометрию в инструмент** — `pct.manifold_report()` (dim/curvature/null-controls) + varimax-базис с метками R1-R6 + main-sequence координата; перенести картинки в `manuscript_figures/`.
2. **Phase 3 (release)** — SCOPE.md (Factor 2 значится «не реализован»), CHANGELOG/ROADMAP (переименовать волну в v0.4), recipe-notebook (magnitude + coordinate-invariance), edge-cases/Tabula Sapiens; тег + Zenodo.
3. **Нормальный атлас** (Tabula Sapiens по органам) → чистые сдвиги для всех 19.
