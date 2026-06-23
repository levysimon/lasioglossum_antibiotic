# activity-and-behavior

R Markdown analysis pipeline for social behavior and locomotor activity under antibiotic treatment and gut reinoculation.

---

## Experimental design

42 arenas across two years (2024–2025), in a **3 × 2 factorial design**:

| Factor | Levels |
|---|---|
| Antibiotic dose | zero, low, high |
| Gut reinoculation | no, yes |

A **reduced 2 × 2 model** is also fitted by collapsing `low` and `zero` antibiotic doses into a single `zero` category.

---

## Data

| File | Description |
|---|---|
| `data/behavior_csv_raw/YYYY_petridishN.csv` | Raw per-arena behavioral event logs (30 min, classified by event type) |
| `data/animal_ta/inderindividual_dist.csv` | Mean inter-individual distance per arena [cm] |
| `data/animal_ta/moving.csv` | Proportion of time moving and average speed per arena |

---

## Analysis pipeline (`activity_and_behavior.Rmd`)
### 1. Data pre-processing
### 2. Behavioral analysis (4 metrics)
### 3. Inter-individual distance
### 4. Activity (proportion time moving, average speed)
### 5. Multivariate analysis (PCA on 7 metrics)
### 6. Combined figure
