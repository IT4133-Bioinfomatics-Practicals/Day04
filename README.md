# 📊 Day 04: Gene Expression Data Visualization Guide

A complete guide to visualizing cancer gene expression data with 7 different plot types.

---

## ⚙️ Quick Setup

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

# Load data
data = pd.read_csv("cancer_expression.csv")

# Identify columns
control_cols = [col for col in data.columns if 'Control' in col]
treatment_cols = [col for col in data.columns if 'Treatment' in col]
```

---

## 🌋 1. Volcano Plot - Find Significant Genes

**Purpose:** Identify genes with significant expression changes between control and treatment.

**Key Concept:** Shows fold change (x-axis) vs p-value (y-axis)

```python
# Calculate means
data['Mean_Control'] = data[control_cols].mean(axis=1)
data['Mean_Treatment'] = data[treatment_cols].mean(axis=1)

# Calculate log2 fold change
data['log2FoldChange'] = np.log2(data['Mean_Treatment'] + 1e-6) - np.log2(data['Mean_Control'] + 1e-6)

# Perform t-tests
ttest_results = ttest_ind(data[treatment_cols].values, data[control_cols].values, axis=1, equal_var=False)
data['p-value'] = ttest_results.pvalue
data['-log10(p-value)'] = -np.log10(data['p-value'])

# Classify significance
fc_thresh = 1.5
p_thresh = 0.05
data['Significance'] = 'Not Significant'
data.loc[(data['log2FoldChange'] > fc_thresh) & (data['p-value'] < p_thresh), 'Significance'] = 'Significant'

# Plot
plt.figure(figsize=(10,6))
sns.scatterplot(data=data, x='log2FoldChange', y='-log10(p-value)', hue='Significance',
                palette={'Not Significant':'grey','Significant':'red'})
plt.axvline(fc_thresh, color='blue', linestyle='--')
plt.axvline(-fc_thresh, color='blue', linestyle='--')
plt.axhline(-np.log10(p_thresh), color='green', linestyle='--')
plt.title('Volcano Plot')
plt.xlabel('Log2 (Fold Change)')
plt.ylabel('-Log10(p-value)')
plt.show()
```

**What to Remember:** 🔴 Red dots = significant genes | 🔵 Blue dashed lines = thresholds

---

## 🔥 2. Heatmap - Visual Expression Patterns

**Purpose:** Display expression levels across samples as a color grid.

**Key Concept:** Red = high expression, Blue = low expression

```python
# Prepare data
heatmap_data = data.set_index('Gene')
heatmap_data = heatmap_data[control_cols + treatment_cols].head(10)

# Plot
plt.figure(figsize=(12,8))
sns.heatmap(heatmap_data, cmap="coolwarm", annot=True)
plt.title('Heatmap of Gene Expression (Top 10 Genes)')
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.tight_layout()
plt.show()
```

**What to Remember:** Each cell = expression level | 🎨 Color intensity = magnitude

---

## 🕷️ 3. Radar Plot (Spider Net) - Multi-Dimensional View

**Purpose:** Compare expression across multiple samples for a single gene.

**Key Concept:** Each axis = one sample, distance from center = expression level

```python
from math import pi

# Choose gene
gene = "Gene10"
gene_row = data[data["Gene"] == gene].iloc[0]
values = gene_row[control_cols + treatment_cols].values

# Setup angles
labels = control_cols + treatment_cols
num_vars = len(labels)
angles = [n / float(num_vars) * 2 * pi for n in range(num_vars)]
angles += angles[:1]
values = np.concatenate((values, [values[0]]))

# Plot
plt.figure(figsize=(8,8))
ax = plt.subplot(111, polar=True)
plt.xticks(angles[:-1], labels, color='grey', size=8)
ax.plot(angles, values, color='red', linewidth=2, linestyle='solid')
ax.fill(angles, values, color='skyblue', alpha=0.4)
plt.title(f'Radar Plot: {gene}')
plt.tight_layout()
plt.show()
```

**What to Remember:** ⭕ Circular plot | ⬆️ Larger area = higher overall expression

---

## ⬆️ 4. UpSet Plot - Intersection Analysis

**Purpose:** Show which genes have high expression in specific sample combinations.

**Key Concept:** Identifies overlaps between conditions

```python
from upsetplot import UpSet, from_memberships

# Define threshold
threshold = 10
memberships = []

# Categorize genes
for _, row in data.iterrows():
    sets = []
    if row['Control1'] > threshold:
        sets.append('Control1')
    if row['Control2'] > threshold:
        sets.append('Control2')
    if row['Treatment1'] > threshold:
        sets.append('Treatment1')
    memberships.append(sets)

# Create data
upset_data = from_memberships(memberships)

# Plot
plt.figure(figsize=(10, 6))
UpSet(upset_data, subset_size='count').plot()
plt.title("UpSet Plot: Gene Expression Intersections")
plt.show()
```

**What to Remember:** Shows 🔄 Venn diagram intersections in bar format

---

## 🎻 5. Violin Plot - Expression Distribution

**Purpose:** Show distribution shape of expression values.

**Key Concept:** Width = frequency, compares Control vs Treatment distributions

```python
gene = "Gene2"
gene_row = data[data["Gene"] == gene].iloc[0]

# Get control and treatment values
control_vals = gene_row[control_cols].values
treatment_vals = gene_row[treatment_cols].values

# Plot
data_violin = [control_vals, treatment_vals]
sns.violinplot(data=data_violin)
plt.xticks([0, 1], ['Control', 'Treatment'])
plt.ylabel('Expression Level')
plt.title(f'Violin Plot: {gene} Distribution')
plt.show()
```

**What to Remember:** 📋 Wider = more genes at that expression level

---

## 📊 6. Bar Chart - Compare Averages

**Purpose:** Simple comparison of mean expression between conditions.

**Key Concept:** Height = average expression value

```python
gene = "Gene10"
gene_data = data[data["Gene"] == gene]

# Calculate means
control_mean = gene_data[control_cols].mean(axis=1).values[0]
treatment_mean = gene_data[treatment_cols].mean(axis=1).values[0]

# Plot
plt.figure(figsize=(8, 6))
plt.bar(["Control", "Treatment"], [control_mean, treatment_mean],
        color=["blue", "red"])
plt.ylabel('Average Expression Level')
plt.title(f'Bar Chart: {gene} Average Expression')
plt.show()
```

**What to Remember:** ⚡ Quick visual comparison | 👁️ Easy to interpret

---

## 💫 7. Scatter Plot - Individual Sample Values

**Purpose:** Show individual data points across samples.

**Key Concept:** Each point = one sample, x-axis = sample index

```python
gene = "Gene10"

# Get values
control_values = data.loc[data['Gene']==gene, control_cols].values.flatten()
treatment_values = data.loc[data['Gene']==gene, treatment_cols].values.flatten()

# Plot
plt.figure(figsize=(10,6))
plt.scatter(range(len(control_values)), control_values,
            color='blue', label='Control', s=100, alpha=0.7)
plt.scatter(range(len(treatment_values)), treatment_values,
            color='red', label='Treatment', s=100, alpha=0.7)
plt.xlabel('Sample Index')
plt.ylabel('Expression Level')
plt.title(f'Scatter Plot: {gene}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

**What to Remember:** 🎯 All individual data points visible | 👀 See outliers easily

---

## 🧠 Memory Tips

| Plot Type      | Use When                             | Quick Code        |
| -------------- | ------------------------------------ | ----------------- |
| 🌋 **Volcano** | Finding significant genes            | log2FC vs p-value |
| 🔥 **Heatmap** | Seeing patterns across genes/samples | Color grid        |
| 🕷️ **Radar**   | Comparing one gene across samples    | Circular plot     |
| ⬆️ **UpSet**   | Finding overlaps between conditions  | Intersection bars |
| 🎻 **Violin**  | Understanding distribution shapes    | Density curves    |
| 📊 **Bar**     | Simple comparison of means           | Two bars          |
| 💫 **Scatter** | Showing individual points            | Dots on grid      |

---

## Common Errors & Fixes

**❌ Error:** `'list' object has no attribute 'set_index'`

- ✅ **Fix:** Don't overwrite `data` with a list. Use different variable names like `violin_data`

**❌ Error:** `NameError: name 'angles' is not defined`

- ✅ **Fix:** Calculate angles: `angles = [n / float(num_vars) * 2 * pi for n in range(num_vars)]`

**❌ Error:** FutureWarnings from upsetplot

- ✅ **Fix:** Add `warnings.filterwarnings('ignore', category=FutureWarning)` at top

---

## 💡 Essential One-Liners

```python
# Get column groups
control_cols = [col for col in data.columns if 'Control' in col]
treatment_cols = [col for col in data.columns if 'Treatment' in col]

# Calculate mean
gene_data[control_cols].mean(axis=1).values[0]

# T-test
ttest_results = ttest_ind(data[treatment_cols].values, data[control_cols].values, axis=1, equal_var=False)

# Log2 fold change
np.log2(treatment/control + 1e-6)

# -log10 p-value
-np.log10(p_value)
```

---

## 🚀 Workflow Summary

1. 📦 **Load & Prepare** → Load CSV, identify Control/Treatment columns
2. 🧮 **Calculate Stats** → Fold changes, p-values, t-tests
3. 🏷️ **Classify** → Mark significant genes
4. 📈 **Visualize** → Choose appropriate plot for question
5. 🔍 **Interpret** → Red/high values = significant changes

---

**Created for IT 4133 - Bioinformatics and Computational Biology**
