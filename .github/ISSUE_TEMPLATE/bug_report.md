---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Standalone code to reproduce the error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Actual behavior**
Please include the full traceback of any errors

**System information:**

Output of `magic.__version__`:

```
If you are running MAGIC in R or Python, please run magic.__version__ and paste the results here.

You can do this with `python -c 'import magic; print(magic.__version__)'`
```

Output of `pd.show_versions()`:

<details>

```
If you are running MAGIC in R or Python, please run pd.show_versions() and paste the results here.

You can do this with `python -c 'import pandas as pd; pd.show_versions()'`
```

</details>

Output of `sessionInfo()`:

<details>

```
If you are running MAGIC in R, please run sessionInfo() and paste the results here.

You can do this with `R -e 'library(Rmagic); sessionInfo()'`
```

</details>

Output of `reticulate::py_discover_config(required_module = "magic")`:

<details>

```
If you are running MAGIC in R, please run `reticulate::py_discover_config(required_module = "magic")` and paste the results here.

You can do this with `R -e 'reticulate::py_discover_config(required_module = "magic")'`
```

</details>

**Additional context**
Add any other context about the problem here.
