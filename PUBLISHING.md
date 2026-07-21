# Publishing pywbgt to PyPI

`pywbgt` publishes to PyPI automatically whenever a version tag (`v*`) is pushed.
Publishing uses **PyPI Trusted Publishing (OIDC)** — there are **no API tokens
or secrets** to manage.

## Files involved

| File | Purpose |
|------|---------|
| `.github/workflows/release.yml` | Verifies the tag, builds wheels + sdist, and publishes — **only on tag push**. |
| `.github/workflows/ci.yml` | Builds and runs the test suite on PRs / pushes to `main`. Never publishes. |
| `setup.py` | Compiles the C/Cython extensions with cross-platform OpenMP flags. |
| `pyproject.toml` | Project metadata + `[tool.cibuildwheel]` build config. |
| `MANIFEST.in` | Ensures the sdist ships all C/Cython sources so a source install compiles. |

## One-time setup

### 1. Register the Trusted Publisher on PyPI

Because `pywbgt` isn't on PyPI yet, use a **pending publisher**:

1. Log in to <https://pypi.org> (verify your email first).
2. Go to **Account settings → Publishing → Add a pending publisher**.
3. Fill in exactly:
   - **PyPI Project Name:** `pywbgt`
   - **Owner:** `kwodzicki`
   - **Repository name:** `pywbgt`
   - **Workflow name:** `release.yml`
   - **Environment name:** `pypi`
4. Save.

### 2. Create the GitHub environment

In the repo: **Settings → Environments → New environment**, name it `pypi` (must
match step 1 and the `environment: pypi` in `release.yml`). Optionally add
protection rules (e.g. a required reviewer) so a release can't fire by accident.

> Note: publishing pushes to the GitHub repo at `github.com/kwodzicki/pywbgt`
> (the `personal` git remote). Make sure your release tags are pushed there.

## Install modes (committed C vs. re-cythonize)

The package ships both the Cython sources (`.pyx`) and the pre-generated C
(`.c`). Which gets used depends on how it's installed:

| Command | What happens |
|---------|--------------|
| `pip install pywbgt` | Installs a prebuilt wheel if one matches the platform — **no compilation at all**. If no wheel matches, builds from the sdist using the committed `.c` (needs a C compiler + OpenMP, **not** Cython). |
| `pip install --no-binary pywbgt` | Forces a source build even when a wheel exists, compiling the committed `.c`. Reproducible; still no Cython. |
| `PYWBGT_CYTHONIZE=1 pip install --no-binary pywbgt` | Forces a source build that **regenerates C from the `.pyx`** sources. Requires Cython (pip installs it automatically in the isolated build env). |

The switch is the `PYWBGT_CYTHONIZE` environment variable, read by `setup.py`:
unset/`0` → compile committed `.c` (default); `1` → re-cythonize from `.pyx`.

CI sets `PYWBGT_CYTHONIZE=1` (in `[tool.cibuildwheel].environment` and in
`ci.yml`) so published wheels and PR test builds always track the current
`.pyx`, making `.pyx` the source of truth and the committed `.c` a convenience
snapshot for the minimal-dependency source path.

## Cutting a release

1. Bump `version` in `pyproject.toml` (e.g. `3.0.5` → `3.0.6`). It must equal the
   tag **without** the leading `v` — `release.yml` fails fast otherwise.
2. Commit, tag, and push:
   ```bash
   git commit -am "Release 3.0.6"
   git tag v3.0.6
   git push personal main --tags     # remote that maps to github.com/kwodzicki/pywbgt
   ```
3. Watch the run under the repo's **Actions** tab. On success the new version
   appears on <https://pypi.org/project/pywbgt/> and installs with:
   ```bash
   pip install pywbgt
   ```

## Notes & gotchas

- **Only tags trigger publishing.** `release.yml` listens solely to
  `push.tags: ["v*"]`; ordinary commits and PRs never publish.
- **Version is single-sourced in `pyproject.toml`** and edited by hand; the
  tag-vs-version check guards against mismatches.
- **A version can only be uploaded once.** If a publish fails after some files
  uploaded, bump the version and re-tag — PyPI rejects re-uploads.
- **Regenerate C sources after editing any `.pyx`:** run
  `cythonize src/pywbgt/*.pyx` (or `PYWBGT_CYTHONIZE=1 python -m build`) and
  commit the updated `.c` files, so the default (committed-C) install path
  doesn't ship stale code. CI wheels re-cythonize regardless.
- **macOS is the riskiest build.** CI uses Apple clang + Homebrew `libomp`
  (your local setup used Homebrew GCC). `ci.yml` exercises this path on every
  PR so problems surface before you tag.
- **Test on TestPyPI first (optional).** Add
  `with: { repository-url: https://test.pypi.org/legacy/ }` to the publish step
  and register a matching pending publisher on <https://test.pypi.org>.
