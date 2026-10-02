#!/usr/bin/env bash
# Fetch what tools/plsim.py needs to run PrairieLearn questions offline, into tools/.cache/
# (git-ignored). Safe to re-run; it only downloads what is missing.
#
#   bash tools/setup_harness.sh            # pinned PrairieLearn commit (known to work)
#   PL_REF=master bash tools/setup_harness.sh   # newest PrairieLearn instead
#
# Needs git and python3 (3.12+; PrairieLearn itself runs 3.13). Downloads only from GitHub:
#   - PrairieLearn (sparse): apps/prairielearn/{elements,python,src/schemas} + docs/{elements,question,assessment}
#   - pure-Python packages, only if `import` fails: sympy 1.14.0, mpmath 1.3.0, chevron,
#     text-unidecode, coloraide 8.12.1; plus a tiny stand-in for pint (units are not used)
#   - Bootstrap 5.3.8 CSS (for tools/preview.py)
# numpy, lxml, pandas, networkx and jsonschema must already be importable (pip install them if not;
# scipy is needed by the answer checkers).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE="${PL_HARNESS_CACHE:-$HERE/.cache}"
PL_REF="${PL_REF:-08f657e309e305b2d19ea12bdeaf25c367ba94be}"   # PrairieLearn master, 2026-10-01
mkdir -p "$CACHE/pylib"
cd "$CACHE"

# ---------------------------------------------------------------- PrairieLearn (sparse, shallow)
if [ ! -d PrairieLearn/apps/prairielearn/elements ]; then
  echo ">> fetching PrairieLearn ($PL_REF)"
  rm -rf PrairieLearn
  git init -q PrairieLearn
  git -C PrairieLearn remote add origin https://github.com/PrairieLearn/PrairieLearn.git
  git -C PrairieLearn sparse-checkout set --no-cone \
    /apps/prairielearn/elements/ /apps/prairielearn/python/ /apps/prairielearn/src/schemas/ \
    /docs/elements/ /docs/question/ /docs/assessment/ /docs/course/
  git -C PrairieLearn fetch -q --depth 1 --filter=blob:none origin "$PL_REF"
  git -C PrairieLearn -c advice.detachedHead=false checkout -q FETCH_HEAD
else
  echo ">> PrairieLearn already present ($(git -C PrairieLearn rev-parse --short HEAD))"
fi

# ---------------------------------------------------------------- pure-Python dependencies
need() {  # need <import-name>  -> true if python cannot import it from the system or the cache
  ! PYTHONPATH="$CACHE/pylib" python3 -c "import $1" 2>/dev/null
}
clone_pkg() {  # clone_pkg <repo-url> <tag-or-empty> <dir-in-repo> <import-name>
  local url="$1" tag="$2" sub="$3" name="$4" dst="src-$4"
  if need "$name"; then
    echo ">> fetching $name"
    rm -rf "$dst"
    if [ -n "$tag" ]; then git -c advice.detachedHead=false clone -q --depth 1 --branch "$tag" "$url" "$dst"; else git clone -q --depth 1 "$url" "$dst"; fi
    ln -sfn "../$dst/$sub" "pylib/$name"
  fi
}
clone_pkg https://github.com/mpmath/mpmath.git 1.3.0 mpmath mpmath
clone_pkg https://github.com/sympy/sympy.git sympy-1.14.0 sympy sympy
clone_pkg https://github.com/noahmorrison/chevron.git "" chevron chevron
clone_pkg https://github.com/kmike/text-unidecode.git "" src/text_unidecode text_unidecode
clone_pkg https://github.com/facelessuser/coloraide.git 8.12.1 coloraide coloraide
if need pint; then
  mkdir -p pylib/pint
  cat > pylib/pint/__init__.py <<'EOF'
"""Stand-in for pint so prairielearn.misc_utils imports; units (pl-units-input) are unsupported."""


class UnitRegistry:  # pragma: no cover
    def __init__(self, *args, **kwargs):
        raise RuntimeError("pint is stubbed in the offline harness")
EOF
fi

# ---------------------------------------------------------------- Bootstrap CSS for previews
if [ ! -f bootstrap.min.css ]; then
  echo ">> fetching Bootstrap CSS"
  git -c advice.detachedHead=false clone -q --depth 1 --branch v5.3.8 --filter=blob:none --sparse https://github.com/twbs/bootstrap.git src-bootstrap
  git -C src-bootstrap sparse-checkout set dist/css >/dev/null
  cp src-bootstrap/dist/css/bootstrap.min.css bootstrap.min.css
  rm -rf src-bootstrap
fi

# ---------------------------------------------------------------- sanity check
for m in numpy lxml pandas networkx jsonschema scipy; do
  python3 -c "import $m" 2>/dev/null || echo "!! python module '$m' is missing: pip install $m"
done
PYTHONPATH="$CACHE/pylib:$CACHE/PrairieLearn/apps/prairielearn/python" python3 -c \
  "import sympy, chevron, prairielearn; print('>> harness ready: sympy', sympy.__version__)"
