#!/usr/bin/env bash
# Build all KS4 atlas pages, update gallery.html, commit the result to the
# branch "atlas-update", push it and print the link to merge it on GitHub.
# Your current branch and working files are left as they are, so the images
# stay on disk as a cache for the next run.
#
#   bash tools/update_atlas.sh <HiPS base URL or local directory>
#
# or set KS4_HIPS_URL once (e.g. in ~/.bashrc) and run without an argument.
# Extra arguments after the URL go to make_atlas.py (e.g. --dec-limit -30).
set -euo pipefail

HIPS="${1:-${KS4_HIPS_URL:-}}"
[ $# -gt 0 ] && shift
if [ -z "$HIPS" ]; then
  echo "Usage: bash tools/update_atlas.sh <HiPS base URL or directory>  (or set KS4_HIPS_URL)" >&2
  exit 1
fi

cd "$(git rev-parse --show-toplevel)"
BASE=$(git rev-parse --abbrev-ref HEAD)
BRANCH=atlas-update
GEN=()
for c in peculiar ngc messier-caldwell planetary nebulae clusters; do
  GEN+=("$c.html" "assets/$c")
done

if ! git diff --quiet -- gallery.html; then
  echo "gallery.html has uncommitted changes. Commit or discard them first." >&2
  exit 1
fi

echo "== Updating $BASE from GitHub"
git fetch origin "$BASE"
# Pages and images are generated, and GitHub has the versions from the last merged run.
# Let those win so the update never stops on them: local edits to tracked generated files
# are reset, and untracked generated files that GitHub also has are removed (both are
# written again below; images that are unchanged come back from GitHub, not from KS4).
git diff --name-only -- "${GEN[@]}" | while IFS= read -r f; do git checkout -- "$f"; done
git ls-files --others --exclude-standard -- "${GEN[@]}" | while IFS= read -r f; do
  if git cat-file -e "origin/$BASE:$f" 2>/dev/null; then rm -f "$f"; fi
done
git merge --ff-only "origin/$BASE"

echo "== Building atlas pages"
python tools/make_atlas.py all --hips-url "$HIPS" "$@"

echo "== Committing to $BRANCH"
PATHS=(gallery.html)
for p in "${GEN[@]}"; do [ -e "$p" ] && PATHS+=("$p"); done
TMP_INDEX=$(mktemp -u)
trap 'rm -f "$TMP_INDEX"' EXIT
export GIT_INDEX_FILE="$TMP_INDEX"
git read-tree HEAD
git add -- "${PATHS[@]}"
TREE=$(git write-tree)
unset GIT_INDEX_FILE
git checkout -- gallery.html   # the new gallery.html lives on the branch; keep $BASE clean
if [ "$TREE" = "$(git rev-parse 'HEAD^{tree}')" ]; then
  echo "Nothing changed."
  exit 0
fi
COMMIT=$(git commit-tree "$TREE" -p HEAD -m "Update KS4 atlas pages")
git push -f origin "$COMMIT:refs/heads/$BRANCH"

echo
echo "Merge on GitHub:"
echo "  https://github.com/seowony/seowony.github.io/compare/$BASE...$BRANCH?expand=1"
