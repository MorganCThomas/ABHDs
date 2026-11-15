#!/usr/bin/env bash
# Build/update the local `public` branch from `main` for a free public GitHub repository.
#
# What it does: takes the tree of SRC (default main), drops raw data, large generated files, archives and unpublished working material,
# applies text redactions, checks the result, and adds ONE commit on top of the existing `public` branch (a fast-forward, so a plain push works).
# It uses git plumbing and a temporary index: your working tree and current branch are never touched.
#
# Usage:  scripts/publish_public.sh [-n] [-m "message"] [SRC]
#   -n   dry run: report what would be published, make no commit
# Then:   git push origin public:main        (the script never pushes)
set -euo pipefail

DRY=0; MSG="Update public snapshot"
while getopts "nm:" o; do case $o in n) DRY=1;; m) MSG=$OPTARG;; *) exit 2;; esac; done
shift $((OPTIND - 1)); SRC=${1:-main}; BR=public

# ---- rules (edit here) -------------------------------------------------------------------------------------------------------
EXCLUDE_RE='^(archive/|pymol/|data/raw/|docs/ABHDs\.pptx$|figures/[^/]+\.html$)|\.(npz|pt|pse)$'
KEEP_RE='^data/raw/uniprot/accessions\.txt$'            # tiny frozen input list that `make data` starts from
# Redactions live in an untracked local file so that this script never contains the strings it removes.
# Format, one per line, tab-separated:  FILE <TAB> TEXT_TO_REPLACE <TAB> REPLACEMENT      (file: .publish_redactions, git-ignored)
REDACTIONS=".publish_redactions"
FORBID_RE='christina/|\.pptx$|\.pse$|\.npz$'            # must not appear in the result
MAX_BYTES=$((20 * 1024 * 1024))
# -------------------------------------------------------------------------------------------------------------------------------

cd "$(git rev-parse --show-toplevel)"
git rev-parse -q --verify "$SRC^{commit}" >/dev/null || { echo "unknown source: $SRC" >&2; exit 1; }
export GIT_INDEX_FILE; GIT_INDEX_FILE=$(mktemp -u); trap 'rm -f "$GIT_INDEX_FILE"' EXIT

git ls-tree -r "$SRC" | awk -F'\t' -v ex="$EXCLUDE_RE" -v keep="$KEEP_RE" '{ if ($2 ~ keep || $2 !~ ex) print }' \
  | awk -F'\t' '{ split($1, a, " "); print a[1] " " a[3] "\t" $2 }' | git update-index --add --index-info

# redactions (index entries replaced by new blobs; the working tree is untouched)
SENSITIVE=()
if [ -f "$REDACTIONS" ]; then
  while IFS=$'\t' read -r file from to; do
    [ -z "${file:-}" ] || [ "${file:0:1}" = "#" ] && continue
    SENSITIVE+=("$from")
    git cat-file -e "$SRC:$file" 2>/dev/null || continue
    new=$(git show "$SRC:$file" | python3 -c "import sys; s=sys.stdin.read(); f,t=sys.argv[1:3]; sys.stdout.write(s.replace(f,t))" "$from" "$to" | git hash-object -w --stdin)
    git update-index --cacheinfo "100644,$new,$file"
  done < "$REDACTIONS"
else
  echo "note: no $REDACTIONS file; no redactions applied" >&2
fi

# checks
bad=$(git ls-files | grep -E "$FORBID_RE" || true); [ -z "$bad" ] || { echo "forbidden paths in snapshot:"; echo "$bad"; exit 1; }
USERS="/Use""rs/"                                      # generic check for local absolute paths
args=(-e "$USERS"); for t in "${SENSITIVE[@]:-}"; do [ -n "$t" ] && args+=(-e "$t"); done
if git grep --cached -lIF "${args[@]}" -- . ':!uv.lock' | grep -q .; then echo "sensitive string or local path found in:"; git grep --cached -lIF "${args[@]}" -- . ':!uv.lock'; exit 1; fi
big=$(git ls-files -s | while read -r mode sha stage path; do s=$(git cat-file -s "$sha"); [ "$s" -gt "$MAX_BYTES" ] && echo "$s $path"; done || true)
[ -z "$big" ] || { echo "files over 20 MB:"; echo "$big"; exit 1; }

tree=$(git write-tree); n=$(git ls-files | wc -l | tr -d ' ')
total=$(git ls-files -s | awk '{print $2}' | git cat-file --batch-check='%(objectsize)' | awk '{s+=$1} END {printf "%.1f", s/1048576}')
parent=$(git rev-parse -q --verify "refs/heads/$BR" || true)
echo "snapshot of $SRC: $n files, $total MB"
if [ -n "$parent" ] && [ "$(git rev-parse "$parent^{tree}")" = "$tree" ]; then echo "no changes since $(git rev-parse --short "$parent")"; exit 0; fi
if [ "$DRY" = 1 ]; then echo "dry run: would commit on $BR"; [ -n "$parent" ] && git diff-tree --stat -r "$parent" "$tree" | tail -5; exit 0; fi
commit=$(printf '%s\n' "$MSG" | git commit-tree "$tree" ${parent:+-p "$parent"})
git update-ref "refs/heads/$BR" "$commit"
echo "$BR -> $(git rev-parse --short "$commit"). Review with: git show --stat $BR; publish with: git push origin $BR:main"
