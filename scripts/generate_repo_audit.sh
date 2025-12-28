#!/usr/bin/env bash
out=docs/REPO_AUDIT.md
echo "# Repository Genesets Audit" > "$out"
echo "Date: $(date -u +%Y-%m-%d)" >> "$out"
echo >> "$out"
echo "## Inventory" >> "$out"
echo "See docs/file_inventory.txt for full list." >> "$out"
echo >> "$out"
while read -r f; do
  echo "### $f" >> "$out"
  echo "- size: $(stat -c%s "$f") bytes" >> "$out"
  lines=$(wc -l < "$f" 2>/dev/null || echo "?")
  echo "- lines: $lines" >> "$out"
  case "$f" in
    *.csv|*.tsv)
      hdr=$(head -n1 "$f" | tr -d '\r')
      echo "- header: \"$hdr\"" >> "$out"
      # first-column symbols (CSV/TSV): handle comma or tab
      delim=','
      [[ "$f" == *.tsv ]] && delim=$'\t'
      uniqs=$(awk -F"$delim" 'NR>1{if($1!="") print $1}' "$f" | sed '/^$/d' | wc -l)
      dups=$(awk -F"$delim" 'NR>1{if($1!="") print $1}' "$f" | sed '/^$/d' | sort | uniq -d | wc -l)
      echo "- first-column entries: $uniqs (duplicates: $dups)" >> "$out"
      if [ "$dups" -gt 0 ]; then
        echo "  - sample duplicates:" >> "$out"
        awk -F"$delim" 'NR>1{if($1!="") print $1}' "$f" | sed '/^$/d' | sort | uniq -d | head -n5 | sed 's/^/    - /' >> "$out"
      fi
      ;;
    *.gmt)
      echo "- GMT file (gene sets):" >> "$out"
      # count sets and basic gene count stats
      sets=$(wc -l < "$f")
      echo "  - sets: $sets" >> "$out"
      # compute min/max/median genes per set
      awk -F"\t" '{print NF-2}' "$f" | sort -n > /tmp/gene_counts.txt
      if [ -s /tmp/gene_counts.txt ]; then
        min=$(head -n1 /tmp/gene_counts.txt)
        max=$(tail -n1 /tmp/gene_counts.txt)
        median=$(awk '{a[NR]=$1} END{if(NR%2==1) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}' /tmp/gene_counts.txt)
        echo "  - genes per set: min=$min, max=$max, median=$median" >> "$out"
      fi
      ;;
    *)
      echo "- unknown file type" >> "$out"
      ;;
  esac
  echo >> "$out"
done < docs/file_inventory.txt

# Summary recommendations
cat >> "$out" <<'EOF'

## Summary Recommendations
- Consolidate canonical files under `data/genesets/sources/` and `data/genesets/custom/` for project-specific sets.
- Add per-file metadata and checksums (manifest present but expand entries).
- Run downstream tools against `manifest.yaml` to centralize paths.
- Add CI checks for parseability and checksum validation (optional networked tests separate).
EOF

chmod +x scripts/generate_repo_audit.sh
./scripts/generate_repo_audit.sh
