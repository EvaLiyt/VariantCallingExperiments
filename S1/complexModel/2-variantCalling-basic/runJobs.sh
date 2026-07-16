#!/bin/bash

for dir in D*V*; do
  [ -d "$dir" ] || continue   # skip if not a directory
  echo "Entering $dir"

  cd "$dir" || exit

  for file in jobs/run_r*sh; do
    [ -f "$file" ] || continue
    chmod +x "$file"
    ./"$file"
    echo "done $dir/$file"
  done

  cd - > /dev/null
done
