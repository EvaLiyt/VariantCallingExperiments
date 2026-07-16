for script in D*/1-generateScript.py; do
    dir=$(dirname "$script")
    echo "Running $script"
    (cd "$dir" && python3 1-generateScript.py)
done
echo "All done"
