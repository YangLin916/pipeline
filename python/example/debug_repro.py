import sys
import pathlib

# Simulate the environment where 'pipeline' is installed (in sys.path)
# We need to construct the expected path based on where this script is running
# Script location: /Users/atlab_yang/Documents/atlab/pipeline/python/example/debug_repro.py
# Expected pkg_dir: /Users/atlab_yang/Documents/atlab/pipeline/python

current_dir = pathlib.Path(__file__).parent.resolve()
pkg_dir_expected = current_dir.parent.resolve()
print(f"Expected package dir: {pkg_dir_expected}")

# Add to sys.path to simulate "installed"
# Note: sys.path usually contains absolute paths
if str(pkg_dir_expected) not in sys.path:
    sys.path.append(str(pkg_dir_expected))

# Original Code Logic
candidates = [pathlib.Path.cwd(), *pathlib.Path.cwd().parents]
pkg_dir = None
for c in candidates:
    candidate = c / 'python'
    if (candidate / 'pipeline').exists():
        pkg_dir = candidate
        break

print(f"Found pkg_dir: {pkg_dir}")
if pkg_dir:
    print(f"Is pkg_dir in sys.path? {str(pkg_dir) in sys.path}")
else:
    print("pkg_dir not found")

if pkg_dir and str(pkg_dir) not in sys.path:
    sys.path.insert(0, str(pkg_dir))
    print('Added to sys.path:', pkg_dir)
else:
    print('Could not find pipeline package; set PYTHONPATH or pip install -e .')
