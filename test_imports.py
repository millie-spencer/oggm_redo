# test_imports.py
try:
    import esmf
    import esmpy
    print("ESMF and ESMPy imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")