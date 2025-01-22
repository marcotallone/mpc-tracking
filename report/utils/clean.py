import os

# Specify the extensions of the files you want to delete
extensions = [
    ".aux", ".log", ".out", ".toc", ".lof", ".lot",
    ".fls", ".blg", ".bbl", ".bcf", ".fdb_latexmk", 
    ".synctex.gz",".listing","-blx.bib", ".run.xml",
    ".nav", ".snm", ".vrb", ".synctex.gz", ".synctex.gz(busy)"
]

# Walk through the previous directory and its subdirectories
for dirpath, dirnames, filenames in os.walk('.'):
    # For each file in the current directory
    for filename in filenames:
        # If the file has one of the specified extensions
        if any(filename.endswith(ext) for ext in extensions):
            # Construct the full file path
            file_path = os.path.join(dirpath, filename)
            # Delete the file
            os.remove(file_path)
            print(f"Deleted {file_path}")