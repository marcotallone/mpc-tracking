import re

# Open the source file
with open('main.tex', 'r') as source_file:
    source_content = source_file.read()

# Find all \begin{info} ... \end{info} blocks
info_blocks = re.findall(r'\\begin{info}.*?\\end{info}', source_content, re.DOTALL)

# Open the destination file
with open('info_summary.tex', 'w') as destination_file:
    for block in info_blocks:
        # Write each block to the destination file
        destination_file.write(block + '\n\n')