import os
import re


def remove_large_insertions(filename):
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                description = line.strip()
            elif i == 1:
                sequence = line.strip()
            elif i == 2:
                structure = line.strip()
    match = re.finditer(r'(\.{200,})', structure)
    if not match:
        return
    new_sequence = sequence
    new_structure = structure
    for span in match:
        new_sequence = new_sequence[:span.start()] + '@' * (span.end() - span.start()) + new_sequence[span.end():]
        new_structure = new_structure[:span.start()] + '@' * (span.end() - span.start()) + new_structure[span.end():]
    new_sequence = re.sub(r'@+', 'XXXX', new_sequence)
    new_structure = re.sub(r'@+', '....', new_structure)
    os.system('mv {0} {0}.insertion'.format(filename))
    with open(filename, 'w') as f:
        f.write(description + '\n')
        f.write(new_sequence + '\n')
        f.write(new_structure + '\n')
