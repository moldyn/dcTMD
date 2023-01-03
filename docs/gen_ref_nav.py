"""Generate the code reference pages and navigation.

This file is taken and adopted from:
https://github.com/mkdocstrings/mkdocstrings

"""

from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()

for path in sorted(Path('src').glob('**/*.py')):
    module_path = path.relative_to('src').with_suffix('')
    doc_path = path.relative_to('src', 'dcTMD').with_suffix('.md')
    full_doc_path = Path('reference', doc_path)

    parts = list(module_path.parts)
    if parts[-1] == '__init__':
        parts = parts[:-1]
        doc_path = doc_path.with_name('index.md')
        full_doc_path = full_doc_path.with_name('index.md')
    elif parts[-1].startswith('_'):
        continue

    nav[parts] = doc_path

    module_string = '.'.join(parts)
    with mkdocs_gen_files.open(full_doc_path, 'w') as fd:
        fd.writelines('::: {m}'.format(m=module_string))
    mkdocs_gen_files.set_edit_path(full_doc_path, Path('../') / path)

with mkdocs_gen_files.open('reference/index.md', 'w') as nav_file:
    nav_file.writelines(nav.build_literate_nav())
