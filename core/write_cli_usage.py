"""
Write a markdown documentation file for command line scripts.
"""
import importlib
from glob import glob
import os.path
import sys


FOLDER_DOC = os.path.join('..', '..', 'TatooineMesher.wiki')

FOLDER_SCRIPTS = os.path.join('..')

URL_WIKI = 'https://github.com/CNR-Engineering/TatooineMesher/wiki'


sys.path.append(FOLDER_SCRIPTS)  # dirty method to import modules easily


class CommandLineScript:
    def __init__(self, path):
        self.path = path
        basename = os.path.basename(self.path)
        self.name = os.path.splitext(basename)[0]

    def help_msg(self):
        """Returns help message with description and usage"""
        mod = importlib.import_module('%s' % self.name)
        return getattr(mod, 'parser').format_help()


# Build sorted list of CLI scripts
with open(os.path.join(FOLDER_DOC, '_Sidebar.md'), 'w') as out_sidebar:
    for file_path in sorted(glob(os.path.join(FOLDER_SCRIPTS, '*.py'))):
        script_name = os.path.splitext(os.path.basename(file_path))[0]

        if not script_name.startswith('_'):
            out_sidebar.write('* [[%s]]\n' % script_name)

            script = CommandLineScript(file_path)

            # Write a markdown file (to be integrated within github wiki)
            with open(os.path.join(FOLDER_DOC, script_name + '.md'), 'w') as fileout:
                fileout.write('```\n')
                fileout.write(script.help_msg())
                fileout.write('```\n')
                fileout.write('\n')
