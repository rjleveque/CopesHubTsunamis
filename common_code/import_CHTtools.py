"""
Code to include to import CHTtools module, which now exists in a
different repository named CHTuser.
"""


from clawpack.clawutil.util import fullpath_import

try:
    CHTuser = os.environ['CHTuser']
except:
    raise Exception("*** Must first set CHTuser environment variable")

CHTtools = fullpath_import(f'{CHTuser}/src/CHTuser/CHTtools.py')
