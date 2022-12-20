
from ._version import get_versions
__version__ = get_versions()['version']
__date__ = get_versions()['date']
__commit__ = get_versions()['full-revisionid']
del get_versions
