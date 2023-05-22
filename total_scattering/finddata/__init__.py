# finddata version
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
# everything else
from .publish_plot import publish_plot
