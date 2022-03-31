import json
import os
import sys

import ssl

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context

try:
    import urllib
    fetcher = urllib.urlretrieve
except AttributeError:
    import urllib.request
    fetcher = urllib.request.urlretrieve


this_module_path = sys.modules[__name__].__file__


def fetch_remote_files():
    r"""Fetch remote files before any test runs"""
    data_dir = os.path.dirname(this_module_path)
    remote_info_file = os.path.join(data_dir, 'remote_files.json')
    remote_info = json.load(open(remote_info_file, 'r'))

    remote_address = remote_info['address']
    remote_files = remote_info['files']

    for basename, md5 in remote_files.items():
        file_path = os.path.join(data_dir, basename)
        if not os.path.isfile(file_path):
            print('Fetching data file ' + basename)
            fetcher(os.path.join(remote_address, md5), file_path)


if __name__ == '__main__':
    fetch_remote_files()
