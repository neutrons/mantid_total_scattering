# Get python version
FROM python:2.7-slim

# Install mantid image
FROM mantidproject/mantid:nightly_ubuntu16.04

# Copy git content from current branch
COPY . /root/mantid_total_scattering

# Add Mantid to python path
ENV MANTIDPATH         /opt/mantidnightly/bin
ENV TSREPO             /root/mantid_total_scattering
ENV PYTHONPATH         ${MANTIDPATH}:${TSREPO}:${PYTHONPATH}

# Install python dependencies
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install python-pip -y && \
    pip install pytest

# Copy Insturment geometry caches
RUN mkdir -p /root/.mantid/instrument/geometryCache

COPY ./dockerfiles/POLARIS9fbf7121b4274c833043ae8933ec643ff7b9313d.vtp /root/.mantid/instrument/geometryCache/POLARIS9fbf7121b4274c833043ae8933ec643ff7b9313d.vtp

# Run unit tests
CMD ["pytest", "/root/mantid_total_scattering/"]
