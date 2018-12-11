# Get python version
FROM python:3.6-slim

# Install mantid image
FROM mantidproject/mantid:nightly_ubuntu16.04

# Pass in git path
#ARG BRANCH=master

# TEST
#RUN echo "Building git branch: $BRANCH"

# Copy git content from current branch
COPY . /root/mantid_total_scattering

# Install pytest
RUN pip3 install pytest

# Copy Insturment geometry caches
RUN mkdir /root/.mantid && \
    mkdir /root/.mantid/instrument/ && \
    mkdir /root/.mantid/instrument/geometryCache

COPY ./dockerfiles/POLARIS9fbf7121b4274c833043ae8933ec643ff7b9313d.vtp /root/.mantid/instrument/geometryCache/POLARIS9fbf7121b4274c833043ae8933ec643ff7b9313d.vtp

# Run unit tests
CMD ["python3", "-m", "pytest", "/root/mantid_total_scattering/"]
