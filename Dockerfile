# Get python version
FROM python:3.6-slim

# Install mantid image
FROM marshallmcdonnell/mantid:latest

# Install mantid_total_scattering
RUN cd /root && \
    git clone https://github.com/marshallmcdonnell/mantid_total_scattering.git


# Install pytest
RUN pip3 install pytest

RUN mkdir /root/.mantid && \
    mkdir /root/.mantid/instrument/ && \
    mkdir /root/.mantid/instrument/geometryCache

# Copy Instrument geomtry cache
COPY ./dockerfiles/POLARIS9fbf7121b4274c833043ae8933ec643ff7b9313d.vtp /root/.mantid/instrument/geometryCache/POLARIS9fbf7121b4274c833043ae8933ec643ff7b9313d.vtp

# Run unit tests
CMD ["python3", "-m", "pytest", "/root/mantid_total_scattering/file_handling/test/"]
