FROM marshallmcdonnell/mantid:latest

# Install mantid_total_scattering
RUN cd /root && \
    git clone https://github.com/marshallmcdonnell/mantid_total_scattering.git

