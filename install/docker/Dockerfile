# Use the official Python image from the Docker Hub
FROM python:3.6.15-slim

# Set the working directory
WORKDIR /app

# Copy the requirements file into the container
COPY requirements.txt .

# Install the dependencies
RUN apt-get update && \
        apt-get install -y --no-install-recommends \
        curl \
        gcc \
        build-essential \
        libffi-dev \
        libssl-dev \
        git \
        samtools \
        zlib1g-dev \
        libpng-dev \
        && \
        pip install --no-cache-dir -r requirements.txt && \
        rm -rf /var/lib/apt/lists/*

# Install TRF
RUN git clone https://github.com/Benson-Genomics-Lab/TRF.git /opt/TRF
# Set up build directory
WORKDIR /opt/TRF/build
# Configure and build TRF
RUN ../configure && make
# Make it system-wide available
RUN cp src/trf /usr/local/bin/trf

# Install Otter
RUN git clone --branch development --recursive https://github.com/holstegelab/otter.git /opt/otter2
# Set up build directory
WORKDIR /opt/otter2/include/WFA2-lib
RUN make clean setup lib_wfa
WORKDIR /opt/otter2
RUN mkdir build
RUN make
# Make it system-wide available
RUN cp build/otter /usr/local/bin/otter

# Clone the TREAT GitHub repository
RUN git clone https://github.com/holstegelab/treat.git /app/treat

# Change working directory to the cloned repository
WORKDIR /app/treat

# Clean image by removing unnecessary programs
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libssl-dev \
        libffi-dev \
        && apt-get remove -y build-essential \
        && apt-get autoremove -y \
        && rm -rf /var/lib/apt/lists/*

# Copy the rest of your application files
# COPY . .

# Define the default command to run your application
ENTRYPOINT ["python", "bin/TREAT.py"]
#ENTRYPOINT ["bash"]
