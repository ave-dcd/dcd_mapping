FROM python:3.12-slim

# Install required packages
RUN apt-get update && apt-get install -y wget postgresql-client && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /usr/src/app

# Copy the script to the container
COPY load_gene_normalizer.sh .

# Make the script executable
RUN chmod +x load_gene_normalizer.sh
