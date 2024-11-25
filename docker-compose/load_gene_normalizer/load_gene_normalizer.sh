#!/bin/sh

psql $POSTGRES_URL -c "CREATE DATABASE gene_normalizer;"
wget -q -O gene_normalizer_dump.tar.gz https://vicc-normalizers.s3.us-east-2.amazonaws.com/gene_normalization/postgresql/gene_norm_20240529154335.sql.tar.gz
tar -xzvf gene_normalizer_dump.tar.gz
psql $POSTGRES_URL -f gene_norm_20240529154335.sql
echo "Loaded Gene Normalizer data."
