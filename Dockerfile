FROM quay.io/biocontainers/intarna:3.1.3--h176a8bc_0
COPY ./app /app
RUN chmod -R 755 /app && mv /app/tRFtarget /usr/local/bin && mkdir /data
RUN pip install Bio pandas