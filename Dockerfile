FROM quay.io/biocontainers/intarna:3.3.1--pl5321h7ff8a90_1
COPY ./app /app
RUN chmod -R 755 /app && mv /app/tRFtarget /usr/local/bin && mkdir /data
RUN pip install Bio pandas