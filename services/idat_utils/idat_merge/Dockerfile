FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
COPY . .

RUN chmod +x /app/exec/plink2 && \
    mv /app/exec/plink2 /usr/local/bin/

RUN chmod +x /app/exec/plink1.9 && \
    mv /app/exec/plink1.9 /usr/local/bin/

RUN chmod +x /app/exec/iaap-cli && \
    mv /app/exec/iaap-cli /usr/local/bin/

RUN pip install --no-cache-dir -r requirements.txt

CMD ["batch_idat.py"]
