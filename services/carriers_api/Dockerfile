FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
COPY . .

RUN chmod +x /app/exec/plink2 && \
    mv /app/exec/plink2 /usr/local/bin/

RUN pip install --no-cache-dir -r requirements.txt

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]