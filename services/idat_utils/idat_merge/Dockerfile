FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
COPY . ./

RUN pip install --no-cache-dir -r requirements.txt

RUN chmod +x /app/exec/plink2_module.sh && \
    mv /app/exec/plink2_module.sh /usr/local/bin/

RUN chmod +x /app/exec/plink_module.sh && \
    mv /app/exec/plink_module.sh /usr/local/bin/

RUN chmod +x /app/exec/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli && \
    mv /app/exec/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli /usr/local/bin/

RUN apt-get update && apt-get -y install gcsfuse
RUN gcsfuse --dir-mode 777 --file-mode 777 --implicit-dirs gp2_genotools_data ~/gp2_genotools_data

CMD ["run_batch_idat.py", "key_path", "study_name"]
