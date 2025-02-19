# gcloud artifacts repositories create carriers-repo \
#     --repository-format=docker \
#     --location=us-central1 \
#     --description="Repository for carriers API"

gcloud builds submit --tag europe-west4-docker.pkg.dev/gp2-release-terra/genotools/carriers-api

gcloud run deploy carriers-api \
    --image europe-west4-docker.pkg.dev/gp2-release-terra/genotools/carriers-api \
    --platform managed \
    --region europe-west4 \
    --allow-unauthenticated \
    --memory 16Gi \
    --cpu 4 \
    --timeout 3600 \
    --execution-environment gen2
