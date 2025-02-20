# gcloud artifacts repositories create carriers-repo \
#     --repository-format=docker \
#     --location=us-central1 \
#     --description="Repository for carriers API"

gcloud builds submit --tag europe-west4-docker.pkg.dev/gp2-release-terra/genotools/carriers-api

python scripts/generate_key.py | \
gcloud secrets create carriers-api-key \
    --data-file=- \
    --replication-policy="automatic"
    
gcloud secrets add-iam-policy-binding carriers-api-key \
    --member="serviceAccount:776926281950-compute@developer.gserviceaccount.com" \
    --role="roles/secretmanager.secretAccessor"


gcloud run deploy carriers-api \
    --image europe-west4-docker.pkg.dev/gp2-release-terra/genotools/carriers-api \
    --platform managed \
    --region europe-west4 \
    --allow-unauthenticated \
    --memory 16Gi \
    --cpu 4 \
    --timeout 3600 \
    --execution-environment gen2 \
    --set-secrets "API_KEY=carriers-api-key:latest"

